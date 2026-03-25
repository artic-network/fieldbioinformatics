#!/usr/bin/env python
from Bio import SeqIO
import itertools
import numpy as np
import os
import pysam
import sys


def collect_depths(bamfile, refName, minDepth, ignoreDeletions, warnRGcov):
    """Collect read depth of coverage per reference position in a BAM file.

    Parameters
    ----------
    bamfile : string
        The BAM file that needs processing

    refName : string
        The name of the reference sequence to collect the depths for

    minDepth : int
        The minimum depth to report coverage for (0 will be reported if coverage < minDepth at a given position)

    ignoreDeletions : bool
        If true, positional depth counts will ignore reads with reference deletions

    warnRGcov : bool
        If true, a warning will be issued if the BAM file has pileup regions where coverage for each readgroup < minDepth && combined coverage is > minDepth"

    Returns
    -------
    numpy.ndarray
        Index is the reference position, value is the corresponding coverage depth
    dict
        Key is readgroup, value is a numpy array of coverage depths where index is the reference position
    """
    # check the BAM file exists
    if not os.path.exists(bamfile):
        raise Exception("bamfile doesn't exist (%s)" % bamfile)

    # open the BAM file
    bamFile = pysam.AlignmentFile(bamfile, "rb")

    # get the TID for the reference
    tid = bamFile.get_tid(refName)
    if tid == -1:
        raise Exception("bamfile does not contain specified reference (%s)" % refName)

    ref_len = bamFile.get_reference_length(refName)

    # create a depth vector to hold the depths at each reference position
    depths = np.zeros(ref_len, dtype=np.int32)

    # create the dict to hold the depths for each readgroup
    rgDepths = {rg["ID"]: np.zeros(ref_len, dtype=np.int32) for rg in bamFile.header["RG"]}

    # iterate reads once — each read's RG tag is looked up a single time
    # and numpy slice assignment updates all covered positions in C, avoiding the
    # O(ref_length × coverage) Python loop that pileup() would require.
    # use fetch(refName) when an index is available (faster for multi-contig refs);
    # fall back to a full sequential scan with reference filtering when no index is present.
    try:
        _reads = bamFile.fetch(refName)
    except ValueError:
        _reads = (r for r in bamFile if r.reference_name == refName)

    for read in _reads:
        if read.is_unmapped or read.cigartuples is None:
            continue

        rg = read.get_tag("RG")
        assert rg in rgDepths, "alignment readgroup not in BAM header: %s" % rg

        rg_arr = rgDepths[rg]
        ref_pos = read.reference_start

        for op, length in read.cigartuples:
            if op in (0, 7, 8):
                # M, =, X — match/mismatch: consumes reference, count depth
                depths[ref_pos:ref_pos + length] += 1
                rg_arr[ref_pos:ref_pos + length] += 1
                ref_pos += length
            elif op == 2:
                # D — deletion: consumes reference, count unless ignoreDeletions
                if not ignoreDeletions:
                    depths[ref_pos:ref_pos + length] += 1
                    rg_arr[ref_pos:ref_pos + length] += 1
                ref_pos += length
            elif op == 3:
                # N — reference skip (splice/amplicon gap): consumes reference, do not count
                ref_pos += length
            elif op in (1, 4, 5, 6):
                # I, S, H, P — do not consume reference, skip
                pass

    # vectorized post-processing: apply depth masking and per-RG coverage check

    # stack all per-RG arrays into a 2D matrix for vectorized operations
    rg_matrix = np.stack(list(rgDepths.values()))  # shape: (num_rg, ref_len)

    # mask positions where combined depth < minDepth
    low_combined = depths < minDepth
    depths[low_combined] = 0

    # mask positions where combined depth is adequate but no single RG has >= minDepth
    any_rg_ok = np.any(rg_matrix >= minDepth, axis=0)
    low_rg_mask = ~low_combined & ~any_rg_ok

    lowRGcov = False
    lowRGvec = []
    if np.any(low_rg_mask):
        depths[low_rg_mask] = 0
        lowRGvec = np.where(low_rg_mask)[0].tolist()
        lowRGcov = True

    # if requested, warn if there are regions with low readgroup coverage that pass the combined depth threshold
    if warnRGcov and lowRGcov:
        regions = list(intervals_extract(lowRGvec))
        sys.stderr.write(
            "alignment has unmasked regions where individual readgroup depth < {}: {}\n".format(
                minDepth, bamfile
            )
        )
        for region in regions:
            sys.stderr.write("region: %s\n" % str(region).strip("[]"))

    # close file and return depth vector
    bamFile.close()
    return depths, rgDepths


# from https://www.geeksforgeeks.org/python-make-a-list-of-intervals-with-sequential-numbers/
def intervals_extract(iterable):
    iterable = sorted(set(iterable))
    for _, group in itertools.groupby(enumerate(iterable), lambda t: t[1] - t[0]):
        group = list(group)
        yield [group[0][1], group[-1][1]]


def go(args):

    # open the reference sequence and collect the sequence header and sequence length of the first record
    records = [x for x in SeqIO.parse(args.reference, "fasta")]
    intervals = []

    for record in records:
        seqID = record.id
        seqLength = len(record.seq)

        # collect the depths from the pileup, replacing any depth<minDepth with 0
        depths, rgDepths = collect_depths(
            args.bamfile,
            seqID,
            args.depth,
            args.ignore_deletions,
            args.warn_rg_coverage,
        )

        # check the number of positions in the reported depths matches the reference sequence
        if len(depths) != seqLength:
            print("pileup length did not match expected reference sequence length")

        # print the readgroup depths to individual files if requested
        if args.store_rg_depths:
            # rg is the readgroup and rgd is the depths per position for this readgroup
            for rg, rgd in rgDepths.items():
                fh = open(args.outfile + "." + rg + ".depths", "a")
                for pos, depth in enumerate(rgd):
                    fh.write("%s\t%s\t%d\t%d\n" % (seqID, rg, pos, depth))
                fh.close()

        # create a mask_vector that records reference positions where depth < minDepth
        mask_vector = np.where(depths == 0)[0].tolist()

        # get the intervals from the mask_vector
        intervals = list(intervals_extract(mask_vector))

        # create the mask outfile
        maskfh = open(args.outfile, "a")
        for i in intervals:
            maskfh.write("%s\t%s\t%s\n" % (seqID, i[0] + 1, i[1] + 1))
        maskfh.close()


def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--depth", type=int, default=20)
    parser.add_argument(
        "--warn-rg-coverage",
        action="store_true",
        default=False,
        help="if set, a warning will be issued if the BAM file has pileup regions where coverage for each readgroup < min. depth BUT the combined coverage is > min. depth",
    )
    parser.add_argument(
        "--ignore-deletions",
        action="store_true",
        default=False,
        help="if set, positional depth counts will ignore reads with reference deletions (i.e. evaluates positional depths on ref matches, not read span",
    )
    parser.add_argument(
        "--store-rg-depths",
        action="store_true",
        default=False,
        help="if set, positional depth counts for each readgroup written to file (filename = <outfile>.<readgroup>.depths)",
    )
    parser.add_argument("reference")
    parser.add_argument("bamfile")
    parser.add_argument("outfile")

    args = parser.parse_args()
    go(args)


if __name__ == "__main__":
    main()
