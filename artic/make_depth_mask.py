#!/usr/bin/env python
from Bio import SeqIO
import itertools
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
    list
        Index is the reference position, value is the corresponding coverage depth
    dict
        Key is readgroup, value is a list of coverage depths where index is the reference position
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

    # create a depth vector to hold the depths at each reference position
    depths = [0] * bamFile.get_reference_length(refName)

    # create the dict to hold the depths for each readgroup
    rgDepths = {}

    # get the read groups and init the depth vectors
    for rg in bamFile.header["RG"]:
        if rg["ID"] == "unmatched":
            continue
        rgDepths[rg["ID"]] = [0] * bamFile.get_reference_length(refName)

    # flag to state if BAM file has low readgroup coverage
    lowRGcov = False

    # vector to keep track of low readgroup coverage regions
    lowRGvec = []

    # generate the pileup
    for pileupcolumn in bamFile.pileup(
        refName, max_depth=10000, truncate=False, min_base_quality=0
    ):

        # process the pileup column
        for pileupread in pileupcolumn.pileups:

            # get the read group for this pileup read and check it's in the BAM header
            rg = pileupread.alignment.get_tag("RG")
            assert rg in rgDepths, "alignment readgroup not in BAM header: %s" % rg

            # process the pileup read
            if pileupread.is_refskip:
                continue
            if pileupread.is_del:
                if not ignoreDeletions:
                    depths[pileupcolumn.pos] += 1
                    rgDepths[rg][pileupcolumn.pos] += 1
            elif not pileupread.is_del:
                depths[pileupcolumn.pos] += 1
                rgDepths[rg][pileupcolumn.pos] += 1
            else:
                raise Exception("unhandled pileup read encountered")

        # if final depth for pileup column < minDepth, report 0 and update the mask_vector
        if depths[pileupcolumn.pos] < minDepth:
            depths[pileupcolumn.pos] = 0

        # if pileupcolumn depth is okay, check that at least one of the readgroups > minDepth
        else:
            rgCovCheck = 0
            for rg in rgDepths:
                if rgDepths[rg][pileupcolumn.pos] >= minDepth:
                    rgCovCheck += 1
            if rgCovCheck == 0:

                # mask the region if it has low coverage in all readgroups
                depths[pileupcolumn.pos] = 0

                # also record the region so that we can report it if requested
                lowRGcov = True
                lowRGvec.append(pileupcolumn.pos)

    # if requested, warn if there are regions with low readgroup coverage that pass the combined depth threshold
    if warnRGcov and lowRGcov:

        # get the regions and print warning
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
    record = list(SeqIO.parse(args.reference, "fasta"))[0]
    seqID = record.id
    seqLength = len(record.seq)

    # collect the depths from the pileup, replacing any depth<minDepth with 0
    try:
        depths, rgDepths = collect_depths(
            args.bamfile,
            seqID,
            args.depth,
            args.ignore_deletions,
            args.warn_rg_coverage,
        )
    except Exception as e:
        print(e)
        raise SystemExit(1)

    # check the number of positions in the reported depths matches the reference sequence
    if len(depths) != seqLength:
        print("pileup length did not match expected reference sequence length")

    # print the readgroup depths to individual files if requested
    if args.store_rg_depths:

        # rg is the readgroup and rgd is the depths per position for this readgroup
        for rg, rgd in rgDepths.items():
            fh = open(args.outfile + "." + rg + ".depths", "w")
            for pos, depth in enumerate(rgd):
                fh.write("%s\t%s\t%d\t%d\n" % (seqID, rg, pos, depth))
            fh.close()

    # create a mask_vector that records reference positions where depth < minDepth
    mask_vector = []
    for pos, depth in enumerate(depths):
        if depth == 0:
            mask_vector.append(pos)

    # get the intervals from the mask_vector
    intervals = list(intervals_extract(mask_vector))

    # create the mask outfile
    maskfh = open(args.outfile, "w")
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
