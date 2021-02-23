#!/usr/bin/env python

from copy import copy
from collections import defaultdict, namedtuple
import sys

import numpy as np
import pysam

from artic.vcftagprimersites import read_bed_file

PassRead = namedtuple(
    'PassRead',
    ('qname', 'amplicon', 'coverage', 'left', 'right'))
# consumesReference lookup for if a CIGAR operation consumes the reference sequence
consumesReference = [True, False, True, True, False, False, False, True]
# consumesQuery lookup for if a CIGAR operation consumes the query sequence
consumesQuery = [True, True, False, False, True, False, False, True]


def overlap_size(o1, o2):
    s1, e1 = sorted(o1)
    s2, e2 = sorted(o2)
    return min(e1, e2) - max(s1, s2)


def trim(segment, primer_pos, end, debug):
    """Soft mask an alignment to fit within primer start/end sites.
    Parameters
    ----------
    segment : pysam.AlignedSegment
        The aligned segment to mask
    primer_pos : int
        The position in the reference to soft mask up to (equates to the start/end position of the primer in the reference)
    end : bool
        If True, the segment is being masked from the end (i.e. for the reverse primer)
    debug : bool
        If True, will print soft masking info during trimming
    """
    # get a copy of the cigar tuples to work with
    cigar = copy(segment.cigartuples)

    # get the segment position in the reference (depends on if start or end of the segment is being processed)
    if not end:
        pos = segment.pos
    else:
        pos = segment.reference_end

    # process the CIGAR to determine how much softmasking is required
    eaten = 0
    while 1:

        # chomp CIGAR operations from the start/end of the CIGAR
        try:
            if end:
                flag, length = cigar.pop()
            else:
                flag, length = cigar.pop(0)
            if debug:
                print("Chomped a %s, %s" % (flag, length), file=sys.stderr)
        except IndexError:
            print(
                "Ran out of cigar during soft masking - completely masked read will be ignored", file=sys.stderr)
            break

        # if the CIGAR operation consumes the reference sequence, increment/decrement the position by the CIGAR operation length
        if (consumesReference[flag]):
            if not end:
                pos += length
            else:
                pos -= length

        # if the CIGAR operation consumes the query sequence, increment the number of CIGAR operations eaten by the CIGAR operation length
        if (consumesQuery[flag]):
            eaten += length

        # stop processing the CIGAR if we've gone far enough to mask the primer
        if not end and pos >= primer_pos and flag == 0:
            break
        if end and pos <= primer_pos and flag == 0:
            break

    # calculate how many extra matches are needed in the CIGAR
    extra = abs(pos - primer_pos)
    if debug:
        print("extra %s" % (extra), file=sys.stderr)
    if extra:
        if debug:
            print("Inserted a %s, %s" % (0, extra), file=sys.stderr)
        if end:
            cigar.append((0, extra))
        else:
            cigar.insert(0, (0, extra))
        eaten -= extra

    # softmask the left primer
    if not end:

        # update the position of the leftmost mappinng base
        segment.pos = pos - extra
        if debug:
            print("New pos: %s" % (segment.pos), file=sys.stderr)

        # if proposed softmask leads straight into a deletion, shuffle leftmost mapping base along and ignore the deletion
        if cigar[0][0] == 2:
            if debug:
                print(
                    "softmask created a leading deletion in the CIGAR, shuffling the alignment", file=sys.stderr)
            while 1:
                if cigar[0][0] != 2:
                    break
                _, length = cigar.pop(0)
                segment.pos += length

        # now add the leading softmask
        cigar.insert(0, (4, eaten))

    # softmask the right primer
    else:
        cigar.append((4, eaten))

    # check the new CIGAR and replace the old one
    if cigar[0][1] <= 0 or cigar[-1][1] <= 0:
        raise ("invalid cigar operation created - possibly due to INDEL in primer")
    segment.cigartuples = cigar
    return



def overlap_trim(args):
    """Detect to which amplicon a read is derived and according to trim primers."""

    if args.report:
        reportfh = open(args.report, "w")
        print(
            "QueryName\tReferenceStart\tReferenceEnd\t"
            "PrimerPair\t"
            "Primer1\tPrimer1Start\t"
            "Primer2\tPrimer2Start\t"
            "IsSecondary\tIsSupplementary\t"
            "Start\tEnd\tCorrectlyPaired", file=reportfh)

    # open the primer scheme and get the pools
    bed = read_bed_file(args.bedfile)
    pools = set(row['PoolName'] for row in bed)
    pools.add('unmatched')
    primer_pairs = defaultdict(dict)
    for b in bed:
       scheme, pair, side = b['Primer_ID'].split('_')
       primer_pairs[pair][side] = b
    # this data structure is more useful for searching...
    amplicons = np.fromiter((
        (k, v['LEFT']['PoolName'],
            v['LEFT']['end'], v['RIGHT']['start'],  # just insert
            v['LEFT']['start'], v['RIGHT']['end'],  # contains primers
            v['LEFT']['Primer_ID'], v['RIGHT']['Primer_ID'])
            for k, v in primer_pairs.items()),
        dtype=[
            ('name', int), ('pool', int),
            ('insert_start', int), ('insert_end', int),
            ('start', int), ('end', int),
            ('left_primer', 'U20'), ('right_primer', 'U20')])

    # iterate over the alignment segments in the input SAM file
    passing_reads = defaultdict(list)  # by (amplicon, is_reverse)
    with pysam.AlignmentFile(args.bamfile, "rb") as bam:
        bam_header = bam.header.copy().to_dict()
        if not args.no_read_groups:
            bam_header['RG'] = []
            for pool in pools:
                read_group = {}
                read_group['ID'] = pool
                bam_header['RG'].append(read_group)

        for segment in bam:
            # filter out unmapped and supplementary alignment segments
            if segment.is_unmapped:
                print("%s skipped as unmapped" %
                      (segment.query_name), file=sys.stderr)
                continue
            if segment.is_supplementary:
                print("%s skipped as supplementary" %
                      (segment.query_name), file=sys.stderr)
                continue

            # determine the amplicon by largest overlap (and find second best)
            overlaps = np.zeros(len(amplicons), dtype=int)
            for i, amp in enumerate(amplicons):
                overlaps[i] = segment.get_overlap(amp['start'], amp['end'])
            best, second = np.argpartition(overlaps, -2)[:-3:-1]
            if overlaps[best] < args.min_overlap:
                print("%s skipped as no good overlap" % segment.qname, file=sys.stderr)
                continue

            # don't take reads if the second best overlap is a large proportion
            # of the mutual overlap of the amplicons
            best_amplicon = amplicons[best]
            second_amplicon = amplicons[second]
            mutual = overlap_size(
                *((amp['start'], amp['end'])
                for amp in (best_amplicon, second_amplicon)))
            if overlaps[second] > args.max_mutual_overlap * mutual:
                print("%s skipped as large secondary overlap" % segment.qname, file=sys.stderr)
                continue
            overlap = overlaps[best]
            amplicon = amplicons[best]

            # check whether the alignment extends to the primer at each end
            extends = (
                segment.reference_start <= amplicon['insert_start'],
                segment.reference_end >= amplicon['insert_end'])

            # if both primers, we call that "correctly paired"
            if args.enforce_amplicon_span and not all(extends):
                print("%s skipped as does not span: %s" %
                    (segment.query_name, extends), file=sys.stderr)
                continue
            passing_reads[amplicon['name'], segment.is_reverse].append(
                PassRead(segment.qname, amplicon['name'], overlap, *extends))
    # end first pass

    # filter alignments
    print("Reads before filtering: {}".format(sum(len(x) for x in passing_reads.values())), file=sys.stderr)
    chosen_reads = list()
    if args.normalise:
        for (amp, is_reverse), reads in passing_reads.items():
            reads = sorted(reads, key=lambda x: x.coverage, reverse=True)
            chosen_reads.extend(reads[0:args.normalise])
    else:
        for reads in passing_reads.values():
            chosen_reads.extend(reads)
    print("Reads after filtering: {}".format(len(chosen_reads)), file=sys.stderr)
    chosen_reads = {r.qname:r for r in chosen_reads}

    with \
            pysam.AlignmentFile(args.bamfile, "rb") as bam, \
            pysam.AlignmentFile("-", "wh", header=bam_header) as outfile:
        for segment in bam:
            wanted = segment.qname in chosen_reads
            if not wanted or segment.is_unmapped or segment.is_supplementary:
                continue
            chosen = chosen_reads[segment.qname]
            amplicon = amplicons[amplicons['name'] == chosen.amplicon]
            if not args.no_read_groups:
                segment.set_tag('RG', str(amplicon['pool'][0]))

            if len(amplicon) > 1:
                raise IndexError("Found more than one amplicon matching: {}".format(chosen))
            else:
                amplicon = amplicon[0]
            if args.start:
                trim_start, trim_end = amplicon['start'], amplicon['end']
            else:
                trim_start, trim_end = amplicon['insert_start'], amplicon['insert_end']

            # softmask the alignment if left primer start/end inside alignment
            if segment.reference_start < trim_start:
                try:
                    trim(segment, trim_start, False, args.verbose)
                    if args.verbose:
                        print("ref start %s < primer_position %s" %
                            (segment.reference_start, trim_start),
                            file=sys.stderr)
                except Exception as e:
                    print(
                        "problem soft masking left primer in {} (error: {}), skipping".format(
                            segment.query_name, e),
                        file=sys.stderr)
                    continue

            # softmask the alignment if right primer start/end inside alignment
            if segment.reference_end > trim_end:
                try:
                    trim(segment, trim_end, True, args.verbose)
                    if args.verbose:
                        print("ref end %s > primer_position %s" %
                            (segment.reference_start, p2_position), file=sys.stderr)
                except Exception as e:
                    print(
                        "problem soft masking right primer in {} (error: {}), skipping".format(
                            segment.query_name, e),
                        file=sys.stderr)
                    continue

            # check the the alignment still contains bases matching the reference
            if 'M' not in segment.cigarstring:
                print("%s dropped as does not match reference post masking" %
                      (segment.query_name), file=sys.stderr)
                continue

            # current alignment segment has passed filters, send it to the outfile
            # update the report with this alignment segment + primer details
            if args.report or args.verbose:
                #"QueryName\tReferenceStart\tReferenceEnd\t"
                #"PrimerPair\t"
                #"Primer1\tPrimer1Start\t"
                #"Primer2\tPrimer2Start\t"
                #"IsSecondary\tIsSupplementary\t"
                #"Start\tEnd\tCorrectlyPaired", file=reportfh)
                report = '\t'.join(
                    str(x) for x in (
                        segment.query_name, segment.reference_start, segment.reference_end,
                        '{}_{}'.format(amplicon['left_primer'], amplicon['right_primer']),
                        amplicon['left_primer'], amplicon['start'],
                        amplicon['right_primer'], amplicon['insert_end'],
                        segment.is_secondary, segment.is_supplementary,
                        amplicon['start'], amplicon['end'], int(all(extends))))
                if args.report:
                    print(report, file=reportfh)
                if args.verbose:
                    print(report, file=sys.stderr)
            outfile.write(segment)


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description='Trim alignments from an amplicon scheme.')
    parser.add_argument('bamfile',
                        help='alignment file.')
    parser.add_argument('bedfile',
                        help='BED file containing the amplicon scheme.')
    parser.add_argument('--min-overlap', type=int, default=0, dest='min_overlap',
                        help='Shortest allowed overlap to amplicon region')
    parser.add_argument('--max-mutual-overlap', type=int, default=2, dest='max_mutual_overlap',
                        help='Maximum permissable overlap to second amplicon region as proportion of amplicon mutual overlap')
    parser.add_argument('--normalise', type=int,
                        help='Subsample to n coverage per strand')
    parser.add_argument('--report', type=str,
                        help='Output report to file')
    parser.add_argument('--start', action='store_true',
                        help='Trim to start of primers instead of ends')
    parser.add_argument('--no-read-groups', dest='no_read_groups', action='store_true',
                        help='Do not divide reads into groups in SAM output')
    parser.add_argument('--verbose', action='store_true',
                        help='Debug mode')
    parser.add_argument('--enforce-amplicon-span', action='store_true', dest='enforce_amplicon_span',
                        help='Discard reads that do not cover the entire amplicon')

    args = parser.parse_args()
    overlap_trim(args)


if __name__ == "__main__":
    main()
