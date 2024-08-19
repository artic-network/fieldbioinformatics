#!/usr/bin/env python

#Written by Nick Loman
#Part of the ZiBRA pipeline (zibraproject.org)

import pysam
import sys
from copy import copy
from .align_trim import trim
from collections import defaultdict

def go(args):
    infile = pysam.AlignmentFile("-", "rb")
    outfile = pysam.AlignmentFile("-", "wh", template=infile)
    for s in infile:
        cigar = copy(s.cigartuples)

        if s.is_unmapped:
            print("%s skipped as unmapped" % (s.query_name), file=sys.stderr)
            continue

        if s.is_supplementary:
            print("%s skipped as supplementary" % (s.query_name), file=sys.stderr)
            continue

        if s.is_secondary:
            print("%s skipped as secondary" % (s.query_name), file=sys.stderr)
            continue

        if s.reference_start + args.nbases < s.reference_end:
            trim(args, cigar, s, s.reference_start + args.nbases, 0)
        if s.reference_end - args.nbases > s.reference_start:
            trim(args, cigar, s, s.reference_end - args.nbases, 1)

        outfile.write(s)

def main():
    import argparse

    parser = argparse.ArgumentParser(description='Trim alignments from an amplicon scheme.')
    parser.add_argument('nbases', type=int, help='Number of bases to trim from ends of full length alignments')
    parser.add_argument('--verbose', action='store_true', help='Debug mode')

    args = parser.parse_args()
    go(args)


if __name__ == "__main__":
    main()
