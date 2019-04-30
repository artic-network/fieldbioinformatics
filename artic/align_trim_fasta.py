#!/usr/bin/env python

#Written by Nick Loman
#Part of the ZiBRA pipeline (zibraproject.org)

import argparse
import pysam
import sys
from .align_trim import find_primer
from .vcftagprimersites import read_bed_file

def find_query_pos(alignment, reference_pos):
    nearest = -1
    for qry, ref in alignment.get_aligned_pairs():
        if qry == None or ref == None: continue

        nearest = qry
        if ref >= reference_pos:
            return nearest
    return nearest

def go(args):
    bed = read_bed_file(args.bedfile)

    infile = pysam.AlignmentFile(args.alignment, "rb")
    for s in infile:
        #print s.get_aligned_pairs()
        #print ">%s\n%s" % (s.query_name, s.query_alignment_sequence)

        p1 = find_primer(bed, s.reference_start, '+')
        p2 = find_primer(bed, s.reference_end, '-')

        primer_start = p1[2]['start']
        # start is the 5' 
        primer_end = p2[2]['start']

        query_align_start = find_query_pos(s, primer_start)
        query_align_end = find_query_pos(s, primer_end)

        print("%s\t%s\t%s\t%s" % (primer_start, primer_end, primer_end - primer_start, s.query_length), sys.stderr)

        startpos = max(0, query_align_start - 40)
        endpos = min(query_align_end+40, s.query_length)

        print(">%s\n%s" % (s.query_name, s.query_sequence[startpos:endpos]))
        #query_align_end + 30])

parser = argparse.ArgumentParser(description='Trim alignments from an amplicon scheme.')
parser.add_argument('alignment', help='BAM alignment')
parser.add_argument('bedfile', help='BED file')

args = parser.parse_args()
go(args)

