#!/usr/bin/env python

from Bio import SeqIO
import sys

recs = list(SeqIO.parse(open(sys.argv[1], "r"), "fasta"))
if len (recs) != 1:
    print("FASTA has more than one sequence", file=sys.stderr)
    raise SystemExit

print("%s:%d-%d" % (recs[0].id, 1, len(recs[0])+1))
