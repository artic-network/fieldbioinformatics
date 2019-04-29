#!/usr/bin/env python

from Bio import SeqIO
import sys

START = int(sys.argv[2])
END = int(sys.argv[3])

for rec in SeqIO.parse(open(sys.argv[1]), "fasta"):
    SeqIO.write([rec[START:END]], sys.stdout, "fasta")

