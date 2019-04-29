#!/usr/bin/env python

import sys
from Bio import SeqIO

for rec in SeqIO.parse(sys.stdin, "fasta"): print(rec.id, len(rec))
