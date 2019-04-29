#!/usr/bin/env python
from Bio import SeqIO
import sys

for rec in SeqIO.parse(sys.argv[1], "fastq"):
    if rec.id in sys.argv[2:]:
        SeqIO.write([rec], sys.stdout, "fastq")
