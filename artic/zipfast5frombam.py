#!/usr/bin/env python

# Written by Nick Loman
# zipfast5frombam.py bamfile fastafile zipfile



import pysam
import sys

import zipfile
from Bio import SeqIO

infile = pysam.AlignmentFile(sys.argv[1], "rb")
headers_to_get = [rec.query_name for rec in infile]
fastafile = SeqIO.parse(open(sys.argv[2]), "fasta")
files_to_get = [rec.description.split(' ')[2] for rec in fastafile if rec.id in headers_to_get]
zipf = zipfile.ZipFile(sys.argv[3], 'w', zipfile.ZIP_DEFLATED)
for fn in files_to_get:
    print("Adding %s" % (fn,))
    zipf.write(fn)
zipf.close()
