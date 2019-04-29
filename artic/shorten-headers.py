#!/usr/bin/env python

from Bio import SeqIO
import sys

for record in SeqIO.parse(sys.stdin, 'fasta'):
    if record.id.startswith('EBOV'):
        record.id = record.id[5:].split('|')[0]
        record.description = record.id
    else:
        record.id = record.id.split('|')[0]
        record.description = record.id
    SeqIO.write(record, sys.stdout, 'fasta')
