#!/usr/bin/env python

import sys
import csv
from Bio import SeqIO
#from Bio.Seq import Seq
#from Bio.SeqRecord import SeqRecord
#from collections import defaultdict

def main():
    records = SeqIO.to_dict(SeqIO.parse(open(sys.argv[1]), 'fasta'))

    reader = csv.DictReader(sys.stdin, dialect="excel-tab")
    clusters = list(reader)

    groups = set([c['group'] for c in clusters])

    for group in groups:
        print("cluster%s\t%s-cluster%s" % (group, sys.argv[1], group))
        with open('%s-cluster%s' %(sys.argv[1], group), 'w') as fout:
            SeqIO.write([records[i['node']] for i in clusters if i['group'] == group], fout, 'fasta')

if __name__ == '__main__':
        main()
