#!/usr/bin/env python

import sys
import csv

import sqlite3
from Bio import SeqIO

con = sqlite3.connect(sys.argv[1])
con.row_factory = sqlite3.Row
cur = con.cursor()

def get_runs(dataset):
    if dataset == 'all':
            cur.execute("select * from runs where include = 'T'")
    else:
        cur.execute("select * from runs where runs.dataset = ? and include = 'T'", (dataset,))
        return cur.fetchall()

runs = get_runs(sys.argv[2])
for row in runs:
        if row['Include'] != 'T': continue
   
#        print "%s\t%s" % (
#            'EM_079517_mut30_2.mutations.txt',
#            '%s_hq_EM_079517_mut30_2_marginAlign.tagged.vcf' % (row['Batch'])
#        )
#        for refnum in [1,2,3,4,5]:
        for refnum in [2]:
            print("%s_hq\tnp-new\t%s\t%s" % (
                row['Batch'],
                '../refs/EM_079517_mut30_%s.mutations.txt' % (refnum),
                'np_EM_079517_mut30_%s_%s_hq.vcf' % (refnum, row['Batch'])
            ))
#    print "%s\tnp-new-filter\t%s\t%s" % (
#            row['Batch'],
#            'EM_079517_mut30_2.mutations.txt',
#            '%s_hq_EM_079517_mut30_2_np_primer.filtered.vcf' % (row['Batch'])
#        )
        print("%s\tnp-new-filter075-30\t%s\t%s" % (
            row['Batch'],
            '../refs/EM_079517_mut30_2.mutations.txt',
            '%s_hq_EM_079517_mut30_2_np_primer.filtered075_30.vcf' % (row['Batch'])
        ))
        print("%s\tnp-new-filter_qual200-50\t%s\t%s" % (
            row['Batch'],
            '../refs/EM_079517_mut30_2.mutations.txt',
            '%s_hq_EM_079517_mut30_2_np_primer.filtered_qual200.vcf' % (row['Batch'])
        ))


