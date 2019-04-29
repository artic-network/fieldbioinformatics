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
        cur.execute("select ROWID, * from runs")
    else:
        cur.execute("select ROWID, * from runs where runs.dataset = ? AND include = 'T'", (dataset,))
    return cur.fetchall()

runs = get_runs(sys.argv[2])
print("sample    aln    query    read_type    read_len    align_len    unalign_len    matches    mismatches    insertions    deletions    tot_errors")

for row in runs:
    fh = open("EM_079517_%s_hq_marginalign.idystats.txt" % (row['batch']))
    fh.readline()
    for ln in fh:
        print("%s\tma\t%s" % (row['batch'], ln), end=' ') 
    fh.close()

#    fh = open("EM_079517_%s_hq_bwa.idystats.txt" % (row['batch']))
#    fh.readline()
#    for ln in fh:
#        print "%s\tbwa\t%s" % (row['batch'], ln), 
#    fh.close()

