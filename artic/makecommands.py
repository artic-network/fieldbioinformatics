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
        cur.execute("select * from runs")
    else:
        cur.execute("select * from runs where runs.dataset = ?", (dataset,))
        return cur.fetchall()

refs = ['EM_079517', 'Makona-G3729']
refs = ['EM_079517', 'EM_079517_mut30_1', 'EM_079517_mut30_2', 'EM_079517_mut30_3', 'EM_079517_mut30_4', 'EM_079517_mut30_5']
refs = ['EM_079517', 'EM_079517_mut30_2']

runs = get_runs(sys.argv[2])

for row in runs:
    for ref in refs:    
        batch2 = row['batch2'] if row['batch2'] else 'na'
        if len(sys.argv) > 3 and sys.argv[3] == 'consensus':
            print("consensus.sh ")
        elif len(sys.argv) > 3:
            print(sys.argv[3])
        else:
            print("align.sh ")
        if len(sys.argv) > 4:
            print("%s %s %s %s_hq %s hq %s" % \
                  (ref, row['batch'], row['batch'], row['batch'], batch2, " ".join(sys.argv[4:])))
        else:
            print("%s %s %s %s_hq %s hq" % \
                  (ref, row['batch'], row['batch'], row['batch'], batch2))
