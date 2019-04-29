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
        cur.execute("select ROWID, * from runs where runs.dataset = ?", (dataset,))
        return cur.fetchall()

runs = get_runs(sys.argv[2])

def import_batch(batch, batch_number):
    try:
        fh = open("times/%s.times.txt" % (batch))
    except:
        return

    cols = fh.readline().split("\t")

    if batch_number == 1:
        batchprefix =''
    else:
        batchprefix = '_second_batch'

    cur.execute(
        "UPDATE runs SET start%s = ?, end%s = ?, duration%s = ?, num_reads_pass%s = ?, num_reads_fail%s = ? WHERE ROWID=?" % (batchprefix, batchprefix, batchprefix, batchprefix, batchprefix),
        (cols[2], cols[3], cols[4], cols[0], cols[1], row['ROWID'])
    )


for row in runs:
    import_batch(row['batch'], 1)
    import_batch(row['batch2'], 2)

con.commit()
