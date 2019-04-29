#!/usr/bin/env python

import sqlite3
import re
import sys
import os.path
from Bio import SeqIO
import json
import subprocess
from io import StringIO

""" 
go through the runsamples
look up sample
get fasta
spit it out into a bin
"""

from collections import defaultdict

samples = json.load(open('samples.json'))

def find_sample(sampleid):
    for s in samples['data']:
        if s['sample_id'] == sampleid:
            return s
    return None

goodfh = open("good.fasta", "w")
partialfh = open("partial.fasta", "w")
badfh = open("bad.fasta", "w")

processed = {}

runsamples = json.load(open('runsamples.json'))
for sample in runsamples['data']:
    cols = sample['Run::run_name'].split("\n")
    run_name = cols[0]
    fn = '%s/%s.vcf' % (run_name, sample['sample_id'])
    if not os.path.exists(fn):
        print("No vcf for %s" % (fn,))
        continue

    if fn in processed:
        continue
    processed[fn] = True

    cmd = "margin_cons.py refs/Zika_FP.fasta %s/%s.vcf %s/%s.primertrimmed.sorted.bam" % (run_name, sample['sample_id'], run_name, sample['sample_id'])
    print(cmd)

    p = subprocess.Popen([cmd], shell=True, stdout=subprocess.PIPE)
    out, err = p.communicate()
    del p
    
    rec = list(SeqIO.parse(StringIO(out), "fasta"))[0]
    print(rec)

    metadata = find_sample(sample['sample_id'])
    print(metadata)
    
    """        
{u'pregnancy_week': u'', u'municipality': u'murici', u'patient_sex': u'male', u'host_species': u'human', u'lab_internal_sample_id': u'', u'sample_id': u'ZBRD103', u'minion_barcodes': u'', u'ct': u'29.09', u'lab_id_lacen': u'150101004197', u'collection_date': u'2015-08-20', u'amplicon_concentration_pool_1': u'', u'pregnancy_trimester': u'', u'sample_number': u'103', u'symptoms': u'', u'creation_persistent_id': u'9EDCA6E1F234B3A6E160D5E819D8918D', u'state': u'alagoas', u'extraction_date': u'2016-06-13', u'creation_host_timestamp': u'09/08/2016 21:06:44', u'rt_positive': u'1', u'patient_age': u'25', u'modification_account_name': u'Admin', u'modification_persistent_id': u'9EDCA6E1F234B3A6E160D5E819D8918D', u'lab': u'lacen_maceio', u'onset_date': u'2015-08-18', u'microcephaly': u'', u'sample_type': u'', u'creation_account_name': u'Admin', u'modification_host_timestamp': u'', u'country': u'brazil', u'notes': u'', u'pregnant': u''}
"""

    rec.id = "%s|%s|%s|%s|%s|%s" % (metadata['lab_id_lacen'], metadata['sample_id'], run_name, metadata['municipality'], metadata['state'], metadata['collection_date'])

    if rec.seq.count('N') < 3000:
        SeqIO.write([rec], goodfh, "fasta")
    elif rec.seq.count('N') < 5500:
        SeqIO.write([rec], partialfh, "fasta")
    else:
        SeqIO.write([rec], badfh, "fasta")

    """
    con = sqlite3.connect(sys.argv[1])
    con.row_factory = sqlite3.Row
    cur = con.cursor()

def lookup_sample(sample):
    cur.execute("select * from samples, runs where runs.Batch = ? and runs.sample_fk = samples.rowid", (sample,))
    row = cur.fetchone()
    return row

for rec in SeqIO.parse(sys.stdin, "fasta"):
    m = re.search(r'EM_079517_(.*)_hq', rec.id)
    sample = m.group(1)

    metadata = lookup_sample(sample)

    if metadata['prefecture'] == 'Kambia':
        country = 'SLE'
    else:
        country = 'GUI' 
    loc = '-'.join([metadata['prefecture'], metadata['sousprefecture'], metadata['village'].encode('ascii', 'ignore')])
    rec.id = '|'.join(('EBOV', metadata['LabID'], 'MinION', country, loc, metadata['date_sample_taken'])).replace(' ', '_')
    rec.description = ""
    SeqIO.write([rec], sys.stdout, "fasta")
"""
