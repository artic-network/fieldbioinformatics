#!/usr/bin/env python

import vcf
import sys
import subprocess
import csv
from collections import defaultdict

def read_bed_file(fn):
    bedfile = []
    with open(fn) as csvfile:
        reader = csv.reader(csvfile, dialect='excel-tab')
        for row in reader:
            bedrow = {}
            bedrow['Primer_ID'] = row[3]
            
            if len(row) >= 6:
                # new style bed
                bedrow['direction'] = row[5]
            elif len(row) == 5:
                # old style without directory
                if 'LEFT' in row[3]:
                    bedrow['direction'] = '+'
                elif 'RIGHT' in row[3]:
                    bedrow['direction'] = '-'
                else:
                    print(("Malformed BED file!", sys.stderr))
                    raise SystemExit

            if bedrow['direction'] == '+':
                bedrow['end'] = int(row[2])
                bedrow['start'] = int(row[1])
            else:
                bedrow['end'] = int(row[1])
                bedrow['start'] = int(row[2])
            bedfile.append(bedrow)
    return bedfile

def overlaps(coords, pos):
    for v in coords:
        if pos >= v['start'] and pos <= v['end']:
            return v
    return False

if __name__ == "__main__":
    if sys.argv[1] not in sets:
        print("Invalid set")
        raise SystemExit

    bedfile = read_bed_file(sys.argv[1])

    vcf_reader = vcf.Reader(filename=sys.argv[2])
    vcf_writer = vcf.Writer(sys.stdout, vcf_reader)
    for record in vcf_reader:
        v = overlaps(bedfile, record.POS)
        if v:
            record.INFO['PRIMER'] = v["Sequence_(5-3')"]

#    PP = list(record.INFO)
#    record.INFO = {}
#    record.INFO['PP'] = PP
#    record.INFO['DEPTH'] = depths[record.CHROM][record.POS]

        vcf_writer.write_record(record)
