#!/usr/bin/env python

import sys
import vcf

MAX_COORD = 18588

def read_truth_set(fn):
    truthset = set()
    for ln in open(fn):
        cols = ln.rstrip().split("\t")
        pos = int(cols[0])
    #    if pos > MAX_COORD: continue

        truthset.add(int(cols[0]))
    return truthset

def read_vcf(fn):
    vcfset = set()
    vcfinfo = {}
    vcf_reader = vcf.Reader(open(fn, 'r'))
    for record in vcf_reader:
        if record.POS > MAX_COORD: continue
        vcfset.add(record.POS)    
        vcfinfo[record.POS] = record.INFO
    return vcfset, vcfinfo

if sys.argv[1].endswith('.vcf'):
    truthset, ignore = read_vcf(sys.argv[1])
else:
    truthset = read_truth_set(sys.argv[1])
for vcffile in sys.argv[2:]:
    vcfset, vcfinfo = read_vcf(vcffile)

    tp = float(len(vcfset & truthset))
    fn = float(len(truthset - vcfset))
    tpr = tp / (tp + fn)

    print("%s\t%s\t%s\t%s\t%s\t%s" % (vcffile, len(vcfset), len(truthset), tp, fn, tpr))
    for sample in vcfset & truthset:
        print(sample, vcfinfo[sample], sys.stderr)
    print("Missing: %s" % (truthset - vcfset), sys.stderr)
    print("Extra:   %s" % (vcfset - truthset), sys.stderr)

