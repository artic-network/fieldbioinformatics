#!/usr/bin/env python

import os.path
import sys
import vcf

#MAX_COORD = 18588
#MAX_COORD = 20000

def read_truth_set(fn):
    truthset = set()
    for ln in open(fn):
        cols = ln.rstrip().split("\t")
        pos = int(cols[0])
#        if pos > MAX_COORD: continue

        truthset.add(int(cols[0]))
    return truthset

def read_vcf(fn):
    vcfset = {}
    vcf_reader = vcf.Reader(open(fn, 'r'))
    for record in vcf_reader:
#        if record.POS > MAX_COORD: continue
        if 'PP' in record.INFO:
            vcfset[record.POS] = max([float(pp) for pp in record.INFO['PP']])
        else:
            vcfset[record.POS] = float(record.QUAL)
    return vcfset

def filter_set(vcfinfo, threshold):
    vcfset = set()
    for k, v in vcfinfo.items():
        if v >= threshold:
            vcfset.add(k)
    return vcfset

print("sample\ttag\tpp-threshold\tvcf\ttotal_calls\tmutations\tTP\tFP\tFN\tTPR")


for ln in open(sys.argv[1]):
    sample, tag, truthset_fn, vcffile_fn = ln.rstrip().split("\t")
    truthset = read_truth_set(truthset_fn)

    try:
        vcfinfo = read_vcf(vcffile_fn)
    except IOError:
        print("Cannot open %s" % (vcffile_fn), sys.stderr)
        continue

    for threshold in [0.0]:
        vcfset = filter_set(vcfinfo, threshold)

        tp = float(len(vcfset & truthset))
        fn = float(len(truthset - vcfset))
#        print truthset-vcfset
        fp = float(len(vcfset - truthset))
        tpr = tp / (tp + fn)

        print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (sample, tag, threshold, vcffile_fn, len(vcfset), len(truthset), tp, fp, fn, tpr))
#        print "FN: %s" % (truthset - vcfset)
#        print "FP: %s" % (vcfset - truthset)

#for sample in vcfset & truthset:
#        print vcfinfo[sample]
#        print "Missing: %s" % (truthset - vcfset)
#    print "Extra:   %s" % (vcfset - truthset)

