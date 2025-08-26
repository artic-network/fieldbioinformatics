#!/usr/bin/env python

from Bio import SeqIO
from cyvcf2 import VCF
import argparse

from primalbedtools.scheme import Scheme
from primalbedtools.bedfiles import merge_primers


def go(args):
    seqs = dict([(rec.id, rec) for rec in SeqIO.parse(open(args.reference), "fasta")])
    cons = {}
    for k in seqs.keys():
        cons[k] = list(seqs[k].seq)

    scheme = Scheme.from_file(args.maskfile)

    scheme.bedlines = merge_primers(scheme.bedlines)

    for region in scheme.bedlines:
        for n in range(region.start, region.end):
            cons[region.chrom][n] = "N"

    vcf_reader = VCF(args.maskvcf)
    for record in vcf_reader:
        for n in range(0, len(record.REF)):
            cons[record.CHROM][record.POS - 1 + n] = "N"

    fh = open(args.output, "w")
    for k in seqs.keys():
        fh.write(">%s\n" % (k))
        fh.write(("".join(cons[k])) + "\n")
    fh.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("reference")
    parser.add_argument("maskfile")
    parser.add_argument("maskvcf")
    parser.add_argument("output")
    args = parser.parse_args()
    go(args)


if __name__ == "__main__":
    main()
