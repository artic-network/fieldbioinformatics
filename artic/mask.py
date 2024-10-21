#!/usr/bin/env python

from Bio import SeqIO
from cyvcf2 import VCF
import argparse
import pandas as pd


def read_3col_bed(fn):
    # read the primer scheme into a pandas dataframe and run type, length and null checks
    bedfile = pd.read_csv(
        fn,
        sep="\t",
        header=None,
        names=["chrom", "start", "end"],
        dtype={"chrom": str, "start": int, "end": int},
        usecols=(0, 1, 2),
        skiprows=0,
    )
    return bedfile


def go(args):
    seqs = dict([(rec.id, rec) for rec in SeqIO.parse(open(args.reference), "fasta")])
    cons = {}
    for k in seqs.keys():
        cons[k] = list(seqs[k].seq)

    bedfile = read_3col_bed(args.maskfile)
    for _, region in bedfile.iterrows():
        for n in range(region["start"], region["end"]):
            cons[region["chrom"]][n] = "N"

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
