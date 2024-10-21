#!/usr/bin/env python
from Bio import SeqIO
import sys
import vcf
import subprocess
from collections import defaultdict
import os.path
import operator
from .vcftagprimersites import read_bed_file
import argparse


def collect_depths(bamfile):
    if not os.path.exists(bamfile):
        raise SystemExit("bamfile %s doesn't exist" % (bamfile,))

    p = subprocess.Popen(["samtools", "depth", bamfile], stdout=subprocess.PIPE)
    out, err = p.communicate()
    depths = defaultdict(dict)
    for ln in out.decode("utf-8").split("\n"):
        if ln:
            contig, pos, depth = ln.split("\t")
            depths[contig][int(pos)] = int(depth)
    return depths


class Reporter:
    def __init__(self, vcffile, depths):
        self.vcffile = vcffile
        self.depths = depths

    def report(self, r, status, allele):
        idfile = os.path.basename(self.vcffile).split(".")[0]
        print("%s\t%s\tstatus\t%s" % (idfile, r.POS, status), file=sys.stderr)
        print("%s\t%s\tallele\t%s" % (idfile, r.POS, allele), file=sys.stderr)
        print("%s\t%s\tref\t%s" % (idfile, r.POS, r.REF), file=sys.stderr)
        print(
            "%s\t%s\tdepth\t%s" % (idfile, r.POS, self.depths[r.CHROM][r.POS]),
            file=sys.stderr,
        )


def go(args):
    MASKED_POSITIONS = defaultdict(set)

    depths = collect_depths(args.bamfile)
    reporter = Reporter(args.vcffile, depths)

    seqs = dict([(rec.id, rec) for rec in SeqIO.parse(open(args.reference), "fasta")])
    cons = {}
    for k in seqs.keys():
        cons[k] = list(seqs[k].seq)

        for n, c in enumerate(cons[k]):
            try:
                depth = depths[k][n + 1]
            except:
                depth = 0

            if depth < args.depth:
                cons[k][n] = "N"

    if args.masked:
        for region in args.masked.split(","):
            contig, positions = region.split(":")
            start, end = positions.split("-")
            start = int(start)
            end = int(end)
            for n in range(start, end):
                cons[contig][n - 1] = "N"
                MASKED_POSITIONS[contig].add(n - 1)

    sett = set()
    vcf_reader = vcf.Reader(open(args.vcffile, "r"))
    for record in vcf_reader:
        if record.ALT[0] != ".":
            # variant call

            if record.POS in MASKED_POSITIONS[record.CHROM]:
                reporter.report(record, "masked_manual", "n")
                cons[record.CHROM][record.POS - 1] = "N"
                continue

            if record.num_het:
                if depths[record.CHROM][record.POS] < args.depth:
                    reporter.report(record, "het_site_low_depth", "y")
                    continue
                else:
                    reporter.report(record, "het_site", "y")
                    cons[record.CHROM][record.POS - 1] = "N"
                    continue

            if "PRIMER" in record.INFO:
                reporter.report(record, "primer_binding_site", "n")
                cons[record.CHROM][record.POS - 1] = "N"
                continue

            # support = float(record.INFO['SupportFraction'])
            # total_reads = int(record.INFO['TotalReads'])
            qual = record.QUAL

            REF = record.REF
            ALT = str(record.ALT[0])

            if len(ALT) > len(REF):
                print(
                    "Skipping insertion at position: %s" % (record.POS), file=sys.stderr
                )
                continue

            if (
                depths[record.CHROM][record.POS] >= args.depth
                and record.QUAL >= args.quality
            ):
                if len(REF) > len(ALT):
                    print(
                        "N-masking confident deletion at %s" % (record.POS),
                        file=sys.stderr,
                    )
                    for n in range(len(REF)):
                        cons[record.CHROM][record.POS - 1 + n] = "N"
                    continue

                reporter.report(record, "variant", ALT)
                sett.add(record.POS)
                if len(REF) > len(ALT):
                    print("deletion", file=sys.stderr)
                    continue

                if len(ALT) > len(REF):
                    print("insertion", file=sys.stderr)
                    continue
                cons[record.CHROM][record.POS - 1] = str(ALT)
            elif len(REF) > len(ALT):
                continue
            else:
                if depths[record.CHROM][record.POS] < args.depth:
                    reporter.report(record, "low_depth_variant", "n")
                else:
                    reporter.report(record, "low_qual_variant", "n")
                # cons[record.CHROM][record.POS-1] = 'N'
                continue

    # print >>sys.stderr, str(sett)

    for k in seqs.keys():
        print(">%s-%s" % (args.bamfile, k))
        print("".join(cons[k]))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--depth", type=int, default=5, help="minimum depth to call a variant"
    )
    parser.add_argument(
        "--quality", type=int, default=0, help="minimum quality to call a variant"
    )
    parser.add_argument(
        "--masked", help="Regions to mask (contig:start-end,contig:start-end)"
    )
    parser.add_argument("reference")
    parser.add_argument("vcffile")
    parser.add_argument("bamfile")
    args = parser.parse_args()
    go(args)


if __name__ == "__main__":
    main()
