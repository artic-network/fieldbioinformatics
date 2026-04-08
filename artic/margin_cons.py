#!/usr/bin/env python
from Bio import SeqIO
import sys
import pysam
import subprocess
from collections import defaultdict
import os.path
import argparse

# MASKED_POSITIONS = [2282]
MASKED_POSITIONS = []


def collect_depths(bamfile):
    if not os.path.exists(bamfile):
        raise SystemExit("bamfile %s doesn't exist" % (bamfile,))

    print(bamfile, file=sys.stderr)

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
        pos = r.pos + 1
        idfile = os.path.basename(self.vcffile).split(".")[0]
        print("%s\t%s\tstatus\t%s" % (idfile, pos, status), file=sys.stderr)
        print(
            "%s\t%s\tdepth\t%s" % (idfile, pos, r.info.get("TotalReads", ["n/a"])),
            file=sys.stderr,
        )
        print(
            "%s\t%s\tbasecalledfrac\t%s"
            % (idfile, pos, r.info.get("BaseCalledFraction", ["n/a"])),
            file=sys.stderr,
        )
        print(
            "%s\t%s\tsupportfrac\t%s"
            % (idfile, pos, r.info.get("SupportFraction", ["n/a"])),
            file=sys.stderr,
        )
        print("%s\t%s\tallele\t%s" % (idfile, pos, allele), file=sys.stderr)
        print("%s\t%s\tref\t%s" % (idfile, pos, r.ref), file=sys.stderr)


def go(args):

    depths = collect_depths(args.bamfile)
    reporter = Reporter(args.vcffile, depths)

    cons = ""

    seq = list(SeqIO.parse(open(sys.argv[1]), "fasta"))[0]
    cons = list(seq.seq)

    for n, c in enumerate(cons):
        try:
            depth = depths[seq.id][n + 1]
        except:
            depth = 0

        if depth < args.depth:
            cons[n] = "N"

    for mask in MASKED_POSITIONS:
        cons[mask - 1] = "N"

    sett = set()
    with pysam.VariantFile(args.vcffile) as vcf_reader:
        for record in vcf_reader:
            pos = record.pos + 1  # convert to 1-based
            alt = (record.alts or (".",))[0]
            if alt != ".":
                # variant call

                if pos in MASKED_POSITIONS:
                    reporter.report(record, "masked_manual", "n")
                    continue

                if "PRIMER" in record.info:
                    reporter.report(record, "primer_binding_site", "n")
                    cons[pos - 1] = "N"
                    continue

                support = float(record.info["SupportFraction"])
                total_reads = int(record.info["TotalReads"])
                qual = record.qual

                REF = record.ref
                ALT = str(alt)

                if len(ALT) > len(REF):
                    print(
                        "Skipping insertion at position: %s" % (pos), file=sys.stderr
                    )
                    continue

                if qual >= 200 and total_reads >= 20:
                    if len(REF) > len(ALT):
                        print(
                            "N-masking confident deletion at %s" % (pos),
                            file=sys.stderr,
                        )
                        for n in range(len(REF)):
                            cons[pos - 1 + n] = "N"
                        continue

                    reporter.report(record, "variant", ALT)
                    sett.add(pos)
                    if len(REF) > len(ALT):
                        print("deletion", file=sys.stderr)
                        continue

                    if len(ALT) > len(REF):
                        print("insertion", file=sys.stderr)
                        continue

                    cons[pos - 1] = str(ALT)
                elif len(REF) > len(ALT):
                    continue
                else:
                    reporter.report(record, "low_qual_variant", "n")
                    cons[pos - 1] = "N"
                    continue

    print(">%s" % (sys.argv[3]))
    print("".join(cons))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--depth", type=int, default=20, help="minimum depth to call a variant"
    )
    parser.add_argument("reference")
    parser.add_argument("vcffile")
    parser.add_argument("bamfile")
    args = parser.parse_args()
    go(args)


if __name__ == "__main__":
    main()
