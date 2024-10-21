#!/usr/bin/env python

import sys
import subprocess
import os
from collections import defaultdict
from cyvcf2 import VCF


def read_vcf(fn):
    vcfinfo = {}
    vcf_reader = VCF(fn)
    for record in vcf_reader:
        vcfinfo[record.POS] = record
    return vcfinfo


def collect_depths(bamfile):
    if not os.path.exists(bamfile):
        raise SystemExit("bamfile %s doesn't exist" % (bamfile,))

    p = subprocess.Popen(["samtools", "depth", bamfile], stdout=subprocess.PIPE)
    out, err = p.communicate()
    depths = defaultdict(dict)
    for ln in out.decode("utf-8").split("\n"):
        if ln:
            contig, pos, depth = ln.split("\t")
            depths[int(pos)] = int(depth)
    return depths


def main():
    positions = {}

    for sample_tag in sys.argv[1:]:
        for vcfset in ["", ".primertrimmed"]:
            vcffn = "%s%s.vcf" % (sample_tag, vcfset)
            if not os.path.exists(vcffn):
                continue

            print(vcffn, file=sys.stderr)

            vcf_reader = VCF(vcffn)
            for record in vcf_reader:
                if len(record.ALT[0]) == 1 and len(record.REF) == 1:
                    positions[record.POS] = "snp"
                else:
                    positions[record.POS] = "indel"

    print("pos\tset\tsample\tvartype\tdepth\tsupportfraction\tbasecalledfrequency")

    # for run, samples in runs.iteritems():
    #    for sample_tag in samples.keys():
    for sample_tag in sys.argv[1:]:
        for vcfset in ["", ".primertrimmed"]:
            vcffn = "%s%s.vcf" % (sample_tag, vcfset)
            if not os.path.exists(vcffn):
                print("%s does not exist" % (vcffn))
                continue

            vcffile = read_vcf(vcffn)
            bamfn = "%s.primertrimmed.sorted.bam" % (sample_tag)
            depths = collect_depths(bamfn)

            # 1-based pyvcf
            for pos, variant_type in positions.items():
                if pos - 1 in depths:
                    depth = depths[pos - 1]
                else:
                    depth = 0

                if pos in vcffile:
                    info = vcffile[pos].INFO
                    print(
                        "%s\t%s\t%s\t%s\t%s\t%s\t%s"
                        % (
                            pos,
                            vcfset,
                            sample_tag,
                            variant_type,
                            depth,
                            info["SupportFraction"],
                            info["BaseCalledFraction"],
                        )
                    )
                else:
                    print(
                        "%s\t%s\t%s\tinvariant\t%s\t0\t0"
                        % (pos, vcfset, sample_tag, depth)
                    )


if __name__ == "__main__":
    main()
