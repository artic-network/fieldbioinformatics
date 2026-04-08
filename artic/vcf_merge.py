import pysam
import sys
from collections import defaultdict

from primalbedtools.scheme import Scheme
from primalbedtools.bedfiles import merge_primers


def _has_variants(fn):
    with pysam.VariantFile(fn) as vcf:
        return next(iter(vcf), None) is not None


def vcf_merge(args):

    try:
        scheme = Scheme.from_file(args.bedfile)
        scheme.bedlines = merge_primers(scheme.bedlines)
    except (ValueError, TypeError) as e:
        print(
            f"Failed to parse primer scheme BED file '{args.bedfile}': {e}",
            file=sys.stderr,
        )
        raise SystemExit(3)

    primer_map = defaultdict(dict)

    for p in scheme.bedlines:
        for n in range(p.start, p.end + 1):
            primer_map[p.pool][n] = p.primername

    template_header = None

    pool_map = {}
    for param in args.vcflist:
        pool_name, file_name = param.split(":")
        pool_map[file_name] = pool_name
        if not template_header:
            if _has_variants(file_name):
                with pysam.VariantFile(file_name) as vcf_reader:
                    template_header = vcf_reader.header.copy()
            else:
                print(
                    f"Not using {file_name} as VCF template since it has no variants",
                    file=sys.stderr,
                )

    template_header.info.add("Pool", 1, "String", "The pool name")

    vcf_writer = pysam.VariantFile(f"{args.prefix}.merged.vcf", "w", header=template_header)
    vcf_writer_primers = pysam.VariantFile(f"{args.prefix}.primers.vcf", "w", header=template_header)

    variants = []
    for file_name, pool_name in pool_map.items():
        if not _has_variants(file_name):
            print(f"Skipping {file_name} as it has no variants", file=sys.stderr)
            continue

        with pysam.VariantFile(file_name) as vcf_reader:
            vcf_reader.header.info.add("Pool", 1, "String", "The pool name")
            for v in vcf_reader:
                v.info["Pool"] = pool_name
                variants.append(v.copy())

    variants.sort(key=lambda v: (v.chrom, v.pos))

    for v in variants:
        if v.pos + 1 in primer_map[v.info["Pool"]]:
            vcf_writer_primers.write(v)
            print(
                "found primer binding site mismatch: %s"
                % (primer_map[v.info["Pool"]][v.pos + 1]),
                file=sys.stderr,
            )
        else:
            vcf_writer.write(v)


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Trim alignments from an amplicon scheme."
    )
    parser.add_argument("prefix")
    parser.add_argument("bedfile")
    parser.add_argument("vcflist", nargs="+")

    args = parser.parse_args()
    vcf_merge(args)


if __name__ == "__main__":
    main()
