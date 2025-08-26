from cyvcf2 import VCF, Writer
from collections import defaultdict


def in_frame(v):
    if len(v.ALT) > 1:
        print("This code does not support multiple genotypes!")
        raise SystemExit
    ref = v.REF
    alt = v.ALT[0]
    bases = len(alt) - len(ref)
    if not bases:
        return True
    if bases % 3 == 0:
        return True
    return False


class Clair3Filter:
    def __init__(self, no_frameshifts, min_depth):
        self.no_frameshifts = no_frameshifts
        self.min_depth = min_depth
        self.min_variant_quality = 10

    def check_filter(self, v):
        qual = v.QUAL

        if qual < self.min_variant_quality:
            return False

        if self.no_frameshifts and not in_frame(v):
            return False

        try:
            # We don't really care about the depth here, just skip it if it isn't there
            depth = v.INFO["DP"]
        except KeyError:
            depth = v.format("DP")[0][0]

        if depth < self.min_depth:
            return False

        return True


def go(args):
    vcf_reader = VCF(args.inputvcf)
    vcf_writer = Writer(args.output_pass_vcf, vcf_reader, "w")
    vcf_writer.write_header()
    vcf_writer_filtered = Writer(args.output_fail_vcf, vcf_reader, "w")
    vcf_writer_filtered.write_header()
    filter = Clair3Filter(args.no_frameshifts, args.min_depth)

    variants = [v for v in vcf_reader]

    group_variants = defaultdict(list)
    for v in variants:
        indx = "%s-%s" % (v.CHROM, v.POS)
        group_variants[indx].append(v)

    for v in variants:

        # quick pre-filter to remove rubbish that we don't want adding to the mask
        try:
            if v.INFO["DP"] <= 1:
                print(f"Suppress variant {v.POS} due to low depth")
                continue

        except KeyError:
            pass

        # now apply the filter to send variants to PASS or FAIL file
        if filter.check_filter(v):
            vcf_writer.write_record(v)
        else:
            variant_passes = False

            indx = "%s-%s" % (v.CHROM, v.POS)
            if len(group_variants[indx]) > 1:
                for check_variant in group_variants[indx]:
                    if filter.check_filter(check_variant):
                        variant_passes = True

            if not variant_passes:
                vcf_writer_filtered.write_record(v)

            else:
                print("Suppress variant %s\n" % (v.POS))


def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--no-frameshifts", action="store_true")
    parser.add_argument("--min-depth", type=int)
    parser.add_argument("inputvcf")
    parser.add_argument("output_pass_vcf")
    parser.add_argument("output_fail_vcf")

    args = parser.parse_args()

    go(args)


if __name__ == "__main__":
    main()
