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


class NanoporeFilter:
    def __init__(self, no_frameshifts):
        self.no_frameshifts = no_frameshifts
        pass

    def check_filter(self, v):
        total_reads = float(v.INFO["TotalReads"])
        qual = v.QUAL
        # strandbias = float(v.INFO["StrandFisherTest"])

        if qual / total_reads < 3:
            return False

        if self.no_frameshifts and not in_frame(v):
            return False

        if v.is_indel:
            strand_fraction_by_strand = v.INFO["SupportFractionByStrand"]
            if float(strand_fraction_by_strand[0]) < 0.5:
                return False

            if float(strand_fraction_by_strand[1]) < 0.5:
                return False

        if total_reads < 20:
            return False

        return True


class MedakaFilter:
    def __init__(self, no_frameshifts):
        self.no_frameshifts = no_frameshifts

    def check_filter(self, v, min_depth):
        try:
            # We don't really care about the depth here, just skip it if it isn't there
            depth = v.INFO["DP"]
        except KeyError:
            depth = v.format("DP")[0][0]
        
        if depth < min_depth:
            return False

        if self.no_frameshifts and not in_frame(v):
            return False

        if v.num_het:
            return False
        return True


def go(args):
    vcf_reader = VCF(args.inputvcf)
    vcf_writer = Writer(args.output_pass_vcf, vcf_reader, "w")
    vcf_writer.write_header()
    vcf_writer_filtered = Writer(args.output_fail_vcf, vcf_reader, "w")
    vcf_writer_filtered.write_header()
    filter = MedakaFilter(args.no_frameshifts)

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

        if v.QUAL < args.min_variant_quality:
            print(f"Suppress variant {v.POS} due to low quality")
            continue

        # now apply the filter to send variants to PASS or FAIL file
        if filter.check_filter(v, args.min_depth):
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
    parser.add_argument("--min-variant-quality", type=int)
    parser.add_argument("--min-depth", type=int)
    parser.add_argument("inputvcf")
    parser.add_argument("output_pass_vcf")
    parser.add_argument("output_fail_vcf")

    args = parser.parse_args()

    go(args)


if __name__ == "__main__":
    main()
