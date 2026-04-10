import pysam
from collections import defaultdict


def in_frame(v):
    alts = v.alts or ()
    if len(alts) > 1:
        print("This code does not support multiple genotypes!")
        raise SystemExit
    ref = v.ref

    if not alts:  # No ALT alleles (e.g. a deletion)
        if len(ref) % 3 == 0:
            return True
        return False

    alt = alts[0]
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
        self.min_frameshift_quality = 50
        self.min_allele_frequency = 0.6

    def check_filter(self, v):
        qual = v.qual

        if qual < self.min_variant_quality:
            return False

        # Filter out low allele frequency variants
        try:
            allele_freq = list(v.samples.values())[0]["AF"]
            if isinstance(allele_freq, tuple):
                allele_freq = allele_freq[0]
        except Exception:
            print(
                f"ERROR: Could not find AF for variant at {v.chrom}:{v.pos}, cannot filter on allele frequency"
            )
            raise SystemExit(1)

        if allele_freq < self.min_allele_frequency:
            return False

        if not in_frame(v):
            if self.no_frameshifts:
                return False
            # require a higher quality for frameshifting indels, they're far more likely to be errors
            if qual < self.min_frameshift_quality:
                return False

        try:
            # We don't really care about the depth here, just skip it if it isn't there
            depth = list(v.samples.values())[0]["DP"]

            if depth < self.min_depth:
                return False

        except Exception:
            pass

        return True


def go(args):
    vcf_reader = pysam.VariantFile(args.inputvcf)
    vcf_writer = pysam.VariantFile(args.output_pass_vcf, "w", header=vcf_reader.header)
    vcf_writer_filtered = pysam.VariantFile(
        args.output_fail_vcf, "w", header=vcf_reader.header
    )
    filter = Clair3Filter(args.no_frameshifts, args.min_depth)

    variants = [v for v in vcf_reader]

    group_variants = defaultdict(list)
    for v in variants:
        indx = "%s-%s" % (v.chrom, v.pos)
        group_variants[indx].append(v)

    for v in variants:

        # quick pre-filter to remove rubbish that we don't want adding to the mask
        try:
            if v.info["DP"] <= 1:
                print(f"Suppress variant {v.pos} due to low depth")
                continue

        except KeyError:
            pass

        # now apply the filter to send variants to PASS or FAIL file
        if filter.check_filter(v):
            vcf_writer.write(v)
        else:
            variant_passes = False

            indx = "%s-%s" % (v.chrom, v.pos)
            if len(group_variants[indx]) > 1:
                for check_variant in group_variants[indx]:
                    if filter.check_filter(check_variant):
                        variant_passes = True

            if not variant_passes:
                vcf_writer_filtered.write(v)

            else:
                print("Suppress variant %s\n" % (v.pos))


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
