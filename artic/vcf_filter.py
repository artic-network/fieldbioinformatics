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
    def __init__(
        self,
        no_frameshifts,
        min_depth,
        min_variant_quality=10,
        min_frameshift_quality=50,
        min_allele_frequency=0.6,
        min_mask_allele_frequency=0.1,
        min_minor_allele_count=5,
    ):
        self.no_frameshifts = no_frameshifts
        self.min_depth = min_depth
        self.min_variant_quality = min_variant_quality
        self.min_frameshift_quality = min_frameshift_quality
        self.min_allele_frequency = min_allele_frequency
        self.min_mask_allele_frequency = min_mask_allele_frequency
        self.min_minor_allele_count = min_minor_allele_count

    def check_filter(self, v):
        """Return 'pass', 'mask', or 'discard'.

        'mask'    — position is ambiguous; replace with N in consensus
                    (low qual, low AF, frameshift, or low depth)
        'discard' — variant is a likely artefact but position is fine;
                    keep reference base (only: low minor allele count)
        'pass'    — variant passes all filters
        """
        qual = v.qual

        if qual < self.min_variant_quality:
            return "mask"

        try:
            allele_freq = list(v.samples.values())[0]["AF"]
            if isinstance(allele_freq, tuple):
                allele_freq = allele_freq[0]
        except Exception:
            print(
                f"ERROR: Could not find AF for variant at {v.chrom}:{v.pos}, cannot filter on allele frequency"
            )
            raise SystemExit(1)

        if allele_freq < self.min_mask_allele_frequency:
            return "discard"

        if allele_freq < self.min_allele_frequency:
            return "mask"

        if not in_frame(v):
            if self.no_frameshifts:
                return "mask"
            # require a higher quality for frameshifting indels, they're far more likely to be errors
            if qual < self.min_frameshift_quality:
                return "mask"

        try:
            depth = list(v.samples.values())[0]["DP"]
            if depth < self.min_depth:
                return "mask"
        except Exception:
            pass

        # Minor allele count: high-quality call with too few supporting reads is a
        # likely artefact — discard rather than mask so the reference base is kept.
        try:
            ad = list(v.samples.values())[0]["AD"]
            if ad[1] < self.min_minor_allele_count:
                return "discard"
        except Exception:
            pass

        return "pass"


def go(args):
    vcf_reader = pysam.VariantFile(args.inputvcf)
    vcf_writer = pysam.VariantFile(args.output_pass_vcf, "w", header=vcf_reader.header)
    vcf_writer_filtered = pysam.VariantFile(
        args.output_fail_vcf, "w", header=vcf_reader.header
    )
    vcf_writer_ignore = pysam.VariantFile(
        args.output_ignore_vcf, "w", header=vcf_reader.header
    )
    filter_obj = Clair3Filter(
        args.no_frameshifts,
        args.min_depth,
        args.min_variant_quality,
        args.min_frameshift_quality,
        args.min_allele_frequency,
        args.min_mask_allele_frequency,
        args.min_minor_allele_count,
    )

    variants = [v for v in vcf_reader]

    group_variants = defaultdict(list)
    for v in variants:
        index = "%s-%s" % (v.chrom, v.pos)
        group_variants[index].append(v)

    for v in variants:

        # Pre-filter: DP <= 1 means essentially no information — mask the position.
        try:
            if v.info["DP"] <= 1:
                print(f"Mask variant {v.pos} due to very low depth (DP<=1)")
                vcf_writer_filtered.write(v)
                continue

        except KeyError:
            pass

        result = filter_obj.check_filter(v)

        if result == "pass":
            vcf_writer.write(v)

        elif result == "mask":
            index = "%s-%s" % (v.chrom, v.pos)
            if len(group_variants[index]) > 1:
                if any(filter_obj.check_filter(ov) == "pass" for ov in group_variants[index]):
                    print("Suppress variant %s\n" % (v.pos))
                    continue
            vcf_writer_filtered.write(v)

        else:  # "discard"
            index = "%s-%s" % (v.chrom, v.pos)
            any_pass = (
                len(group_variants[index]) > 1
                and any(filter_obj.check_filter(ov) == "pass" for ov in group_variants[index])
            )
            if any_pass:
                print("Suppress variant %s\n" % (v.pos))
            else:
                vcf_writer_ignore.write(v)


def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--no-frameshifts", action="store_true")
    parser.add_argument("--min-depth", type=int)
    parser.add_argument("--min-variant-quality", type=int, default=10)
    parser.add_argument("--min-frameshift-quality", type=int, default=50)
    parser.add_argument("--min-allele-frequency", type=float, default=0.6)
    parser.add_argument("--min-mask-allele-frequency", type=float, default=0.1)
    parser.add_argument("--min-minor-allele-count", type=int, default=5)

    parser.add_argument("inputvcf")
    parser.add_argument("output_pass_vcf")
    parser.add_argument("output_fail_vcf")
    parser.add_argument("output_ignore_vcf")

    args = parser.parse_args()

    go(args)


if __name__ == "__main__":
    main()
