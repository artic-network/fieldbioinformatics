import vcf
import sys
from operator import attrgetter
from collections import defaultdict
from .vcftagprimersites import read_bed_file

def in_frame(v):
    if len(v.ALT) > 1:
       print ("This code does not support multiple genotypes!")
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
        total_reads = float(v.INFO['TotalReads'])
        qual = v.QUAL
        strandbias = float(v.INFO['StrandFisherTest'])

        if qual / total_reads < 3:
            return False

        if self.no_frameshifts and not in_frame(v):
            return False

        if v.is_indel:
            strand_fraction_by_strand = v.INFO['SupportFractionByStrand']
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

    def check_filter(self, v):
        depth = v.INFO['DP']
        if depth < 20:
            return False

        if self.no_frameshifts and not in_frame(v):
            return False

        if v.num_het:
            return False
        return True

def go(args):
    vcf_reader = vcf.Reader(filename=args.inputvcf)
    vcf_writer = vcf.Writer(open(args.output_pass_vcf, 'w'), vcf_reader)
    vcf_writer_filtered = vcf.Writer(open(args.output_fail_vcf, 'w'), vcf_reader)
    if args.nanopolish:
        filter = NanoporeFilter(args.no_frameshifts)
    elif args.medaka:
        filter = MedakaFilter(args.no_frameshifts)
    else:
        print("Please specify a VCF type, i.e. --nanopolish or --medaka\n")
        raise SystemExit

    variants = [v for v in vcf_reader]

    group_variants = defaultdict(list)
    for v in variants:
        indx = "%s-%s" % (v.CHROM, v.POS)
        group_variants[indx].append(v)
    
    for v in variants:

        # if using medaka, we need to do a quick pre-filter to remove rubbish that we don't want adding to the mask
        if args.medaka:
            if v.INFO['DP'] <= 1:
                continue
            if v.QUAL < 20:
                continue

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
                print ("Suppress variant %s\n" % (v.POS))

def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--nanopolish', action='store_true')
    parser.add_argument('--medaka', action='store_true')
    parser.add_argument('--no-frameshifts', action='store_true')
    parser.add_argument('inputvcf')
    parser.add_argument('output_pass_vcf')
    parser.add_argument('output_fail_vcf')

    args = parser.parse_args()

    go(args)

if __name__ == "__main__":
    main()


