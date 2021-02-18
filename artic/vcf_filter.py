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
    def __init__(self, no_frameshifts, min_depth, min_qual, qual_cov_ratio):
        self.no_frameshifts = no_frameshifts
        self.min_depth = min_depth
        self.min_qual = min_qual
        self.qual_cov_ratio = qual_cov_ratio
        pass

    def check_filter(self, v):
        total_reads = float(v.INFO['TotalReads'])
        qual = v.QUAL

        if (qual / total_reads) < self.qual_cov_ratio:
            return False

        if self.no_frameshifts and not in_frame(v):
            return False

        if v.is_indel:
            strand_fraction_by_strand = v.INFO['SupportFractionByStrand']
            if float(strand_fraction_by_strand[0]) < 0.5: 
                return False

            if float(strand_fraction_by_strand[1]) < 0.5:
                return False

        if total_reads < self.min_depth:
            return False
        
        if qual < self.min_qual:
            return False

        return True

class MedakaFilter:
    def __init__(self, no_frameshifts, min_depth, min_qual):
        self.no_frameshifts = no_frameshifts
        self.min_depth = min_depth
        self.min_qual = min_qual

    def check_filter(self, v):
        total_reads = float(v.INFO['DP'])
        qual = v.QUAL
        # Medaka genotype quality score (GQ) 
        ## NOTE: this isn't very robust way to retrieve but will do for our use case.
        medaka_score = float(v.samples[0].data.GQ)

        if self.no_frameshifts and not in_frame(v):
            return False

        if v.num_het:
            return False

        if total_reads < self.min_depth:
            return False

        if (medaka_score < self.min_qual or qual < self.min_qual): 
            return False

        return True

def go(args):
    vcf_reader = vcf.Reader(filename=args.inputvcf)
    vcf_writer = vcf.Writer(open(args.output_pass_vcf, 'w'), vcf_reader)
    vcf_writer_filtered = vcf.Writer(open(args.output_fail_vcf, 'w'), vcf_reader)
    if args.nanopolish:
        filter = NanoporeFilter(args.no_frameshifts, args.min_depth, args.min_qual, args.nanopolish_qual_cov_ratio)
    elif args.medaka:
        filter = MedakaFilter(args.no_frameshifts, args.min_depth, args.min_qual)
    else:
        print("Please specify a VCF type, i.e. --nanopolish or --medaka\n")
        raise SystemExit

    variants = [v for v in vcf_reader]

    group_variants = defaultdict(list)
    for v in variants:
        indx = "%s-%s" % (v.CHROM, v.POS)
        group_variants[indx].append(v)
    
    for v in variants:

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
    parser.add_argument('--min-depth', required=False, type=int, default=20)
    parser.add_argument('--min-qual', required=False, type=int, default=20)
    parser.add_argument('--nanopolish-qual-cov-ratio', required=False, type=float, default=3.0)
    parser.add_argument('inputvcf')
    parser.add_argument('output_pass_vcf')
    parser.add_argument('output_fail_vcf')

    args = parser.parse_args()

    go(args)

if __name__ == "__main__":
    main()


