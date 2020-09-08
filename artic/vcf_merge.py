import vcf
import sys
from operator import attrgetter
from collections import defaultdict
from .vcftagprimersites import read_bed_file

def vcf_merge(args):
   bed = read_bed_file(args.bedfile)

   primer_map = defaultdict(dict)

   for p in bed:
      for n in range(p['start'], p['end']+1):
         primer_map[p['PoolName']][n] = p['Primer_ID']

   first_vcf = None

   pool_map = {}
   for param in args.vcflist:
      pool_name, file_name = param.split(":")
      pool_map[file_name] = pool_name
      if not first_vcf: 
         first_vcf = file_name

   vcf_reader = vcf.Reader(filename=first_vcf)
   vcf_reader.infos["Pool"] = vcf.parser._Format("Pool", 1, "String", "The pool name")
   vcf_writer = vcf.Writer(open(args.prefix+'.merged.vcf', 'w'), vcf_reader)
   vcf_writer_primers = vcf.Writer(open(args.prefix+'.primers.vcf', 'w'), vcf_reader)

   variants = []
   for file_name, pool_name in pool_map.items():
      vcf_reader = vcf.Reader(filename=file_name)
      for v in vcf_reader:
         v.INFO['Pool'] = pool_name
         variants.append(v)

   variants.sort(key=attrgetter('CHROM', 'POS'))

   for v in variants:
      if v.POS in primer_map[v.INFO['Pool']]:
         vcf_writer_primers.write_record(v)
         print("found primer binding site mismatch: %s" % (primer_map[v.INFO['Pool']][v.POS]), file=sys.stderr)
      else:
         vcf_writer.write_record(v)

def main():
    import argparse

    parser = argparse.ArgumentParser(description='Trim alignments from an amplicon scheme.')
    parser.add_argument('prefix')
    parser.add_argument('bedfile')
    parser.add_argument('vcflist', nargs='+')

    args = parser.parse_args()
    vcf_merge(args)

if __name__ == "__main__":
    main()

