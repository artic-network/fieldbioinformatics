#Written by Nick Loman (@pathogenomenick)

import os
import sys
from Bio import SeqIO
from clint.textui import colored, puts, indent

def run(parser, args):
    # scheme
    # sample
    # read_file

    if args.scheme.find('/') != -1:
        scheme_name, scheme_version = args.scheme.split('/')
    else:
        scheme_name = args.scheme
        scheme_version = "V1"

    ref = "%s/%s/%s/%s.reference.fasta" % (args.scheme_directory, scheme_name, scheme_version, scheme_name)
    bed = "%s/%s/%s/%s.scheme.bed" % (args.scheme_directory, scheme_name, scheme_version, scheme_name)

    if args.read_file:
        read_file = args.read_file
    else:
        read_file = "%s.fasta" % (args.sample)

    if not os.path.exists(ref):
        print((colored.red('Scheme reference file not found: ') + ref))
        raise SystemExit
    if not os.path.exists(bed):
        print((colored.red('Scheme BED file not found: ') + bed))
        raise SystemExit

    cmds.append("bwa index %s" % (ref,))
    cmds.append("bwa mem -t %s -x ont2d %s %s | samtools view -bS - | samtools sort -o %s.sorted.bam -" % (args.threads, ref, read_file, args.sample))
    cmds.append("samtools index %s.sorted.bam" % (args.sample,))



