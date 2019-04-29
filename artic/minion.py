#Written by Nick Loman (@pathogenomenick)

import os
import sys
from Bio import SeqIO
from clint.textui import colored, puts, indent

def get_nanopolish_header(ref):
    recs = list(SeqIO.parse(open(ref), "fasta"))
    if len (recs) != 1:
        print(("FASTA has more than one sequence", sys.stderr))
        raise SystemExit

    return  "%s:%d-%d" % (recs[0].id, 1, len(recs[0])+1)

def run(parser, args):
    log = "%s.minion.log.txt" % (args.sample)
    logfh = open(log, 'w')

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

    cmds = []

    nanopolish_header = get_nanopolish_header(ref)

    # 3) index the ref & align with bwa"
    cmds.append("bwa index %s" % (ref,))
    cmds.append("bwa mem -t %s -x ont2d %s %s | samtools view -bS - | samtools sort -o %s.sorted.bam -" % (args.threads, ref, read_file, args.sample))
    cmds.append("samtools index %s.sorted.bam" % (args.sample,))

    # 4) trim the alignments to the primer start sites and normalise the coverage to save time
    if args.normalise:
        normalise_string = '--normalise %d' % (args.normalise)
    else:
        normalise_string = ''
    cmds.append("align_trim --start %s %s --report %s.alignreport.txt < %s.sorted.bam 2> %s.alignreport.er | samtools view -bS - | samtools sort -T %s - -o %s.trimmed.sorted.bam" % (normalise_string, bed, args.sample, args.sample, args.sample, args.sample, args.sample))
    cmds.append("align_trim %s %s --report %s.alignreport.txt < %s.sorted.bam 2> %s.alignreport.er | samtools view -bS - | samtools sort -T %s - -o %s.primertrimmed.sorted.bam" % (normalise_string, bed, args.sample, args.sample, args.sample, args.sample, args.sample))
    cmds.append("samtools index %s.trimmed.sorted.bam" % (args.sample))
    cmds.append("samtools index %s.primertrimmed.sorted.bam" % (args.sample))

    #covplot.R $sample.alignreport.txt

    # 6) do variant calling using the raw signal alignment
    if not args.skip_nanopolish:
        if args.nanopolish_read_file:
                indexed_nanopolish_file = args.nanopolish_read_file
        else:
                indexed_nanopolish_file = read_file

        cmds.append("nanopolish variants --min-flanking-sequence 10 --fix-homopolymers -x %s --progress -t %s --reads %s -o %s.vcf -b %s.trimmed.sorted.bam -g %s -w \"%s\"  --snps --ploidy 1" % (args.max_haplotypes, args.threads, indexed_nanopolish_file, args.sample, args.sample, ref, nanopolish_header))
        cmds.append("nanopolish variants --fix-homopolymers -x %s --progress -t %s --reads %s -o %s.primertrimmed.vcf -b %s.primertrimmed.sorted.bam -g %s -w \"%s\" --snps --ploidy 1" % (args.max_haplotypes, args.threads, indexed_nanopolish_file, args.sample, args.sample, ref, nanopolish_header))

        #python nanopore-scripts/expand-cigar.py --bam "$sample".primertrimmed.sorted.bam --fasta $ref | python nanopore-scripts/count-errors.py /dev/stdin > "$sample".errors.txt

        # 7) do phasing
        #nanopolish phase-reads --reads $sample.fasta --bam $sample.trimmed.sorted.bam --genome $ref $sample.vcf

        # 8) variant frequency plot
        cmds.append("vcfextract %s > %s.variants.tab" % (args.sample, args.sample))

        # 8) filter the variants and produce a consensus
        # here we use the vcf file without primer binding site trimming (to keep nanopolish happy with flanks)
        # but we use the primertrimmed sorted bam file in order that primer binding sites do not count
        # for the depth calculation to determine any low coverage sites that need masking
        cmds.append("margin_cons %s %s.vcf %s.primertrimmed.sorted.bam a > %s.consensus.fasta" % (ref, args.sample, args.sample, args.sample))

    for cmd in cmds:
        print((colored.green("Running: ") + cmd, sys.stderr))
        print((cmd, logfh))
        retval = os.system(cmd)
        if retval != 0:
            print((colored.red('Command failed:' ) + cmd, sys.stderr))

    logfh.close()

