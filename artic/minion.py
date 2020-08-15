#Written by Nick Loman (@pathogenomenick)

from clint.textui import colored, puts, indent
from Bio import SeqIO
import os
import sys
import time
from .vcftagprimersites import read_bed_file

def get_nanopolish_header(ref):
    recs = list(SeqIO.parse(open(ref), "fasta"))
    if len (recs) != 1:
        print("FASTA has more than one sequence", file=sys.stderr)
        raise SystemExit(1)

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
        print(colored.red('Scheme reference file not found: ') + ref)
        raise SystemExit(1)
    if not os.path.exists(bed):
        print(colored.red('Scheme BED file not found: ') + bed)
        raise SystemExit(1)

    pools = set([row['PoolName'] for row in read_bed_file(bed)])

    cmds = []

    nanopolish_header = get_nanopolish_header(ref)

    if not args.medaka and not args.skip_nanopolish:
        if not args.fast5_directory or not args.sequencing_summary:
              print(colored.red('Must specify FAST5 directory and sequencing summary for nanopolish mode.'))
              raise SystemExit(1)
        
        cmds.append("nanopolish index -s %s -d %s %s" % (args.sequencing_summary, args.fast5_directory, args.read_file,))

    # 3) index the ref & align with minimap or bwa
    if not args.bwa:
        cmds.append("minimap2 -a -x map-ont -t %s %s %s | samtools view -bS -F 4 - | samtools sort -o %s.sorted.bam -" % (args.threads, ref, read_file, args.sample))
    else:
        cmds.append("bwa index %s" % (ref,))
        cmds.append("bwa mem -t %s -x ont2d %s %s | samtools view -bS -F 4 - | samtools sort -o %s.sorted.bam -" % (args.threads, ref, read_file, args.sample))
    cmds.append("samtools index %s.sorted.bam" % (args.sample,))

    # 4) trim the alignments to the primer start sites and normalise the coverage to save time
    if args.normalise:
        normalise_string = '--normalise %d' % (args.normalise)
    else:
        normalise_string = ''
    cmds.append("align_trim --start %s %s --report %s.alignreport.txt < %s.sorted.bam 2> %s.alignreport.er | samtools sort -T %s - -o %s.trimmed.rg.sorted.bam" % (normalise_string, bed, args.sample, args.sample, args.sample, args.sample, args.sample))
    cmds.append("align_trim %s %s --remove-incorrect-pairs --report %s.alignreport.txt < %s.sorted.bam 2> %s.alignreport.er | samtools sort -T %s - -o %s.primertrimmed.rg.sorted.bam" % (normalise_string, bed, args.sample, args.sample, args.sample, args.sample, args.sample))
    cmds.append("samtools index %s.trimmed.rg.sorted.bam" % (args.sample))
    cmds.append("samtools index %s.primertrimmed.rg.sorted.bam" % (args.sample))

    # 6) do variant calling on each read group using the raw signal alignment
    if args.medaka:
        for p in pools:
            if os.path.exists("%s.%s.hdf" % (args.sample, p)):
                os.remove("%s.%s.hdf" % (args.sample, p))
            cmds.append("medaka consensus --model %s --threads %s --chunk_len 800 --chunk_ovlp 400 --RG %s %s.trimmed.rg.sorted.bam %s.%s.hdf" % (args.medaka_model, args.threads, p, args.sample, args.sample, p))
            if args.no_indels:
                cmds.append("medaka snp %s %s.%s.hdf %s.%s.medaka.vcf" % (ref, args.sample, p, args.sample, p))
            else:
                cmds.append("medaka variant %s %s.%s.hdf %s.%s.medaka.vcf" % (ref, args.sample, p, args.sample, p))
            
            # annotate VCF with read depth info
            cmds.append("medaka tools annotate --pad 25 --RG %s %s.%s.medaka.vcf %s %s.trimmed.rg.sorted.bam %s.%s.vcf" % (p, args.sample, p, ref, args.sample, args.sample, p))

    else:
        if not args.skip_nanopolish:
            indexed_nanopolish_file = read_file
            if args.no_indels:
                nanopolish_extra_args = " --snps"
            else:
                nanopolish_extra_args = ""
            for p in pools:
                cmds.append("nanopolish variants --min-flanking-sequence 10 -x %s --progress -t %s --reads %s -o %s.%s.vcf -b %s.trimmed.rg.sorted.bam -g %s -w \"%s\" --ploidy 1 -m 0.15 --read-group %s %s" % (args.max_haplotypes, args.threads, indexed_nanopolish_file, args.sample, p, args.sample, ref, nanopolish_header, p, nanopolish_extra_args))

    # merge the called variants for each read group
    merge_vcf_cmd = "artic_vcf_merge %s %s" % (args.sample, bed)
    for p in pools:
        merge_vcf_cmd += " %s:%s.%s.vcf" % (p, args.sample, p)
    cmds.append(merge_vcf_cmd)

    # add in longshot here for medaka annotate results, just to apply a filter to get tests passing
    if args.medaka and not args.no_longshot:
        cmds.append("bgzip -f %s.merged.vcf" % (args.sample))
        cmds.append("tabix -p vcf %s.merged.vcf.gz" % (args.sample))
        cmds.append("longshot -P 0 -F -A --no_haps --bam %s.primertrimmed.rg.sorted.bam --ref %s --out %s.merged.vcf --potential_variants %s.merged.vcf.gz" % (args.sample, ref, args.sample, args.sample))

    # set up some name holder vars for ease
    if args.medaka:
        method = 'medaka'
    else:
        method = 'nanopolish'
    vcf_file = "%s.pass.vcf" % (args.sample)

    # filter the variants to produce PASS and FAIL lists, then index them
    cmds.append("artic_vcf_filter --%s %s %s.merged.vcf %s.pass.vcf %s.fail.vcf" % (method, bed, args.sample, args.sample, args.sample))
    cmds.append("bgzip -f %s" % (vcf_file))
    cmds.append("tabix -p vcf %s.gz" % (vcf_file))

    # get the depth of coverage for each readgroup and create a coverage mask and plots
    cmds.append("artic_make_depth_mask --store-rg-depths %s %s.primertrimmed.rg.sorted.bam %s.coverage_mask.txt" % (ref, args.sample, args.sample))
    cmds.append("artic_plot_amplicon_depth --primerScheme %s --sampleID %s --outFilePrefix %s %s*.depths" % (bed, args.sample, args.sample, args.sample))

    # add the failed variants to the coverage mask
    # artic_mask must be run before bcftools consensus
    cmds.append("artic_mask %s %s.coverage_mask.txt %s.fail.vcf %s.preconsensus.fasta" % (ref, args.sample, args.sample, args.sample))
    cmds.append("bcftools consensus -f %s.preconsensus.fasta %s.gz -m %s.coverage_mask.txt -o %s.consensus.fasta" % (args.sample, vcf_file, args.sample, args.sample))

    # apply header to the consensus sequence and run alignment
    fasta_header = "%s/ARTIC/%s" % (args.sample, method)
    cmds.append("artic_fasta_header %s.consensus.fasta \"%s\"" % (args.sample, fasta_header))
    cmds.append("cat %s.consensus.fasta %s > %s.muscle.in.fasta" % (args.sample, ref, args.sample))
    cmds.append("muscle -in %s.muscle.in.fasta -out %s.muscle.out.fasta" % (args.sample, args.sample))

    for cmd in cmds:
        print(colored.green("Running: ") + cmd, file=sys.stderr)
        if not args.dry_run:
            timerStart = time.perf_counter()
            retval = os.system(cmd)
            if retval != 0:
                print(colored.red('Command failed:' ) + cmd, file=sys.stderr)
                raise SystemExit(20)
            timerStop = time.perf_counter()

            # print the executed command and the runtime to the log file
            print("{}\t{}" .format(cmd, timerStop-timerStart), file=logfh)
        
        # if it's a dry run, print just the command
        else:
            print(cmd, file=logfh)
    logfh.close()
