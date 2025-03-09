# Written by Nick Loman (@pathogenomenick)

from clint.textui import colored
import os
import sys
import time
from artic.utils import read_bed_file, get_scheme, choose_model


def run(parser, args):

    if any([args.bed, args.ref]) and any(
        [
            args.scheme_name,
            args.scheme_version,
            args.scheme_length,
        ]
    ):
        print(
            colored.red(
                "You cannot mix 'Remote Scheme Options' with 'Local Scheme Options', for more information run 'artic minion --help'"
            ),
            file=sys.stderr,
        )
        raise SystemExit(1)

    if any([args.bed, args.ref]) and not all([args.bed, args.ref]):
        print(
            colored.red(
                "You must provide both a BED file and a reference sequence for local scheme options"
            ),
            file=sys.stderr,
        )
        raise SystemExit(1)

    if any([args.scheme_name, args.scheme_version, args.scheme_length]) and not all(
        [args.scheme_name, args.scheme_version]
    ):
        print(
            colored.red(
                "You must provide a scheme name and a scheme version at minimum for remote scheme options"
            ),
            file=sys.stderr,
        )
        raise SystemExit(1)

    if args.normalise > 10000:
        print(
            colored.red(
                "--normalise appears to be an aribtrarily high value, if you wish to not normalise your BAM files set '--normalise 0' to disable normalisation entirely"
            ),
            file=sys.stderr,
        )
        raise SystemExit(1)

    if args.model_dir:
        model_path = args.model_dir
    else:
        if not os.getenv("CONDA_PREFIX"):
            print(
                "CONDA_PREFIX is not set, this probably means you are not running this inside a conda environment, please provide a model path argument '--model-path'",
                file=sys.stderr,
            )
            raise SystemExit(1)

        model_path = f"{os.getenv('CONDA_PREFIX')}/bin/models/"

    # check for model
    if not args.model:
        model = choose_model(args.read_file)
        full_model_path = f"{model_path}/{model['name']}"

    else:
        full_model_path = f"{model_path}/{args.model}"

    if not os.path.exists(full_model_path):
        print(
            colored.red(
                f"Model '{model['name']}' not found in '{model_path}', please run 'artic_get_models' to download the clair3 models from ONT"
            ),
            file=sys.stderr,
        )
        raise SystemExit(1)

    # check the parameters and set up the filenames
    ## find the primer scheme, reference sequence and confirm scheme version

    if args.bed and args.ref:
        bed = args.bed
        ref = args.ref
    else:
        bed, ref, _ = get_scheme(
            scheme_name=args.scheme_name,
            scheme_version=args.scheme_version,
            scheme_directory=args.scheme_directory,
            scheme_length=args.scheme_length,
            read_file=args.read_file,
        )

    if not os.path.exists(bed) or not os.path.exists(ref):
        print(
            colored.red(
                "Failed to find primer scheme or reference sequence: {} {}".format(
                    bed, ref
                )
            ),
            file=sys.stderr,
        )
        raise SystemExit(1)

    ## set up the read file
    if args.read_file:
        read_file = args.read_file
    else:
        read_file = "%s.fasta" % (args.sample)
    if not os.path.exists(read_file):
        print(
            colored.red("failed to find read-file: {}".format(read_file)),
            file=sys.stderr,
        )
        raise SystemExit(1)

    ## collect the primer pools
    pools = set([row["PoolName"] for row in read_bed_file(bed)])

    ## create a holder to keep the pipeline commands in
    cmds = []

    # 2) linearise the reference if desired, check the reference fasta has an index and create one if not
    if args.linearise_fasta:
        linearised_ref = f"{args.sample}.linearised.fasta"
        # No f-strings here to reduce amount of escaping
        cmds.append(
            """awk '!/^>/ { printf "%s", $0; n = "\\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' """
            + ref
            + " > "
            + linearised_ref
        )
        ref = linearised_ref

    if not os.path.exists("%s.fai" % (ref)):
        cmds.append("samtools faidx %s" % (ref))

    # 3) index the ref & align with minimap
    cmds.append(
        f"minimap2 -a -x map-ont -t {args.threads} {ref} {read_file} | samtools view -bS -F 4 - | samtools sort -o {args.sample}.sorted.bam -"
    )

    cmds.append(f"samtools index {args.sample}.sorted.bam")

    # 4) trim the alignments to the primer start sites and normalise the coverage to save time
    if args.normalise:
        normalise_string = f"--normalise {args.normalise}"
    else:
        normalise_string = ""

    cmds.append(
        f"align_trim {normalise_string} {bed} --primer-match-threshold {args.primer_match_threshold} --remove-incorrect-pairs --min-mapq {args.min_mapq} --report {args.sample}.alignreport.csv < {args.sample}.sorted.bam > {args.sample}.trimmed.rg.sam"
    )

    cmds.append(
        f"samtools sort -T {args.sample} {args.sample}.trimmed.rg.sam -o {args.sample}.trimmed.rg.sorted.bam"
    )
    cmds.append(f"rm {args.sample}.trimmed.rg.sam")

    cmds.append(
        f"align_trim {normalise_string} {bed} --primer-match-threshold {args.primer_match_threshold} --min-mapq {args.min_mapq} --remove-incorrect-pairs --trim-primers --report {args.sample}.alignreport.csv --amp-depth-report {args.sample}.amplicon_depths.tsv < {args.sample}.sorted.bam > {args.sample}.primertrimmed.rg.sam"
    )

    cmds.append(
        f"samtools sort -T {args.sample} {args.sample}.primertrimmed.rg.sam -o {args.sample}.primertrimmed.rg.sorted.bam"
    )
    cmds.append(f"rm {args.sample}.primertrimmed.rg.sam")

    cmds.append(f"samtools index {args.sample}.trimmed.rg.sorted.bam")
    cmds.append(f"samtools index {args.sample}.primertrimmed.rg.sorted.bam")

    # 6) do variant calling on each read group
    for p in pools:
        if os.path.exists("%s.%s.hdf" % (args.sample, p)):
            os.remove("%s.%s.hdf" % (args.sample, p))

        # Split the BAM by read group
        for p in pools:
            cmds.append(
                f"samtools view -b -r {p} {args.sample}.trimmed.rg.sorted.bam -o {args.sample}.{p}.trimmed.rg.sorted.bam"
            )

            cmds.append(f"samtools index {args.sample}.{p}.trimmed.rg.sorted.bam")

            cmds.append(
                f"run_clair3.sh --enable_long_indel --chunk_size=10000 --haploid_precise --no_phasing_for_fa --bam_fn='{args.sample}.{p}.trimmed.rg.sorted.bam' --ref_fn='{ref}' --output='{args.sample}_rg_{p}' --threads='{args.threads}' --platform='ont' --model_path='{full_model_path}' --include_all_ctgs"
            )

            cmds.append(
                f"bgzip -dc {args.sample}_rg_{p}/merge_output.vcf.gz > {args.sample}.{p}.vcf"
            )

    # 7) merge the called variants for each read group
    merge_vcf_cmd = "artic_vcf_merge %s %s 2> %s.primersitereport.txt" % (
        args.sample,
        bed,
        args.sample,
    )
    for p in pools:
        merge_vcf_cmd += " %s:%s.%s.vcf" % (p, args.sample, p)

    cmds.append(merge_vcf_cmd)

    pre_filter_vcf = f"{args.sample}.merged.vcf"
    cmds.append(f"bgzip -kf {pre_filter_vcf}")
    cmds.append(f"tabix -f -p vcf {pre_filter_vcf}.gz")

    fs_str = "--no-frameshifts" if args.no_frameshifts else ""
    indel_str = "--no-indels" if args.no_indels else ""
    cmds.append(
        f"artic_vcf_filter {fs_str} {indel_str} --min-depth {args.min_depth} {pre_filter_vcf}.gz {args.sample}.pass.vcf {args.sample}.fail.vcf"
    )

    # 9) get the depth of coverage for each readgroup, create a coverage mask and plots, and add failed variants to the coverage mask (artic_mask must be run before bcftools consensus)
    cmds.append(
        f"artic_make_depth_mask --depth {args.min_depth} --store-rg-depths {ref} {args.sample}.primertrimmed.rg.sorted.bam {args.sample}.coverage_mask.txt"
    )

    cmds.append(
        f"artic_mask {ref} {args.sample}.coverage_mask.txt {args.sample}.fail.vcf {args.sample}.preconsensus.fasta"
    )

    post_filter_vcf_file = f"{args.sample}.pass.vcf"

    # 10) generate the consensus sequence
    cmds.append(f"bgzip -kf {post_filter_vcf_file}")
    cmds.append(f"tabix -f -p vcf {post_filter_vcf_file}.gz")

    # Normalise variants in the pass/fail VCF files
    post_normalisation_vcf_file = f"{args.sample}.normalised.vcf.gz"
    cmds.append(
        f"bcftools norm --check-ref x --fasta-ref {args.sample}.preconsensus.fasta -O z -o {post_normalisation_vcf_file} {post_filter_vcf_file}.gz"
    )
    # cmds.append(f"bgzip -kf {post_normalisation_vcf_file}")
    cmds.append(f"tabix -f -p vcf {post_normalisation_vcf_file}")
    cmds.append(
        f"bcftools consensus -f {args.sample}.preconsensus.fasta {post_normalisation_vcf_file} -m {args.sample}.coverage_mask.txt -o {args.sample}.consensus.fasta"
    )

    # 11) apply the header to the consensus sequence and run alignment against the reference sequence
    caller = "clair3"
    cmds.append(
        f"artic_fasta_header {args.sample}.consensus.fasta {args.sample} {caller}"
    )

    if args.align_consensus:
        cmds.append(
            f"mafft --6merpair --addfragments {args.sample}.consensus.fasta {ref} > {args.sample}.aligned.fasta"
        )

    # 13) setup the log file and run the pipeline commands
    log = "%s.minion.log.txt" % (args.sample)
    logfh = open(log, "w")
    for cmd in cmds:
        print(colored.green("Running: ") + cmd, file=sys.stderr)
        if not args.dry_run:
            timerStart = time.perf_counter()
            retval = os.system(cmd)
            if retval != 0:
                print(colored.red("Command failed:") + cmd, file=sys.stderr)
                raise SystemExit(20)
            timerStop = time.perf_counter()

            ## print the executed command and the runtime to the log file
            print("{}\t{}".format(cmd, timerStop - timerStart), file=logfh)

        ## if it's a dry run, print just the command
        else:
            print(cmd, file=logfh)
    logfh.close()
