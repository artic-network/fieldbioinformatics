# Written by Nick Loman (@pathogenomenick)

from clint.textui import colored
import os
import subprocess
import sys
import time
from artic.utils import get_scheme, choose_model
from primalbedtools.scheme import Scheme


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
                f"Model '{str(args.model) if args.model else model['name']}' not found in '{model_path}', please run 'artic_get_models' to download the clair3 models from ONT"
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

        os.system(f"cp {bed} {args.sample}.primer.bed")
        os.system(f"cp {ref} {args.sample}.reference.fasta")

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
    scheme = Scheme.from_file(bed)
    pools = set([row.pool for row in scheme.bedlines] + ["unmatched"])

    ## create a holder to keep the pipeline commands in
    cmds = []

    # check the reference fasta has an index and create one if not
    if not os.path.exists("%s.fai" % (ref)):
        cmds.append("samtools faidx %s" % (ref))

    # index the ref & align with minimap
    cmds.append(
        f"minimap2 -a -x map-ont -t {args.threads} {ref} {read_file} | samtools view -bS -F 4 - | samtools sort -o {args.sample}.sorted.bam -"
    )

    cmds.append(f"samtools index {args.sample}.sorted.bam")

    # trim the alignments to the primer start sites and normalise the coverage to save time
    normalise_string = f"--normalise {args.normalise}" if args.normalise else ""

    incorrect_pairs_string = (
        "--allow-incorrect-pairs" if args.allow_mismatched_primers else ""
    )

    cmds.append(
        f"align_trim {normalise_string} {bed} --primer-match-threshold {args.primer_match_threshold} --min-mapq {args.min_mapq} {incorrect_pairs_string} --report {args.sample}.alignreport.tsv --amp-depth-report {args.sample}.amplicon_depths.tsv --samfile {args.sample}.sorted.bam -o {args.sample}.primertrimmed.rg.bam"
    )

    cmds.append(
        f"samtools sort -T {args.sample} {args.sample}.primertrimmed.rg.bam -o {args.sample}.primertrimmed.rg.sorted.bam"
    )
    cmds.append(f"rm {args.sample}.primertrimmed.rg.bam")

    cmds.append(f"samtools index {args.sample}.primertrimmed.rg.sorted.bam")

    # do variant calling on each read group
    for p in pools:
        if os.path.exists("%s.%s.hdf" % (args.sample, p)):
            os.remove("%s.%s.hdf" % (args.sample, p))

        # Split the BAM by read group
        cmds.append(
            f"samtools view -b -r {p} {args.sample}.primertrimmed.rg.sorted.bam -o {args.sample}.{p}.primertrimmed.rg.sorted.bam"
        )

        cmds.append(f"samtools index {args.sample}.{p}.primertrimmed.rg.sorted.bam")

        cmds.append(
            f"run_clair3.sh --enable_long_indel --chunk_size=10000 --haploid_sensitive --no_phasing_for_fa --bam_fn='{args.sample}.{p}.primertrimmed.rg.sorted.bam' --ref_fn='{ref}' --output='{args.sample}_rg_{p}' --threads='{args.threads}' --platform='ont' --model_path='{full_model_path}' --include_all_ctgs --enable_variant_calling_at_sequence_head_and_tail"
        )

        cmds.append(
            f"bgzip -dc {args.sample}_rg_{p}/merge_output.vcf.gz > {args.sample}.{p}.vcf"
        )

        cmds.append(f"rm {args.sample}.{p}.primertrimmed.rg.sorted.bam")

    # merge the called variants for each read group
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

    # get the depth of coverage for each readgroup, create a coverage mask and plots, and add failed variants to the coverage mask (artic_mask must be run before bcftools consensus)
    cmds.append(
        f"artic_make_depth_mask --depth {args.min_depth} {ref} {args.sample}.primertrimmed.rg.sorted.bam {args.sample}.coverage_mask.txt"
    )

    cmds.append(
        f"artic_mask {ref} {args.sample}.coverage_mask.txt {args.sample}.fail.vcf {args.sample}.preconsensus.fasta"
    )

    post_filter_vcf_file = f"{args.sample}.pass.vcf"

    # generate the consensus sequence
    cmds.append(f"bgzip -kf {post_filter_vcf_file}")
    cmds.append(f"tabix -f -p vcf {post_filter_vcf_file}.gz")

    # Normalise variants in the pass/fail VCF files
    post_normalisation_vcf_file = f"{args.sample}.normalised.vcf.gz"
    cmds.append(
        f"bcftools norm --check-ref x --fasta-ref {args.sample}.preconsensus.fasta -O z -o {post_normalisation_vcf_file} {post_filter_vcf_file}.gz"
    )
    cmds.append(f"tabix -f -p vcf {post_normalisation_vcf_file}")
    cmds.append(
        f"bcftools consensus -f {args.sample}.preconsensus.fasta {post_normalisation_vcf_file} -m {args.sample}.coverage_mask.txt -o {args.sample}.consensus.fasta"
    )

    # apply the header to the consensus sequence and run alignment against the reference sequence
    caller = "clair3"
    linearise_fasta = "--linearise-fasta" if args.linearise_fasta else ""
    cmds.append(
        f"artic_fasta_header {linearise_fasta} {args.sample}.consensus.fasta {args.sample}"
    )

    if args.align_consensus:
        cmds.append(
            f"mafft --6merpair --addfragments {args.sample}.consensus.fasta {ref} > {args.sample}.aligned.fasta"
        )

    # define anticipated error codes and messages
    general_error_codes = {
        137: {
            "message": "Process ran out of memory, please try running on a machine with more memory available.",
        },
        143: {
            "message": "Process was terminated, this is likely due to a user or system interrupt.",
        },
        127: {
            "message": "A required tool was not found, please ensure all dependencies are installed and available in your PATH.",
        },
    }

    anticipated_tool_specific_error_codes = {
        "align_trim": [
            {
                "input_exit_code": 1,
                "message": "No reads aligned to the reference sequence, this could be due to an incorrect reference / primer scheme, or very low depth input data. Please check your input files.",
                "output_exit_code": 2,
            }
        ]
    }

    # setup the log file and run the pipeline commands
    log = "%s.minion.log.txt" % (args.sample)
    logfh = open(log, "w")
    for cmd in cmds:
        print(colored.green("Running: ") + cmd, file=sys.stderr)
        if not args.dry_run:
            timerStart = time.perf_counter()
            subprocess_return = subprocess.run(cmd, shell=True, capture_output=True)
            if subprocess_return.returncode != 0:
                ## Check for general anticipated errors
                if subprocess_return.returncode in general_error_codes:
                    print(
                        colored.yellow(
                            general_error_codes[subprocess_return.returncode]["message"]
                        ),
                        file=sys.stderr,
                    )
                    print(colored.red("Command failed: ") + cmd, file=sys.stderr)
                    print(
                        colored.red("Command stderr: ")
                        + subprocess_return.stderr.decode()
                    )
                    raise SystemExit(subprocess_return.returncode)

                ## check for anticipated tool-specific errors
                cmd_parts = []
                for cmd_part in cmd.split(" "):
                    if "-" not in cmd_part:
                        cmd_parts.append(cmd_part)
                        continue

                    break

                if cmd_parts:
                    base_cmd = " ".join(cmd_parts)

                    if base_cmd in anticipated_tool_specific_error_codes:
                        for error_case in anticipated_tool_specific_error_codes[
                            base_cmd
                        ]:
                            if (
                                subprocess_return.returncode
                                == error_case["input_exit_code"]
                            ):
                                print(
                                    colored.yellow(error_case["message"]),
                                    file=sys.stderr,
                                )
                                print(
                                    colored.red("Command failed: ") + cmd,
                                    file=sys.stderr,
                                )
                                print(
                                    colored.red("Command stderr: ")
                                    + subprocess_return.stderr.decode()
                                )
                                raise SystemExit(error_case["output_exit_code"])

                print(
                    colored.red("Unexpected command failure: ") + cmd, file=sys.stderr
                )
                print(
                    colored.red("Command stderr: ") + subprocess_return.stderr.decode(),
                    file=sys.stderr,
                )

                raise SystemExit(subprocess_return.returncode)
            timerStop = time.perf_counter()

            ## print the executed command and the runtime to the log file
            print("{}\t{}".format(cmd, timerStop - timerStart), file=logfh)

        ## if it's a dry run, print just the command
        else:
            print(cmd, file=logfh)
    logfh.close()
