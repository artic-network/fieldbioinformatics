# Written by Nick Loman (@pathogenomenick)

from clint.textui import colored
import os
import sys
import time
import requests
import hashlib
from artic.utils import read_bed_file


def check_scheme_hashes(filepath, manifest_hash):
    with open(filepath, "rb") as fh:
        data = fh.read()
        hash_sha256 = hashlib.sha256(data).hexdigest()
    if hash_sha256 != manifest_hash:
        print(
            colored.yellow(f"sha256 hash for {str(filepath)} does not match manifest"),
            file=sys.stderr,
        )
        raise SystemExit(1)


# def scheme_fetcher(
#     scheme_name: str, scheme_directory: Path, scheme_version: str, scheme_length: int
# ):

#     pass


def get_scheme(scheme_name, scheme_directory, scheme_version="1"):
    """Get and check the ARTIC primer scheme.
    When determining a version, the original behaviour (parsing the scheme_name and
    separating on /V ) is used over a specified scheme_version. If neither are
    provided, the version defaults to 1.
    If 0 is provided as version, the latest scheme will be downloaded.

    Parameters
    ----------
    scheme_name : str
        The primer scheme name
    scheme_directory : str
        The directory containing the primer scheme and reference sequence
    scheme_version : str
        The primer scheme version (optional)
    Returns
    -------
    str
        The location of the checked primer scheme
    str
        The location of the checked reference sequence
    str
        The version being used
    """
    # try getting the version from the scheme name (old behaviour)
    if scheme_name.find("/V") != -1:
        scheme_name, scheme_version = scheme_name.split("/V")

    # create the filenames and check they exist
    bed = "%s/%s/V%s/%s.scheme.bed" % (
        scheme_directory,
        scheme_name,
        scheme_version,
        scheme_name,
    )
    ref = "%s/%s/V%s/%s.reference.fasta" % (
        scheme_directory,
        scheme_name,
        scheme_version,
        scheme_name,
    )
    if os.path.exists(bed) and os.path.exists(ref):
        return bed, ref, scheme_version

    # if they don't exist, try downloading them to the current directory
    print(
        colored.yellow(
            "could not find primer scheme and reference sequence, downloading"
        ),
        file=sys.stderr,
    )

    try:
        manifest = requests.get(
            "https://raw.githubusercontent.com/artic-network/primer-schemes/master/schemes_manifest.json"
        ).json()
    except requests.exceptions.RequestException as error:
        print("Manifest Exception:", error)
        raise SystemExit(2)

    for scheme, scheme_contents in dict(manifest["schemes"]).items():
        if (
            scheme == scheme_name.lower()
            or scheme_name.lower() in scheme_contents["aliases"]
        ):
            print(
                colored.yellow(
                    f"\tfound requested scheme:\t{scheme} (using alias {scheme_name})"
                ),
                file=sys.stderr,
            )
            if scheme_version == 0:
                print(
                    colored.yellow(
                        f"Latest version for scheme {scheme} is -> {scheme_contents['latest_version']}"
                    ),
                    file=sys.stderr,
                )
                scheme_version = scheme_contents["latest_version"]
            elif scheme_version not in dict(scheme_contents["primer_urls"]).keys():
                print(
                    colored.yellow(
                        f"Requested scheme version {scheme_version} not found; using latest version ({scheme_contents['latest_version']}) instead"
                    ),
                    file=sys.stderr,
                )
                scheme_version = scheme_contents["latest_version"]
                bed = "%s/%s/V%s/%s.scheme.bed" % (
                    scheme_directory,
                    scheme_name,
                    scheme_version,
                    scheme_name,
                )
                ref = "%s/%s/V%s/%s.reference.fasta" % (
                    scheme_directory,
                    scheme_name,
                    scheme_version,
                    scheme_name,
                )

            os.makedirs(os.path.dirname(bed), exist_ok=True)
            with requests.get(scheme_contents["primer_urls"][scheme_version]) as fh:
                open(bed, "wt").write(fh.text)

            os.makedirs(os.path.dirname(ref), exist_ok=True)
            with requests.get(scheme_contents["reference_urls"][scheme_version]) as fh:
                open(ref, "wt").write(fh.text)

            check_scheme_hashes(
                bed, scheme_contents["primer_sha256_checksums"][scheme_version]
            )
            check_scheme_hashes(
                ref, scheme_contents["reference_sha256_checksums"][scheme_version]
            )

            return bed, ref, scheme_version

    print(
        colored.yellow(
            f"\tRequested scheme:\t{scheme_name} could not be found, exiting"
        ),
        file=sys.stderr,
    )
    raise SystemExit(1)


def run(parser, args):

    # check for medaka-model
    if not args.model:
        print(
            colored.red("Must specify --model for clair3 or medaka variant calling"),
        )
        raise SystemExit(1)

    # 1) check the parameters and set up the filenames
    ## find the primer scheme, reference sequence and confirm scheme version
    bed, ref, _ = get_scheme(args.scheme, args.scheme_directory, args.scheme_version)

    # ## if in strict mode, validate the primer scheme
    # if args.strict:
    #     checkScheme = "artic-tools validate_scheme %s" % (bed)
    #     print(colored.green("Running: "), checkScheme, file=sys.stderr)
    #     if os.system(checkScheme) != 0:
    #         print(colored.red("primer scheme failed strict checking"), file=sys.stderr)
    #         raise SystemExit(1)

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

    # 2) check the reference fasta has and index and create one if not
    if not os.path.exists("%s.fai" % (ref)) and args.clair3:
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
        f"align_trim {normalise_string} {bed} --primer-match-threshold {args.primer_match_threshold} --remove-incorrect-pairs --min-mapq {args.min_mapq} --report {args.sample}.alignreport.csv < {args.sample}.sorted.bam 2> {args.sample}.alignreport.er | samtools sort -T {args.sample} - -o {args.sample}.trimmed.rg.sorted.bam"
    )

    cmds.append(
        f"align_trim {normalise_string} {bed} --primer-match-threshold {args.primer_match_threshold} --min-mapq {args.min_mapq} --remove-incorrect-pairs --trim-primers --report {args.sample}.alignreport.csv < {args.sample}.sorted.bam 2> {args.sample}.alignreport.er | samtools sort -T {args.sample} - -o {args.sample}.primertrimmed.rg.sorted.bam"
    )

    cmds.append(f"samtools index {args.sample}.trimmed.rg.sorted.bam")
    cmds.append(f"samtools index {args.sample}.primertrimmed.rg.sorted.bam")

    # 6) do variant calling on each read group
    for p in pools:
        if os.path.exists("%s.%s.hdf" % (args.sample, p)):
            os.remove("%s.%s.hdf" % (args.sample, p))

        if args.clair3:
            # Use a specific model path if provided else use the default conda path
            if args.model_path:
                model_path = f"'{args.model_path}/{args.model}'"

            else:
                model_path = f"'{os.getenv('CONDA_PREFIX')}/bin/models/{args.model}'"

            # Split the BAM by read group
            for p in pools:
                cmds.append(
                    f"samtools view -b -r {p} {args.sample}.primertrimmed.rg.sorted.bam -o {args.sample}.{p}.primertrimmed.rg.sorted.bam"
                )

                cmds.append(
                    f"samtools index {args.sample}.{p}.primertrimmed.rg.sorted.bam"
                )

                cmds.append(
                    f"run_clair3.sh --chunk_size=10000 --no_phasing_for_fa --bam_fn='{args.sample}.{p}.primertrimmed.rg.sorted.bam' --ref_fn='{ref}' --output='./rg_{p}' --threads='{args.threads}' --platform='ont' --model_path={model_path} --include_all_ctgs"
                )

                cmds.append(
                    f"bgzip -dc ./rg_{p}/merge_output.vcf.gz > {args.sample}.{p}.vcf"
                )

        else:
            cmds.append(
                "medaka consensus --model %s --threads %s --chunk_len 800 --chunk_ovlp 400 --RG %s %s.trimmed.rg.sorted.bam %s.%s.hdf"
                % (args.model, args.threads, p, args.sample, args.sample, p)
            )
            if args.no_indels:
                cmds.append(
                    "medaka snp %s %s.%s.hdf %s.%s.vcf"
                    % (ref, args.sample, p, args.sample, p)
                )
            else:
                cmds.append(
                    "medaka variant %s %s.%s.hdf %s.%s.vcf"
                    % (ref, args.sample, p, args.sample, p)
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

    # 8) check and filter the VCFs
    ## if using strict, run the vcf checker to remove vars present only once in overlap regions (this replaces the original merged vcf from the previous step)
    # if args.strict:
    #     cmds.append("bgzip -f %s.merged.vcf" % (args.sample))
    #     cmds.append("tabix -p vcf %s.merged.vcf.gz" % (args.sample))
    #     cmds.append(
    #         "artic-tools check_vcf --dropPrimerVars --dropOverlapFails --vcfOut %s.merged.filtered.vcf %s.merged.vcf.gz %s 2> %s.vcfreport.txt"
    #         % (args.sample, args.sample, bed, args.sample)
    #     )
    #     cmds.append(
    #         "mv %s.merged.filtered.vcf %s.merged.vcf" % (args.sample, args.sample)
    #     )
    pre_filter_vcf = f"{args.sample}.merged.vcf"
    cmds.append(f"bgzip -kf {pre_filter_vcf}")
    cmds.append(f"tabix -f -p vcf {pre_filter_vcf}.gz")

    # if doing the medaka workflow and longshot required, do it on the merged VCF
    if not args.no_longshot:
        pre_filter_vcf = f"{args.sample}.longshot.vcf"
        cmds.append(
            f"longshot --min_mapq {args.min_mapq} -P 0 -F --max_cov 500 --no_haps --bam {args.sample}.primertrimmed.rg.sorted.bam --ref {ref} --out {pre_filter_vcf} --potential_variants {args.sample}.merged.vcf.gz"
        )
        cmds.append(f"bgzip -kf {pre_filter_vcf}")
        cmds.append(f"tabix -f -p vcf {pre_filter_vcf}.gz")

    ## filter the variants to produce PASS and FAIL lists, then index them
    fs_str = "--no-frameshifts" if args.no_frameshifts else ""
    indel_str = "--no-indels" if args.no_indels else ""
    cmds.append(
        f"artic_vcf_filter {fs_str} {indel_str} {pre_filter_vcf}.gz {args.sample}.pass.vcf {args.sample}.fail.vcf"
    )

    # 9) get the depth of coverage for each readgroup, create a coverage mask and plots, and add failed variants to the coverage mask (artic_mask must be run before bcftools consensus)
    cmds.append(
        f"artic_make_depth_mask --store-rg-depths {ref} {args.sample}.primertrimmed.rg.sorted.bam {args.sample}.coverage_mask.txt"
    )

    cmds.append(
        f"artic_mask {ref} {args.sample}.coverage_mask.txt {args.sample}.fail.vcf {args.sample}.preconsensus.fasta"
    )

    post_filter_vcf_file = f"{args.sample}.pass.vcf"

    # 10) generate the consensus sequence
    cmds.append(f"bgzip -kf {post_filter_vcf_file}")
    cmds.append(f"tabix -f -p vcf {post_filter_vcf_file}.gz")

    # Normalise variants in the pass/fail VCF files
    post_normalisation_vcf_file = f"{args.sample}.normalised.vcf"
    cmds.append(
        f"bcftools norm --check-ref x --fasta-ref {args.sample}.preconsensus.fasta -O z -o {post_normalisation_vcf_file} {post_filter_vcf_file}.gz"
    )
    cmds.append(f"bgzip -kf {post_normalisation_vcf_file}")
    cmds.append(f"tabix -f -p vcf {post_normalisation_vcf_file}.gz")
    cmds.append(
        f"bcftools consensus -f {args.sample}.preconsensus.fasta {post_normalisation_vcf_file}.gz -m {args.sample}.coverage_mask.txt -o {args.sample}.consensus.fasta"
    )

    # 11) apply the header to the consensus sequence and run alignment against the reference sequence
    fasta_header = "%s/ARTIC/medaka" % (args.sample)
    cmds.append(
        'artic_fasta_header %s.consensus.fasta "%s"' % (args.sample, fasta_header)
    )

    if args.use_muscle:
        cmds.append(
            "cat %s.consensus.fasta %s > %s.muscle.in.fasta"
            % (args.sample, ref, args.sample)
        )
        cmds.append(
            "muscle -in %s.muscle.in.fasta -out %s.muscle.out.fasta"
            % (args.sample, args.sample)
        )

    # 12) get some QC stats
    if args.strict:
        cmds.append(
            "artic_get_stats --scheme {} --align-report {}.alignreport.txt --vcf-report {}.vcfreport.txt {}".format(
                bed, args.sample, args.sample, args.sample
            )
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
