#Written by Nick Loman (@pathogenomenick)

from clint.textui import colored, puts, indent
from Bio import SeqIO
import os
import sys
import time
import requests
import hashlib
from .vcftagprimersites import read_bed_file

def get_nanopolish_header(ref):
    """Checks the reference sequence for a single fasta entry.

    Parameters
    ----------
    ref : str
        The fasta file containing the reference sequence

    Returns
    -------
    str
        A formatted header string for nanopolish
    """
    recs = list(SeqIO.parse(open(ref), "fasta"))
    if len (recs) != 1:
        print("FASTA has more than one sequence", file=sys.stderr)
        raise SystemExit(1)
    return  "%s:%d-%d" % (recs[0].id, 1, len(recs[0])+1)

def check_scheme_hashes(filepath, manifest_hash):
    with open(filepath, "rb") as fh:
        data = fh.read()
        hash_sha256 = hashlib.sha256(data).hexdigest()
    if hash_sha256 != manifest_hash:
        print(
            colored.yellow(
                f"sha256 hash for {str(filepath)} does not match manifest"
            ),
            file=sys.stderr,
        )
        raise SystemExit(1)

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
    if scheme_name.find('/V') != -1:
        scheme_name, scheme_version = scheme_name.split('/V')

    # create the filenames and check they exist
    bed = "%s/%s/V%s/%s.scheme.bed" % (scheme_directory, scheme_name, scheme_version, scheme_name)
    ref = "%s/%s/V%s/%s.reference.fasta" % (scheme_directory, scheme_name, scheme_version, scheme_name)
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
        manifest = requests.get("https://raw.githubusercontent.com/artic-network/primer-schemes/master/schemes_manifest.json").json()
    except requests.exceptions.RequestException as error:
        print("Manifest Exception:", error)
        raise SystemExit(2)

    for scheme, scheme_contents in dict(manifest["schemes"]).items():
        if scheme == scheme_name.lower() or scheme_name.lower() in scheme_contents["aliases"]:
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
                bed = "%s/%s/V%s/%s.scheme.bed" % (scheme_directory, scheme_name, scheme_version, scheme_name)
                ref = "%s/%s/V%s/%s.reference.fasta" % (scheme_directory, scheme_name, scheme_version, scheme_name)
            

            os.makedirs(os.path.dirname(bed), exist_ok=True)
            with requests.get(scheme_contents["primer_urls"][scheme_version]) as fh:
                open(bed, 'wt').write(fh.text)
            
            os.makedirs(os.path.dirname(ref), exist_ok=True)
            with requests.get(scheme_contents["reference_urls"][scheme_version]) as fh:
                open(ref, 'wt').write(fh.text)
            
            check_scheme_hashes(bed, scheme_contents["primer_sha256_checksums"][scheme_version])
            check_scheme_hashes(ref, scheme_contents["reference_sha256_checksums"][scheme_version])

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
    if args.medaka and (args.medaka_model is None):
        print(colored.red('Must specify --medaka-model if using the --medaka workflow.'))
        raise SystemExit(1)

    # 1) check the parameters and set up the filenames
    ## find the primer scheme, reference sequence and confirm scheme version
    bed, ref, _ = get_scheme(args.scheme, args.scheme_directory, args.scheme_version)

    ## if in strict mode, validate the primer scheme
    if args.strict:
        checkScheme = "artic-tools validate_scheme %s" % (bed)
        print(colored.green("Running: "), checkScheme, file=sys.stderr)
        if (os.system(checkScheme) != 0):
            print(colored.red("primer scheme failed strict checking"), file=sys.stderr)
            raise SystemExit(1)

    ## set up the read file
    if args.read_file:
        read_file = args.read_file
    else:
        read_file = "%s.fasta" % (args.sample)
    if not os.path.exists(read_file):
        print(colored.red("failed to find read-file: {}" .format(read_file)), file=sys.stderr)
        raise SystemExit(1)

    ## collect the primer pools
    pools = set([row['PoolName'] for row in read_bed_file(bed)])

    ## create a holder to keep the pipeline commands in
    cmds = []

    # 2) if using nanopolish, set up the reference header and run the nanopolish indexing
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
    if args.strict:
        strict_string = '--enforce-amplicon-span'
    else:
        strict_string = ''
    cmds.append("align_trim %s %s %s --start --remove-incorrect-pairs --report %s.alignreport.txt < %s.sorted.bam 2> %s.alignreport.er | samtools sort -T %s - -o %s.trimmed.rg.sorted.bam" % (strict_string, normalise_string, bed, args.sample, args.sample, args.sample, args.sample, args.sample))
    cmds.append("align_trim %s %s %s --remove-incorrect-pairs --report %s.alignreport.txt < %s.sorted.bam 2> %s.alignreport.er | samtools sort -T %s - -o %s.primertrimmed.rg.sorted.bam" % (strict_string, normalise_string, bed, args.sample, args.sample, args.sample, args.sample, args.sample))
    cmds.append("samtools index %s.trimmed.rg.sorted.bam" % (args.sample))
    cmds.append("samtools index %s.primertrimmed.rg.sorted.bam" % (args.sample))

    # 6) do variant calling on each read group, either using the medaka or nanopolish workflow
    if args.medaka:
        for p in pools:
            if os.path.exists("%s.%s.hdf" % (args.sample, p)):
                os.remove("%s.%s.hdf" % (args.sample, p))
            cmds.append("medaka consensus --model %s --threads %s --chunk_len 800 --chunk_ovlp 400 --RG %s %s.trimmed.rg.sorted.bam %s.%s.hdf" % (args.medaka_model, args.threads, p, args.sample, args.sample, p))
            if args.no_indels:
                cmds.append("medaka snp %s %s.%s.hdf %s.%s.vcf" % (ref, args.sample, p, args.sample, p))
            else:
                cmds.append("medaka variant %s %s.%s.hdf %s.%s.vcf" % (ref, args.sample, p, args.sample, p))
            # note: medaka annotate can't overwrite files so we need to write to a tmp file and then overwrite it ourselves
            cmds.append("medaka tools annotate --pad 25 --RG %s %s.%s.vcf %s %s.trimmed.rg.sorted.bam %s.medaka-annotate.vcf" % (p, args.sample, p, ref, args.sample, args.sample))
            cmds.append("mv %s.medaka-annotate.vcf %s.%s.vcf" % (args.sample, args.sample, p))

    else:
        if not args.skip_nanopolish:
            indexed_nanopolish_file = read_file
            if args.no_indels:
                nanopolish_extra_args = " --snps"
            else:
                nanopolish_extra_args = ""
            for p in pools:
                cmds.append("nanopolish variants --min-flanking-sequence 10 -x %s --progress -t %s --reads %s -o %s.%s.vcf -b %s.trimmed.rg.sorted.bam -g %s -w \"%s\" --ploidy 1 -m 0.15 --read-group %s %s" % (args.max_haplotypes, args.threads, indexed_nanopolish_file, args.sample, p, args.sample, ref, nanopolish_header, p, nanopolish_extra_args))

    # 7) merge the called variants for each read group and index the file
    merge_vcf_cmd = "artic_vcf_merge %s %s 2> %s.primersitereport.txt" % (args.sample, bed, args.sample)
    for p in pools:
        merge_vcf_cmd += " %s:%s.%s.vcf" % (p, args.sample, p)
    cmds.append(merge_vcf_cmd)
    cmds.append("bgzip -f %s.merged.vcf" % (args.sample))
    cmds.append("tabix -p vcf %s.merged.vcf.gz" % (args.sample))

    # 8) check and filter the VCFs
    ## if using strict, run the vcf checker to remove vars present only once in overlap regions (this replaces the original merged vcf from the previous step)
    if args.strict:
        cmds.append("artic-tools check_vcf --summaryOut {}.vcfreport.txt --vcfOut {}.merged.filtered.vcf {}.merged.vcf.gz {} 2> {}.vcfcheck.log" .format(args.sample, args.sample, args.sample, bed, args.sample))
        cmds.append("mv %s.merged.filtered.vcf %s.merged.vcf" % (args.sample, args.sample))
        cmds.append("bgzip -f %s.merged.vcf" % (args.sample))
        cmds.append("tabix -f -p vcf %s.merged.vcf.gz" % (args.sample))
    else:
        cmds.append("artic-tools check_vcf --summaryOut {}.vcfreport.txt {}.merged.vcf.gz {} 2> {}.vcfcheck.log" .format(args.sample, args.sample, bed, args.sample))

    ## set up some name holder vars for ease
    if args.medaka:
        method = 'medaka'
    else:
        method = 'nanopolish'
    vcf_file = "%s.pass.vcf" % (args.sample)

    ## filter the variants to produce PASS and FAIL lists, then index them
    if args.no_frameshifts and not args.no_indels:
        cmds.append("artic_vcf_filter --%s --min-depth %s --min-qual %s --nanopolish-qual-cov-ratio %s --no-frameshifts %s.merged.vcf.gz %s.pass.vcf %s.fail.vcf" % (method, args.min_depth, args.min_qual, args.nanopolish_qual_cov_ratio, args.sample, args.sample, args.sample))
    else:
        cmds.append("artic_vcf_filter --%s --min-depth %s --min-qual %s --nanopolish-qual-cov-ratio %s %s.merged.vcf.gz %s.pass.vcf %s.fail.vcf" % (method, args.min_depth, args.min_qual, args.nanopolish_qual_cov_ratio, args.sample, args.sample, args.sample))

    # 9) get the depth of coverage for each readgroup, create a coverage mask and plots, and add failed variants to the coverage mask (artic_mask must be run before bcftools consensus)
    cmds.append("artic_make_depth_mask --depth %s %s %s.primertrimmed.rg.sorted.bam %s.coverage_mask.txt" % (args.min_depth, ref, args.sample, args.sample))
    cmds.append("artic_mask %s %s.coverage_mask.txt %s.fail.vcf %s.preconsensus.fasta" % (ref, args.sample, args.sample, args.sample))

    # 10) generate the consensus sequence
    cmds.append("bcftools norm --check-ref x --fasta-ref %s.preconsensus.fasta -O z -o %s.gz %s" %(args.sample, vcf_file, vcf_file))
    cmds.append("tabix -p vcf %s.gz" % (vcf_file))
    cmds.append("bcftools consensus -f %s.preconsensus.fasta %s.gz -m %s.coverage_mask.txt -o %s.consensus.fasta" % (args.sample, vcf_file, args.sample, args.sample))

    # 11) apply the header to the consensus sequence and run alignment against the reference sequence
    fasta_header = "%s/ARTIC/%s" % (args.sample, method)
    cmds.append("artic_fasta_header %s.consensus.fasta \"%s\"" % (args.sample, fasta_header))
    cmds.append("cat %s.consensus.fasta %s > %s.muscle.in.fasta" % (args.sample, ref, args.sample))
    cmds.append("muscle -in %s.muscle.in.fasta -out %s.muscle.out.fasta" % (args.sample, args.sample))

    # 12) get some QC stats
    cmds.append("artic_get_stats --min-depth {} --scheme {} --align-report {}.alignreport.txt --vcf-report {}.vcfreport.txt {}" .format(args.min_depth, bed, args.sample, args.sample, args.sample))

    # 13) setup the log file and run the pipeline commands
    log = "%s.minion.log.txt" % (args.sample)
    logfh = open(log, 'w')
    for cmd in cmds:
        print(colored.green("Running: ") + cmd, file=sys.stderr)
        if not args.dry_run:
            timerStart = time.perf_counter()
            retval = os.system(cmd)
            if retval != 0:
                print(colored.red('Command failed:' ) + cmd, file=sys.stderr)
                raise SystemExit(20)
            timerStop = time.perf_counter()

            ## print the executed command and the runtime to the log file
            print("{}\t{}" .format(cmd, timerStop-timerStart), file=logfh)
        
        ## if it's a dry run, print just the command
        else:
            print(cmd, file=logfh)
    logfh.close()
