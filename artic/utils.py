#!/usr/bin/env python

import sys
import os
import requests
import hashlib
import re
import pandas as pd
from clint.textui import colored
import subprocess
from Bio import SeqIO
import gzip
import csv
import json


class clair3_manifest:

    def __init__(self):
        self.models = [
            {
                "name": "r1041_e82_260bps_fast_g632",
                "model_fname": "r1041_e82_260bps_fast_g632.tar.gz",
                "model_url": "https://cdn.oxfordnanoportal.com/software/analysis/models/clair3/r1041_e82_260bps_fast_g632.tar.gz",
                "rerio": True,
            },
            {
                "name": "r1041_e82_400bps_fast_g632",
                "model_fname": "r1041_e82_400bps_fast_g632.tar.gz",
                "model_url": "https://cdn.oxfordnanoportal.com/software/analysis/models/clair3/r1041_e82_400bps_fast_g632.tar.gz",
                "rerio": True,
            },
            {
                "name": "r1041_e82_400bps_sup_g615",
                "model_fname": "r1041_e82_400bps_sup_g615.tar.gz",
                "model_url": "https://cdn.oxfordnanoportal.com/software/analysis/models/clair3/r1041_e82_400bps_sup_g615.tar.gz",
                "rerio": True,
            },
            {
                "name": "r1041_e82_260bps_hac_g632",
                "model_fname": "r1041_e82_260bps_hac_g632.tar.gz",
                "model_url": "https://cdn.oxfordnanoportal.com/software/analysis/models/clair3/r1041_e82_260bps_hac_g632.tar.gz",
                "rerio": True,
            },
            {
                "name": "r1041_e82_400bps_hac_g615",
                "model_fname": "r1041_e82_400bps_hac_g615.tar.gz",
                "model_url": "https://cdn.oxfordnanoportal.com/software/analysis/models/clair3/r1041_e82_400bps_hac_g615.tar.gz",
                "rerio": True,
            },
            {
                "name": "r1041_e82_400bps_sup_v400",
                "model_fname": "r1041_e82_400bps_sup_v400.tar.gz",
                "model_url": "https://cdn.oxfordnanoportal.com/software/analysis/models/clair3/r1041_e82_400bps_sup_v400.tar.gz",
                "rerio": True,
            },
            {
                "name": "r1041_e82_260bps_hac_v400",
                "model_fname": "r1041_e82_260bps_hac_v400.tar.gz",
                "model_url": "https://cdn.oxfordnanoportal.com/software/analysis/models/clair3/r1041_e82_260bps_hac_v400.tar.gz",
                "rerio": True,
            },
            {
                "name": "r1041_e82_400bps_hac_g632",
                "model_fname": "r1041_e82_400bps_hac_g632.tar.gz",
                "model_url": "https://cdn.oxfordnanoportal.com/software/analysis/models/clair3/r1041_e82_400bps_hac_g632.tar.gz",
                "rerio": True,
            },
            {
                "name": "r1041_e82_400bps_sup_v410",
                "model_fname": "r1041_e82_400bps_sup_v410.tar.gz",
                "model_url": "https://cdn.oxfordnanoportal.com/software/analysis/models/clair3/r1041_e82_400bps_sup_v410.tar.gz",
                "rerio": True,
            },
            {
                "name": "r1041_e82_260bps_hac_v410",
                "model_fname": "r1041_e82_260bps_hac_v410.tar.gz",
                "model_url": "https://cdn.oxfordnanoportal.com/software/analysis/models/clair3/r1041_e82_260bps_hac_v410.tar.gz",
                "rerio": True,
            },
            {
                "name": "r1041_e82_400bps_hac_v400",
                "model_fname": "r1041_e82_400bps_hac_v400.tar.gz",
                "model_url": "https://cdn.oxfordnanoportal.com/software/analysis/models/clair3/r1041_e82_400bps_hac_v400.tar.gz",
                "rerio": True,
            },
            {
                "name": "r1041_e82_400bps_sup_v420",
                "model_fname": "r1041_e82_400bps_sup_v420.tar.gz",
                "model_url": "https://cdn.oxfordnanoportal.com/software/analysis/models/clair3/r1041_e82_400bps_sup_v420.tar.gz",
                "rerio": True,
            },
            {
                "name": "r1041_e82_260bps_sup_g632",
                "model_fname": "r1041_e82_260bps_sup_g632.tar.gz",
                "model_url": "https://cdn.oxfordnanoportal.com/software/analysis/models/clair3/r1041_e82_260bps_sup_g632.tar.gz",
                "rerio": True,
            },
            {
                "name": "r1041_e82_400bps_hac_v410",
                "model_fname": "r1041_e82_400bps_hac_v410.tar.gz",
                "model_url": "https://cdn.oxfordnanoportal.com/software/analysis/models/clair3/r1041_e82_400bps_hac_v410.tar.gz",
                "rerio": True,
            },
            {
                "name": "r1041_e82_400bps_sup_v430",
                "model_fname": "r1041_e82_400bps_sup_v430.tar.gz",
                "model_url": "https://cdn.oxfordnanoportal.com/software/analysis/models/clair3/r1041_e82_400bps_sup_v430.tar.gz",
                "rerio": True,
            },
            {
                "name": "r1041_e82_260bps_sup_v400",
                "model_fname": "r1041_e82_260bps_sup_v400.tar.gz",
                "model_url": "https://cdn.oxfordnanoportal.com/software/analysis/models/clair3/r1041_e82_260bps_sup_v400.tar.gz",
                "rerio": True,
            },
            {
                "name": "r1041_e82_400bps_hac_v420",
                "model_fname": "r1041_e82_400bps_hac_v420.tar.gz",
                "model_url": "https://cdn.oxfordnanoportal.com/software/analysis/models/clair3/r1041_e82_400bps_hac_v420.tar.gz",
                "rerio": True,
            },
            {
                "name": "r1041_e82_400bps_sup_v500",
                "model_fname": "r1041_e82_400bps_sup_v500.tar.gz",
                "model_url": "https://cdn.oxfordnanoportal.com/software/analysis/models/clair3/r1041_e82_400bps_sup_v500.tar.gz",
                "rerio": True,
            },
            {
                "name": "r1041_e82_260bps_sup_v410",
                "model_fname": "r1041_e82_260bps_sup_v410.tar.gz",
                "model_url": "https://cdn.oxfordnanoportal.com/software/analysis/models/clair3/r1041_e82_260bps_sup_v410.tar.gz",
                "rerio": True,
            },
            {
                "name": "r1041_e82_400bps_hac_v430",
                "model_fname": "r1041_e82_400bps_hac_v430.tar.gz",
                "model_url": "https://cdn.oxfordnanoportal.com/software/analysis/models/clair3/r1041_e82_400bps_hac_v430.tar.gz",
                "rerio": True,
            },
            {
                "name": "r104_e81_hac_g5015",
                "model_fname": "r104_e81_hac_g5015.tar.gz",
                "model_url": "https://cdn.oxfordnanoportal.com/software/analysis/models/clair3/r104_e81_hac_g5015.tar.gz",
                "rerio": True,
            },
            {
                "name": "r1041_e82_400bps_hac_v500",
                "model_fname": "r1041_e82_400bps_hac_v500.tar.gz",
                "model_url": "https://cdn.oxfordnanoportal.com/software/analysis/models/clair3/r1041_e82_400bps_hac_v500.tar.gz",
                "rerio": True,
            },
            {
                "name": "r104_e81_sup_g5015",
                "model_fname": "r104_e81_sup_g5015.tar.gz",
                "model_url": "https://cdn.oxfordnanoportal.com/software/analysis/models/clair3/r104_e81_sup_g5015.tar.gz",
                "rerio": True,
            },
            {
                "name": "r941_prom_sup_g5014",
                "model_fname": "NA",
                "model_url": "NA",
                "rerio": False,
            },
            {
                "name": "r941_prom_hac_g360+g422",
                "model_fname": "NA",
                "model_url": "NA",
                "rerio": False,
            },
        ]


def getPrimerDirection(primerID):
    """Infer the primer direction based on it's ID containing LEFT/RIGHT

    Parameters
    ----------
    primerID : string
        The primer ID from the 4th field of the primer scheme
    """
    if "LEFT" in primerID:
        return "+"
    elif "RIGHT" in primerID:
        return "-"
    else:
        print("LEFT/RIGHT must be specified in Primer ID", file=sys.stderr)
        raise SystemExit(1)


def merge_sites(canonical, alt):
    """Merges a canonical primer site with an alt site, producing an interval that encompasses both

    Parameters
    ----------
    canonical : dict
        The canonical primer site, provided as a dictionary of the bed file row
    alt : dict
        The alt primer site, provided as a dictionary of the bed file row

    Returns
    -------
    dict
        A dictionary of the merged site, where the dict represents a bed file row
    """
    # base the merged site on the canonical
    mergedSite = canonical

    # check the both the canonical and alt are the same direction
    if canonical["direction"] != alt["direction"]:
        print(
            "could not merge alt with different orientation to canonical",
            file=sys.stderr,
        )
        raise SystemExit(1)

    # merge the start/ends of the alt with the canonical to get the largest window possible
    if alt["start"] < canonical["start"]:
        mergedSite["start"] = alt["start"]
    if alt["end"] > canonical["end"]:
        mergedSite["end"] = alt["end"]
    return mergedSite


def identify_bed_file(bed_file):
    """Identifies the version of the primer ID format used in the bed file

    Parameters
    ----------
    bed_file : str
        The bed file to identify the primer ID format

    Returns
    -------
    int
        The version of the primer ID format used in the bed file
    """

    V1_pattern = re.compile(
        r"[a-zA-Z0-9\-]+_[0-9]+_(LEFT|RIGHT)(_ALT[0-9]*|_alt[0-9]*)*"
    )

    V2_pattern = re.compile(r"[a-zA-Z0-9\-]+_[0-9]+_(LEFT|RIGHT)_[0-9]+")

    version = False

    with open(bed_file, "r") as bed_fh:
        bed = bed_fh.readlines()

    for line in bed:
        if line.startswith("#"):
            continue

        splits = line.strip().split("\t")

        if len(splits) < 6:
            print(
                colored.red(
                    f"Invalid bed file format, only found {len(splits)} columns. For valid formats please see https://github.com/quick-lab/primerschemes"
                )
            )
            raise SystemExit(1)

        if not V1_pattern.match(splits[3]) and not V2_pattern.match(splits[3]):
            print(
                colored.red(
                    f"Invalid primer ID format, {splits[3]} does not match the expected format"
                )
            )
            raise SystemExit(1)

        if V2_pattern.match(splits[3]):
            if len(splits) < 7:
                print(
                    colored.red(
                        f"Invalid bed file format, only found {len(splits)} columns. For valid formats please see https://github.com/ChrisgKent/primal-page"
                    )
                )

            if not version:
                version = 3

            if version != 3:
                print(
                    colored.red(
                        "Scheme BED does not appear to be a consistent scheme version, for primer scheme formats please see https://github.com/ChrisgKent/primal-page"
                    )
                )
                raise SystemExit(1)

        elif V1_pattern.match(splits[3]):

            if len(splits) == 7:
                if not version:
                    version = 2
            elif len(splits) == 6:
                if not version:
                    version = 1
            else:
                print(
                    colored.red(
                        f"Invalid bed file format, found {len(splits)} columns with V1 primer names. For valid formats please see https://github.com/ChrisgKent/primal-page"
                    )
                )

            if version == 3:
                print(
                    colored.red(
                        "Scheme BED mixed primer ID formats, please ensure your BED file is consistent"
                    )
                )
                raise SystemExit(1)

        return version


def read_bed_file(fn):
    """Parses a V2/V3 bed file and collapses primers into canonical primer sites

    Parameters
    ----------
    fn : str
        The bedfile to parse

    Returns
    -------
    list
        A list of dictionaries, where each dictionary contains a row of the parsed bedfile.
        The available dictionary keys are - Primer_ID, direction, start, end
    """

    # read the primer scheme into a pandas dataframe and run type, length and null checks
    version = identify_bed_file(fn)

    if version in (1, 2):
        return read_bed_file_legacy(fn)

    primers = pd.read_csv(
        fn,
        sep="\t",
        comment="#",
        header=None,
        names=["chrom", "start", "end", "Primer_ID", "PoolName", "direction"],
        dtype={
            "chrom": str,
            "start": int,
            "end": int,
            "Primer_ID": str,
            "PoolName": str,
        },
        usecols=(0, 1, 2, 3, 4, 5),
        skiprows=0,
    )
    if len(primers.index) < 1:
        print("primer scheme file is empty", file=sys.stderr)
        raise SystemExit(1)
    if primers.isnull().sum().sum():
        print("malformed primer scheme file", file=sys.stderr)
        raise SystemExit(1)

    canonical_primers = {}
    for _, row in primers.iterrows():
        scheme_name, primer_id, direction, primer_n = row["Primer_ID"].split("_")

        if (row["chrom"], primer_id, direction) not in canonical_primers:
            canonical_primers[(row["chrom"], primer_id, direction)] = row.to_dict()
            continue

        canonical_primers[(row["chrom"], primer_id, direction)] = merge_sites(
            canonical_primers[(row["chrom"], primer_id, direction)], row
        )

    # return the bedFile as a list
    return [value for value in canonical_primers.values()]


def read_bed_file_legacy(fn):
    """Parses a bed file and collapses alts into canonical primer sites

    Parameters
    ----------
    fn : str
        The bedfile to parse

    Returns
    -------
    list
        A list of dictionaries, where each dictionary contains a row of the parsed bedfile.
        The available dictionary keys are - Primer_ID, direction, start, end
    """

    # read the primer scheme into a pandas dataframe and run type, length and null checks
    primers = pd.read_csv(
        fn,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "Primer_ID", "PoolName", "direction"],
        dtype={
            "chrom": str,
            "start": int,
            "end": int,
            "Primer_ID": str,
            "PoolName": str,
        },
        usecols=(0, 1, 2, 3, 4, 5),
        skiprows=0,
    )
    if len(primers.index) < 1:
        print("primer scheme file is empty", file=sys.stderr)
        raise SystemExit(1)
    if primers.isnull().sum().sum():
        print("malformed primer scheme file", file=sys.stderr)
        raise SystemExit(1)

    # separate alt primers into a new dataframe
    altFilter = primers["Primer_ID"].str.contains("_alt")
    alts = pd.DataFrame(
        columns=("chrom", "start", "end", "Primer_ID", "PoolName", "direction")
    )
    alts = pd.concat([alts, primers[altFilter]])
    primers = primers.drop(primers[altFilter].index.values)

    # convert the primers dataframe to dictionary, indexed by Primer_ID
    #  - verify_integrity is used to prevent duplicate Primer_IDs being processed
    bedFile = primers.set_index(
        "Primer_ID", drop=False, verify_integrity=True
    ).T.to_dict()

    # if there were no alts, return the bedfile as a list of dicts
    if len(alts.index) == 0:
        return list(bedFile.values())

    # merge alts
    for _, row in alts.iterrows():
        primerID = row["Primer_ID"].split("_alt")[0]

        # check the bedFile if another version of this primer exists
        if primerID not in bedFile:

            # add to the bed file and continue
            bedFile[primerID] = row.to_dict()
            continue

        # otherwise, we've got a primer ID we've already seen so merge the alt
        mergedSite = merge_sites(bedFile[primerID], row)

        # update the bedFile
        bedFile[primerID] = mergedSite

    # return the bedFile as a list
    return [value for value in bedFile.values()]


def overlaps(coords, pos):
    for v in coords:
        if pos >= v["start"] and pos <= v["end"]:
            return v
    return False


def check_hash(filepath, manifest_hash):
    with open(filepath, "rb") as fh:
        data = fh.read()
        hash_md5 = hashlib.md5(data).hexdigest()
    if hash_md5 != manifest_hash:
        print(
            colored.yellow(
                f"md5 hash for {str(filepath)} does not match manifest, if this scheme has been downloaded previously, please ensure the file is not corrupted or has been tampered with. If this is a new download, please raise this as an issue here: https://github.com/quick-lab/primerschemes/issues"
            ),
            file=sys.stderr,
        )
        raise SystemExit(1)


def get_scheme(
    scheme_name: str,
    scheme_version: str,
    scheme_directory: str,
    scheme_length: int = False,
    read_file: str = False,
):
    """Get the primer scheme and reference fasta file from the manifest

    Args:
        scheme_name (str): Name of the scheme
        scheme_version (str): Version of the scheme
        scheme_directory (str): Directory where schemes are stored
        scheme_length (int, optional): Length of the scheme. Defaults to False.

    Returns:
        str: Path to the primer bed file
        str: Path to the reference fasta file
        str: Version of the scheme
    """

    try:
        response = requests.get(
            "https://raw.githubusercontent.com/quick-lab/primerschemes/main/index.json"
        )
        manifest = response.json()
        os.makedirs(scheme_directory, exist_ok=True)
        with open(os.path.join(scheme_directory, "index.json"), "w") as manifest_file:
            json.dump(manifest, manifest_file)

    except Exception:
        os.makedirs(scheme_directory, exist_ok=True)
        index_path = os.path.join(scheme_directory, "index.json")
        if os.path.exists(index_path):
            print(
                colored.yellow(
                    f"Could not get a remote copy of manifest but found a local one in {scheme_directory}, using this instead. Be warned this may be out of date"
                )
            )
            with open(index_path, "r") as manifest_file:
                manifest = json.load(manifest_file)
        else:
            print(
                colored.red(
                    f"Failed to fetch manifest and no local copy found in {scheme_directory}, please check your internet connection and try again"
                )
            )
            raise SystemExit(1)

    try:
        response = requests.get(
            "https://raw.githubusercontent.com/quick-lab/primerschemes/main/aliases.json"
        )
        aliases = response.json()
        os.makedirs(scheme_directory, exist_ok=True)
        with open(os.path.join(scheme_directory, "aliases.json"), "w") as aliases_file:
            json.dump(aliases, aliases_file)

    except Exception:
        os.makedirs(scheme_directory, exist_ok=True)
        aliases_path = os.path.join(scheme_directory, "aliases.json")
        if os.path.exists(aliases_path):
            print(
                colored.yellow(
                    f"Could not get a remote copy of manifest but found a local one in {scheme_directory}, using this instead. Be warned this may be out of date"
                )
            )
            with open(aliases_path, "r") as aliases_file:
                aliases = json.load(aliases_file)
        else:
            print(
                colored.red(
                    f"Failed to fetch manifest and no local copy found in {scheme_directory}, please check your internet connection and try again"
                )
            )
            raise SystemExit(1)

    scheme_name = scheme_name.lower().rstrip()

    if not manifest["primerschemes"].get(scheme_name):
        if not aliases.get(scheme_name):
            print(
                colored.red(
                    f"Requested scheme: {scheme_name} could not be found, please check https://github.com/quick-lab/primerschemes for available schemes or provide a scheme bed and fasta reference directly using --bed and --ref"
                )
            )
            raise SystemExit(1)

        scheme_name = aliases[scheme_name]

        if not manifest["primerschemes"].get(scheme_name):
            print(
                colored.red(
                    "Scheme name alias does not exist in manifest, this should never happen, please create an issue on https://github.com/quick-lab/primerschemes/issues or https://github.com/artic-network/fieldbioinformatics/issues if you see this message"
                )
            )
            raise SystemExit(1)

    scheme = manifest["primerschemes"][scheme_name]

    if len(scheme.keys()) > 1:
        if not scheme_length:
            print(
                colored.red(
                    f"Multiple lengths of the scheme {scheme_name} found, please provide a scheme length using --scheme-length, available lengths are: {', '.join(scheme.keys())}"
                )
            )
            raise SystemExit(1)

    else:
        scheme_length = list(scheme.keys())[0]

    if not scheme.get(scheme_length):
        print(
            colored.red(
                f"Provided length: {scheme_length} of the scheme {scheme_name} not found, please provide one of the following lengths using --scheme-length: {', '.join(scheme.keys())}"
            )
        )
        raise SystemExit(1)

    scheme_version = scheme_version.lower().rstrip()

    version_pattern = re.compile(r"v\d+\.\d+\.\d+")

    if not version_pattern.match(scheme_version):
        print(
            colored.red(
                "Invalid scheme version format, please provide a version in the format 'vX.X.X', e.g. v1.0.0"
            )
        )
        raise SystemExit(1)

    if not scheme[scheme_length].get(scheme_version):
        print(
            colored.red(
                f"Requested scheme version: {scheme_version} not found, available versions are: {', '.join(scheme[scheme_length].keys())}"
            )
        )
        raise SystemExit(1)

    if scheme[scheme_length][scheme_version].get("refselect"):
        if not read_file:
            print(
                colored.red(
                    f"Reference selection is available for this scheme but reads were not provided. Either provide a read file using --read-file or specify the reference to use by providing the same scheme name with the appropriate suffix, choices are: {', '.join(str(x) for x in scheme[scheme_length].keys() if '-' in x)}"
                )
            )
            raise SystemExit(1)

        print(
            colored.yellow(
                f"Reference selection is available for scheme {scheme_name}, deciding which reference to use based on your reads. If you would prefer to specify the reference to use, provide the same scheme name with the appropriate suffix, choices are: {', '.join(str(x) for x in scheme[scheme_length].keys() if '-' in x)}"
            )
        )

        if len(scheme[scheme_length][scheme_version]["refselect"].keys()) > 1:
            print(
                colored.red(
                    f"Multiple reference selection options found for scheme {scheme_name}, fieldbioinformatics does not currently support multi pathogen schemes"
                )
            )
            raise SystemExit(1)

        msa_data = next(
            iter(scheme[scheme_length][scheme_version]["refselect"].values())
        )

        scheme_select_dir = os.path.join(
            scheme_directory, scheme_name, scheme_length, scheme_version
        )

        suffix = pick_best_ref(
            multi_ref_url=msa_data["url"],
            multi_ref_md5=msa_data["md5"],
            read_file=read_file,
            n_reads=10000,
            scheme_path=scheme_select_dir,
            mm2_threads=4,
        )

        scheme_version = f"{scheme_version}-{suffix}"

        print(colored.yellow(f"Selected reference suffix: {suffix}"))

    scheme = scheme[scheme_length][scheme_version]

    scheme_path = os.path.join(
        scheme_directory, scheme_name, scheme_length, scheme_version
    )

    if not os.path.exists(scheme_path):
        os.makedirs(scheme_path, exist_ok=True)

    bed_url = scheme["primer_bed_url"]

    bed_name = bed_url.split("/")[-1]

    bed_path = os.path.join(scheme_path, bed_name)

    if not os.path.exists(bed_path):
        try:
            response = requests.get(bed_url)
        except requests.exceptions.RequestException as error:
            print(
                colored.red(
                    f"Failed to fetch primer bed file due to Exception: {error}"
                )
            )
            raise SystemExit(1)

        with open(bed_path, "w") as bed_file:
            bed_file.write(response.text)

    check_hash(bed_path, scheme["primer_bed_md5"])

    ref_url = scheme["reference_fasta_url"]

    ref_name = ref_url.split("/")[-1]

    ref_path = os.path.join(scheme_path, ref_name)

    if not os.path.exists(ref_path):
        try:
            response = requests.get(ref_url)
        except requests.exceptions.RequestException as error:
            print(
                colored.red(
                    f"Failed to fetch reference fasta file with Exception: {error}"
                )
            )
            raise SystemExit(1)

        if response.status_code != 200:
            print(
                colored.red(
                    f"Failed to fetch reference fasta file with status code: {response.status_code}"
                )
            )
            raise SystemExit(1)

        with open(ref_path, "w") as ref_file:
            ref_file.write(response.text)

    check_hash(ref_path, scheme["reference_fasta_md5"])

    return bed_path, ref_path, scheme_version


def get_scheme_legacy(scheme_name, scheme_directory, scheme_version="1"):
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

            # check_scheme_hashes(
            #     bed, scheme_contents["primer_sha256_checksums"][scheme_version]
            # )
            # check_scheme_hashes(
            #     ref, scheme_contents["reference_sha256_checksums"][scheme_version]
            # )

            return bed, ref, scheme_version

    print(
        colored.yellow(
            f"\tRequested scheme:\t{scheme_name} could not be found, exiting"
        ),
        file=sys.stderr,
    )
    raise SystemExit(1)


def choose_model(read_file: str) -> dict:
    """
    Choose the appropriate clair3 model based on the `basecall_model_version_id` field in the read header (if it exists)

    Args:
        read_file (str): Path to the fastq file

    Returns:
        dict: The chosen clair3 model as a dictionary
    """

    models_class = clair3_manifest()
    models = models_class.models

    if read_file.endswith(".gz"):
        with gzip.open(
            read_file,
            "rt",
        ) as handle:
            reads = SeqIO.parse(handle, "fastq")

            read = next(reads)
    else:
        reads = SeqIO.parse(read_file, "fastq")
        read = next(reads)

    split_description = read.description.split()
    split_description = [x.split("=") for x in split_description if "=" in x]

    tags = {x[0]: x[1] for x in split_description}

    if "basecall_model_version_id" not in tags:
        print(
            colored.red(
                "Provided fastq does not contain basecall_model_version_id in the read header so clair3 model cannot be chosen automatically, please provide an appropriate model with the --model parameter",
            ),
            file=sys.stderr,
        )
        sys.exit(6)

    if len(tags["basecall_model_version_id"].split("_")) == 4:
        molecule, pore_type, kit_id, model = tags["basecall_model_version_id"].split(
            "_"
        )

        if "@" in model:
            preset, basecaller_version = model.split("@")
            basecaller_version = basecaller_version.replace(".", "")

        pore_type = pore_type.replace(".", "")
        kit_id = kit_id.replace(".", "")

        possible_models = [x for x in models if pore_type in x["name"]]

        if not possible_models:
            print(
                colored.red(
                    f"No model found for basecall model id {tags['basecall_model_version_id']}, please provide a model with the --model parameter",
                ),
                file=sys.stderr,
            )
            sys.exit(6)

        if len(possible_models) == 1:
            print(
                colored.yellow(
                    f"Model {possible_models[0]['name']} chosen for basecall model id {tags['basecall_model_version_id']}, if this is incorrect please provide the correct model with the --model parameter"
                ),
                file=sys.stderr,
            )
            return possible_models[0]

        try:
            possible_models = [x for x in possible_models if preset in x["name"]]
        except NameError:
            possible_models = [x for x in possible_models if model in x["name"]]

        if not possible_models:
            print(
                colored.red(
                    f"No model found for basecall model id {tags['basecall_model_version_id']}, please provide a model with the --model parameter",
                    file=sys.stderr,
                )
            )
            sys.exit(6)

        if len(possible_models) == 1:
            print(
                colored.yellow(
                    f"Model {possible_models[0]['name']} chosen for basecall model id {tags['basecall_model_version_id']}, if this is incorrect please provide the correct model with the --model parameter"
                ),
                file=sys.stderr,
            )
            return possible_models[0]

    elif len(tags["basecall_model_version_id"].split("_")) == 5:
        molecule, pore_type, kit_id, speed, model = tags[
            "basecall_model_version_id"
        ].split("_")

        if "@" in model:
            preset, basecaller_version = model.split("@")
            basecaller_version = basecaller_version.replace(".", "")

        pore_type = pore_type.replace(".", "")
        kit_id = kit_id.replace(".", "")

        search_string = f"{pore_type}_{kit_id}_{speed}"

        possible_models = [x for x in models if search_string in x["name"]]

    if not possible_models:
        print(
            colored.red(
                f"No model found for basecall model id {tags['basecall_model_version_id']}, please provide a model with the --model parameter",
                file=sys.stderr,
            )
        )
        sys.exit(6)

    if len(possible_models) == 1:
        print(
            colored.yellow(
                f"Model {possible_models[0]['name']} chosen for basecall model id {tags['basecall_model_version_id']}, if this is incorrect please provide the correct model with the --model parameter"
            )
        )
        return possible_models[0]

    try:
        possible_models = [x for x in possible_models if preset in x["name"]]
    except NameError:
        possible_models = [x for x in possible_models if model in x["name"]]

    if not possible_models:
        print(
            colored.red(
                f"No model found for basecall model id {tags['basecall_model_version_id']}, please provide a model with the --model parameter",
                file=sys.stderr,
            )
        )
        sys.exit(6)

    if len(possible_models) == 1:
        print(
            colored.yellow(
                f"Model {possible_models[0]['name']} chosen for basecall model id {tags['basecall_model_version_id']}, if this is incorrect please provide the correct model with the --model parameter"
            )
        )
        return possible_models[0]

    try:
        for x in possible_models:
            if basecaller_version in x["name"]:
                print(
                    colored.yellow(
                        f"Model {x['name']} chosen for basecall model id {tags['basecall_model_version_id']}, if this is incorrect please provide the correct model with the --model parameter"
                    ),
                    file=sys.stderr,
                )
                return x

    except NameError:
        pass

    print(
        colored.red(
            f"Multiple potential models found, please provide the appropriate model from the following with the '--model' parameter:\n{' '.join([str(x['name']) for x in possible_models])}"
        ),
        file=sys.stderr,
    )

    sys.exit(6)


def pick_best_ref(
    multi_ref_url: str,
    multi_ref_md5: str,
    read_file: str,
    n_reads: int,
    scheme_path: str,
    mm2_threads: int,
):
    """Pick the best reference from a multi-reference alignment

    Args:
        multi_ref_url (str): URL to the multi-reference alignment from quick-lab/primerschemes
        multi_ref_md5 (str): MD5 hash of the multi-reference alignment
        read_file (str): Path to the read file to test against the references
        n_reads (int): How many reads to sample from the read file for testing (default: 10000)
        scheme_path (str): Path to the scheme directory
        mm2_threads (int): Number of threads to use when aligning using minimap2

    Raises:
        SystemExit: If the ref selection fasta file cannot be fetched
        SystemExit: If the ref selection fasta file get returns a status code other than 200

    Returns:
        str: Primer scheme suffix of the most appropriate reference for the reads
    """

    ref_selection_path = os.path.join(scheme_path, "refselect.fasta")

    if not os.path.exists(scheme_path):
        os.makedirs(scheme_path, exist_ok=True)

    if not os.path.exists(ref_selection_path):
        try:
            response = requests.get(multi_ref_url)
        except requests.exceptions.RequestException as error:
            print(colored.red(f"Failed to ref selection fasta with Exception: {error}"))
            raise SystemExit(1)

        if response.status_code != 200:
            print(
                colored.red(
                    f"Failed to fetch ref selection fasta with status code: {response.status_code}"
                )
            )
            raise SystemExit(1)

        with open(ref_selection_path, "w") as ref_select_fasta:
            ref_select_fasta.write(response.text)

    check_hash(ref_selection_path, multi_ref_md5)

    flat_ref_path = os.path.join(scheme_path, "flattened_references.fasta")

    possible_references = {}

    with open(ref_selection_path, "r") as ref_selection_fh:

        for reference in SeqIO.parse(ref_selection_fh, "fasta"):
            suffix = reference.description.split()[1]
            possible_references[reference.id] = suffix

            # Flatten out the alignment into flat fasta reference
            flattened = str(reference.seq).replace("-", "")

            with open(flat_ref_path, "a") as flat_ref_fh:
                flat_ref_fh.write(f">{reference.description}\n{flattened}\n")

    # fq_it = mappy.fastx_read(fastq_path)
    seqtk = subprocess.run(
        ["seqtk", "sample", str(read_file), str(n_reads)], stdout=subprocess.PIPE
    )

    result = subprocess.run(
        [
            "minimap2",
            "-x",
            "map-ont",
            "-t",
            str(mm2_threads),
            str(flat_ref_path),
            "-",
        ],
        input=seqtk.stdout,
        stdout=subprocess.PIPE,
    )

    reader = csv.DictReader(
        result.stdout.decode("utf-8").split("\n"),
        delimiter="\t",
        fieldnames=[
            "query_name",
            "query_len",
            "query_start",
            "query_end",
            "strand",
            "ctg_name",
            "ctg_len",
            "ref_start",
            "ref_end",
            "n_matches",
            "alignment_len",
            "mapq",
        ],
    )

    read_results = {}

    for alignment in reader:

        if not alignment:
            continue

        identity = int(alignment["n_matches"]) / int(alignment["alignment_len"])

        read_results.setdefault(alignment["query_name"], {})
        read_results[alignment["query_name"]].setdefault(alignment["ctg_name"], 0)

        if identity > read_results[alignment["query_name"]][alignment["ctg_name"]]:
            read_results[alignment["query_name"]][alignment["ctg_name"]] = identity

    ref_results = {ref: 0 for ref in possible_references.keys()}

    for read, details in read_results.items():
        best_ctg = max(details, key=details.get)
        ref_results[best_ctg] += 1

    best_ref = max(ref_results, key=ref_results.get)

    return possible_references[best_ref]
