from artic.utils import get_scheme
from pathlib import Path


def main():
    import argparse
    import os

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--scheme-name",
        type=str,
        required=True,
        help="Name of the scheme you wish to fetch",
    )
    parser.add_argument(
        "--scheme-version",
        type=str,
        required=True,
        help="Version of the scheme you wish to fetch",
    )
    parser.add_argument(
        "--scheme-directory",
        type=Path,
        help="Directory where schemes are stored",
        default=f"{os.getcwd()}/primer-schemes",
    )
    parser.add_argument(
        "--scheme-length",
        type=str,
        help="Length of the scheme to fetch, only required if the scheme has multiple possible lengths.",
    )
    parser.add_argument(
        "--read-file",
        type=Path,
        help="FASTQ file containing reads sequenced with the scheme in question, this is only required when fetching a scheme with reference selection functionality.",
    )
    args = parser.parse_args()

    get_scheme(
        scheme_name=args.scheme_name,
        scheme_version=args.scheme_version,
        scheme_directory=args.scheme_directory,
        scheme_length=args.scheme_length,
        read_file=args.read_file,
    )


if __name__ == "__main__":
    main()
