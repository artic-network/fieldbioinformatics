#!/usr/bin/env python

# Written by Nick Loman (@pathogenomenick)
# Thanks to Aaron Quinlan for the argparse implementation from poretools.

import argparse
import sys
import os

from importlib.metadata import version


def run_subtool(parser, args):
    if args.command == "minion":
        from . import minion as submodule
    if args.command == "guppyplex":
        from . import guppyplex as submodule
    if args.command == "rampart":
        from . import rampart as submodule
    if args.command == "filter":
        from . import filter_reads as submodule
    if args.command == "run":
        from . import run as submodule
    if args.command == "export":
        from . import export as submodule

    # run the chosen submodule.
    submodule.run(parser, args)


class ArgumentParserWithDefaults(argparse.ArgumentParser):
    def __init__(self, *args, **kwargs):
        super(ArgumentParserWithDefaults, self).__init__(*args, **kwargs)
        self.add_argument(
            "-q",
            "--quiet",
            help="Do not output warnings to stderr",
            action="store_true",
            dest="quiet",
        )


def init_pipeline_parser():
    """Wraps the argparse parser initialisation.

    Returns
    -------
    argparse.ArgumentParser
        The initialised argparse Argument Parser for the pipeline
    """
    parser = argparse.ArgumentParser(
        prog="artic", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-v",
        "--version",
        help="Installed Artic version",
        action="version",
        version=f"artic {version('artic')}",
    )
    subparsers = parser.add_subparsers(
        title="[sub-commands]", dest="command", parser_class=ArgumentParserWithDefaults
    )

    # minion
    parser_minion = subparsers.add_parser(
        "minion", help="Run the alignment/variant-call/consensus pipeline"
    )
    parser_minion.add_argument(
        "sample", metavar="sample", help="The name of the sample"
    )
    parser_minion.add_argument(
        "--model",
        help="The model to use for clair3, if not provided the pipeline will try to figure it out the appropriate model from the read fastq",
    )
    parser_minion.add_argument(
        "--model-dir",
        metavar="model_path",
        help="Path containing clair3 models, defaults to models packaged with conda installation (default: $CONDA_PREFIX/bin/models/)",
        type=str,
    )
    parser_minion.add_argument(
        "--min-mapq",
        type=int,
        default=20,
        help="Minimum mapping quality to consider (default: %(default)d)",
    )
    parser_minion.add_argument(
        "--normalise",
        dest="normalise",
        type=int,
        default=100,
        help="Normalise down to moderate coverage to save runtime (default: %(default)d, deactivate with `--normalise 0`)",
    )
    parser_minion.add_argument(
        "--primer-match-threshold",
        type=int,
        default=35,
        help="Allow fuzzy primer matching within this threshold (default: %(default)d)",
    )
    parser_minion.add_argument(
        "--min-depth",
        type=int,
        default=20,
        help="Minimum depth required to call a variant (default: %(default)d)",
    )
    parser_minion.add_argument(
        "--threads",
        type=int,
        default=8,
        help="Number of threads (default: %(default)d)",
    )

    remote_scheme_options = parser_minion.add_argument_group(
        "Remote Scheme Options",
        "Options for automagically fetching primer schemes from the scheme repository (https://github.com/quick-lab/primerschemes)",
    )
    remote_scheme_options.add_argument(
        "--scheme-directory",
        metavar="scheme_directory",
        default=f"{os.getcwd()}/primer-schemes",
        help="Default scheme directory (default: %(default)s)",
    )
    remote_scheme_options.add_argument(
        "--scheme-name",
        metavar="scheme_name",
        type=str,
        help="Name of the scheme, e.g. sars-cov-2, mpxv to fetch from the scheme repository (https://github.com/quick-lab/primerschemes)",
    )
    remote_scheme_options.add_argument(
        "--scheme-length",
        type=int,
        metavar="scheme_length",
        default=False,
        help="Length of the scheme to fetch from the scheme repository (https://github.com/quick-lab/primerschemes). This is only required if the --scheme-name has multiple possible lengths.",
    )
    remote_scheme_options.add_argument(
        "--scheme-version",
        metavar="vX.X.X",
        type=str,
        help="Primer scheme version",
    )

    local_scheme_options = parser_minion.add_argument_group(
        "Local Scheme Options",
        "Options for providing a local primer scheme and reference sequence (BED file and FASTA file) directly to the pipeline",
    )
    local_scheme_options.add_argument(
        "--bed",
        metavar="bed",
        help="BED file containing primer scheme",
    )
    local_scheme_options.add_argument(
        "--ref",
        help="Reference sequence for the scheme",
    )

    parser_minion.add_argument(
        "--read-file",
        metavar="read_file",
        help="FASTQ file containing reads",
        required=True,
    )
    parser_minion.add_argument(
        "--no-indels",
        action="store_true",
        help="Do not report InDels (uses SNP-only mode of nanopolish/medaka)",
    )
    parser_minion.add_argument(
        "--no-frameshifts",
        action="store_true",
        help="Remove variants which induce frameshifts (ignored when --no-indels set)",
    )
    parser_minion.add_argument(
        "--align-consensus",
        action="store_true",
        help="Run a mafft alignment of consensus to reference after generation",
    )
    parser_minion.add_argument("--dry-run", action="store_true")
    parser_minion.set_defaults(func=run_subtool)

    # guppyplex
    # This is a workflow that aggregates the previous gather and demultiplex steps into a single task.
    # This is making an assumption that the results from MinKnow demultiplex are good-enough.
    parser_guppyplex = subparsers.add_parser(
        "guppyplex", help="Aggregate pre-demultiplexed reads from MinKNOW/Guppy"
    )
    parser_guppyplex.add_argument(
        "--directory",
        metavar="directory",
        help="Basecalled and demultiplexed (guppy) results directory",
        required=True,
    )
    parser_guppyplex.add_argument(
        "--max-length",
        type=int,
        metavar="max_length",
        help="remove reads greater than read length",
    )
    parser_guppyplex.add_argument(
        "--min-length",
        type=int,
        metavar="min_length",
        help="remove reads less than read length",
    )
    parser_guppyplex.add_argument(
        "--quality",
        type=float,
        metavar="quality",
        default=7,
        help="remove reads against this quality filter",
    )
    parser_guppyplex.add_argument(
        "--sample",
        type=float,
        metavar="sample",
        default=1,
        help="sampling frequency for random sample of sequence to reduce excess",
    )
    parser_guppyplex.add_argument(
        "--skip-quality-check",
        action="store_true",
        help="Do not filter on quality score (speeds up)",
    )
    parser_guppyplex.add_argument("--prefix", help="Prefix for guppyplex files")
    parser_guppyplex.add_argument(
        "--output", metavar="output", help="FASTQ file to write"
    )
    parser_guppyplex.set_defaults(func=run_subtool)

    # filter
    parser_filter = subparsers.add_parser("filter", help="Filter FASTQ files by length")
    parser_filter.add_argument("filename", metavar="filename", help="FASTQ file.")
    parser_filter.add_argument(
        "--max-length",
        type=int,
        metavar="max_length",
        help="remove reads greater than read length",
    )
    parser_filter.add_argument(
        "--min-length",
        type=int,
        metavar="min_length",
        help="remove reads less than read length",
    )
    parser_filter.set_defaults(func=run_subtool)

    # rampart
    parser_rampart = subparsers.add_parser(
        "rampart", help="Interactive prompts to start RAMPART"
    )
    parser_rampart.add_argument(
        "--protocol-directory",
        metavar="protocol_directory",
        help="The RAMPART protocols directory.",
        default="/home/artic/artic/artic-ebov/rampart",
    )
    parser_rampart.add_argument(
        "--run-directory",
        metavar="run_directory",
        help="The run directory",
        default="/var/lib/MinKNOW/data",
    )
    parser_rampart.set_defaults(func=run_subtool)

    # export
    parser_export = subparsers.add_parser(
        "export", help="Export reads and fAST5 into a neat archive"
    )
    parser_export.add_argument("prefix")
    parser_export.add_argument("bamfile")
    parser_export.add_argument("sequencing_summary")
    parser_export.add_argument("fast5_directory")
    parser_export.add_argument("output_directory")
    parser_export.set_defaults(func=run_subtool)

    # run
    parser_run = subparsers.add_parser(
        "run", help="Process an entire run folder interactively"
    )
    parser_run.set_defaults(func=run_subtool)

    # return the parser
    return parser


def main():

    # init the pipeline parser
    parser = init_pipeline_parser()

    # collect the args
    args = parser.parse_args()

    # if args.quiet:
    #    logger.setLevel(logging.ERROR)

    # run the subcommand or print usage if no subcommand provided
    if args.command:
        args.func(parser, args)
    else:
        parser.print_usage()


if __name__ == "__main__":
    main()
