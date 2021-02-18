#!/usr/bin/env python

# Written by Nick Loman (@pathogenomenick)
# Thanks to Aaron Quinlan for the argparse implementation from poretools.

import argparse
import sys

from . import version


def run_subtool(parser, args):
    if args.command == 'extract':
        from . import extract as submodule
    if args.command == 'basecaller':
        from . import basecaller as submodule
    if args.command == 'demultiplex':
        from . import demultiplex as submodule
    if args.command == 'minion':
        from . import minion as submodule
    if args.command == 'gather':
        from . import gather as submodule
    if args.command == 'guppyplex':
        from . import guppyplex as submodule
    if args.command == 'rampart':
        from . import rampart as submodule
    if args.command == 'filter':
        from . import filter_reads as submodule
    if args.command == 'run':
        from . import run as submodule
    if args.command == 'export':
        from . import export as submodule

    # run the chosen submodule.
    submodule.run(parser, args)


class ArgumentParserWithDefaults(argparse.ArgumentParser):
    def __init__(self, *args, **kwargs):
        super(ArgumentParserWithDefaults, self).__init__(*args, **kwargs)
        self.add_argument("-q", "--quiet", help="Do not output warnings to stderr",
                          action="store_true",
                          dest="quiet")


def init_pipeline_parser():
    """Wraps the argparse parser initialisation.

    Returns
    -------
    argparse.ArgumentParser
        The initialised argparse Argument Parser for the pipeline
    """
    parser = argparse.ArgumentParser(
        prog='artic', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-v", "--version", help="Installed Artic version",
                        action="version",
                        version="%(prog)s " + str(version.__version__))
    subparsers = parser.add_subparsers(
        title='[sub-commands]', dest='command', parser_class=ArgumentParserWithDefaults)

    # extract
    parser_extract = subparsers.add_parser('extract',
                                           help='Create an empty poredb database')
    parser_extract.add_argument('directory', metavar='directory',
                                help='The name of the database.')
    parser_extract.add_argument('--basecaller', metavar='basecaller',
                                default='ONT Albacore Sequencing Software',
                                help='The name of the basecaller')
    parser_extract.set_defaults(func=run_subtool)

    # callers
    parser_extract = subparsers.add_parser(
        'basecaller', help='Display basecallers in files')
    parser_extract.add_argument(
        'directory', metavar='directory', help='Directory of FAST5 files.')
    parser_extract.set_defaults(func=run_subtool)

    # demultiplex
    parser_demultiplex = subparsers.add_parser(
        'demultiplex', help='Run demultiplex')
    parser_demultiplex.add_argument(
        'fasta', metavar='fasta', help='Undemultiplexed FASTA file.')
    parser_demultiplex.add_argument(
        '--threads', type=int, default=8, help='Number of threads')
    parser_demultiplex.add_argument(
        '--prefix', help='Prefix for demultiplexed files')
    parser_demultiplex.add_argument(
        '--no-remove-directory', dest='no_remove_directory', action='store_true')
    parser_demultiplex.set_defaults(func=run_subtool)

    # minion
    parser_minion = subparsers.add_parser('minion', help='Run the alignment/variant-call/consensus pipeline')
    parser_minion.add_argument(
        'scheme', metavar='scheme', help='The name of the scheme')
    parser_minion.add_argument(
        'sample', metavar='sample', help='The name of the sample')
    parser_minion.add_argument('--medaka', dest='medaka', action='store_true',
                               help='Use medaka instead of nanopolish for variants')
    parser_minion.add_argument('--medaka-model', metavar='medaka_model', help='The model to use for medaka (required if using --medaka)')
    parser_minion.add_argument('--minimap2', dest='minimap2', default=True,
                               action='store_true', help='Use minimap2 (default)')
    parser_minion.add_argument(
        '--bwa', dest='bwa', action='store_true', help='Use bwa instead of minimap2')
    parser_minion.add_argument('--normalise', dest='normalise', type=int,
                               default=100, help='Normalise down to moderate coverage to save runtime (default: %(default)d, deactivate with `--normalise 0`)')
    parser_minion.add_argument(
        '--threads', type=int, default=8, help='Number of threads (default: %(default)d)')
    parser_minion.add_argument('--scheme-directory', metavar='scheme_directory',
                               default='./primer-schemes', help='Default scheme directory')
    parser_minion.add_argument('--scheme-version', metavar='scheme_version',
                                default=1, help='Primer scheme version (default: %(default)d)')
    parser_minion.add_argument('--max-haplotypes', type=int, default=1000000,
                               metavar='max_haplotypes', help='max-haplotypes value for nanopolish')
    parser_minion.add_argument('--min-depth', type=int, default=20,
                                metavar='min_depth', help='minimum read depth per-base for consensus masking')
    parser_minion.add_argument('--min-qual', type=int, default=20,
                                metavar='min_qual', help='minimum variant quality for consensus masking')
    parser_minion.add_argument('--nanopolish-qual-cov-ratio', type=float, default=3.0, metavar='nanopolish_qual_cov_ratio', help='the quality/coverage cutoff for filtering nanopolish variants')
    parser_minion.add_argument('--read-file', metavar='read_file',
                               help='Use alternative FASTA/FASTQ file to <sample>.fasta')
    parser_minion.add_argument('--fast5-directory', help='FAST5 Directory')
    parser_minion.add_argument(
        '--sequencing-summary', help='Path to Guppy sequencing summary')
    parser_minion.add_argument('--skip-nanopolish', action='store_true')
    parser_minion.add_argument('--no-indels', action='store_true', help='Do not report InDels (uses SNP-only mode of nanopolish/medaka)')
    parser_minion.add_argument('--no-frameshifts', action='store_true', help='Remove variants which induce frameshifts (ignored when --no-indels set)')
    parser_minion.add_argument('--dry-run', action='store_true')
    parser_minion.add_argument('--strict', action='store_true', help='Run with full-length amplicon alignments & strict filtering of variants against primer scheme')
    parser_minion.set_defaults(func=run_subtool)

    # gather
    parser_gather = subparsers.add_parser(
        'gather', help='Gather up demultiplexed files')
    parser_gather.add_argument('--directory', nargs='+', metavar='directory',
                               help='Basecalled (guppy) results directory or directories.')
    parser_gather.add_argument(
        '--max-length', type=int, metavar='max_length', help='remove reads greater than read length')
    parser_gather.add_argument(
        '--min-length', type=int, metavar='min_length', help='remove reads less than read length')
    parser_gather.add_argument('--prefix', help='Prefix for gathered files')
    parser_gather.add_argument('--prompt-directory', metavar='run_directory',
                               help='The run directory for interactive prompts', default='/var/lib/minknown/data')
    parser_gather.add_argument(
        '--fast5-directory', metavar='fast5_directory', help='The directory with fast5 files')
    parser_gather.add_argument('--no-fast5s', action='store_true',
                               help='Do not use fast5s and nanopolish', default=0)
    parser_gather.add_argument('--limit', type=int, help='Only gather n reads')
    parser_gather.set_defaults(func=run_subtool)

    # guppyplex
    # This is a workflow that aggregates the previous gather and demultiplex steps into a single task.
    # This is making an assumption that the results from MinKnow demultiplex are good-enough.
    parser_guppyplex = subparsers.add_parser(
        'guppyplex', help='Aggregate pre-demultiplexed reads from MinKNOW/Guppy')
    parser_guppyplex.add_argument('--directory', metavar='directory',
                                  help='Basecalled and demultiplexed (guppy) results directory', required=True)
    parser_guppyplex.add_argument(
        '--max-length', type=int, metavar='max_length', help='remove reads greater than read length')
    parser_guppyplex.add_argument(
        '--min-length', type=int, metavar='min_length', help='remove reads less than read length')
    parser_guppyplex.add_argument('--quality', type=float, metavar='quality',
                                  default=7, help='remove reads against this quality filter')
    parser_guppyplex.add_argument('--sample', type=float, metavar='sample', default=1,
                                  help='sampling frequency for random sample of sequence to reduce excess')
    parser_guppyplex.add_argument(
        '--skip-quality-check', action='store_true', help='Do not filter on quality score (speeds up)')
    parser_guppyplex.add_argument(
        '--prefix', help='Prefix for guppyplex files')
    parser_guppyplex.add_argument('--output', metavar='output',
                                  help='FASTQ file to write')
    parser_guppyplex.set_defaults(func=run_subtool)

    # filter
    parser_filter = subparsers.add_parser(
        'filter', help='Filter FASTQ files by length')
    parser_filter.add_argument(
        'filename', metavar='filename', help='FASTQ file.')
    parser_filter.add_argument(
        '--max-length', type=int, metavar='max_length', help='remove reads greater than read length')
    parser_filter.add_argument(
        '--min-length', type=int, metavar='min_length', help='remove reads less than read length')
    parser_filter.set_defaults(func=run_subtool)

    # rampart
    parser_rampart = subparsers.add_parser(
        'rampart', help='Interactive prompts to start RAMPART')
    parser_rampart.add_argument('--protocol-directory', metavar='protocol_directory',
                                help='The RAMPART protocols directory.', default='/home/artic/artic/artic-ebov/rampart')
    parser_rampart.add_argument('--run-directory', metavar='run_directory',
                                help='The run directory', default='/var/lib/MinKNOW/data')
    parser_rampart.set_defaults(func=run_subtool)

    # export
    parser_export = subparsers.add_parser(
        'export', help='Export reads and fAST5 into a neat archive')
    parser_export.add_argument('prefix')
    parser_export.add_argument('bamfile')
    parser_export.add_argument('sequencing_summary')
    parser_export.add_argument('fast5_directory')
    parser_export.add_argument('output_directory')
    parser_export.set_defaults(func=run_subtool)

    # run
    parser_run = subparsers.add_parser(
        'run', help='Process an entire run folder interactively')
    parser_run.set_defaults(func=run_subtool)

    # return the parser
    return parser


def main():

    # init the pipeline parser
    parser = init_pipeline_parser()

    # collect the args
    args = parser.parse_args(sys.argv[1:])

    # if args.quiet:
    #    logger.setLevel(logging.ERROR)

    # run the subcommand or print usage if no subcommand provided
    if args.command:
        args.func(parser, args)
    else:
        parser.print_usage()


if __name__ == "__main__":
    main()
