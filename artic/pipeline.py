#!/usr/bin/env python

#Written by Nick Loman (@pathogenomenick)
#Thanks to Aaron Quinlan for the argparse implementation from poretools.

import sys
import hashlib
import re
import argparse
import sqlite3
from artic import version

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
    if args.command == 'rampart':
        from . import rampart as submodule
    if args.command == 'filter':
        from . import filter as submodule
    if args.command == 'run':
        from . import run as submodule

    # run the chosen submodule.
    submodule.run(parser, args)

class ArgumentParserWithDefaults(argparse.ArgumentParser):
    def __init__(self, *args, **kwargs):
        super(ArgumentParserWithDefaults, self).__init__(*args, **kwargs)
        self.add_argument("-q", "--quiet", help="Do not output warnings to stderr",
                          action="store_true",
                          dest="quiet")

def main():
    parser = argparse.ArgumentParser(prog='artic', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-v", "--version", help="Installed Artic version",
                        action="version",
                        version="%(prog)s " + str(version.__version__))
    subparsers = parser.add_subparsers(title='[sub-commands]', dest='command', parser_class=ArgumentParserWithDefaults)

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
    parser_extract = subparsers.add_parser('basecaller', help='Display basecallers in files')
    parser_extract.add_argument('directory', metavar='directory', help='Directory of FAST5 files.')
    parser_extract.set_defaults(func=run_subtool)

    # demultiplex
    parser_demultiplex = subparsers.add_parser('demultiplex', help='Run demultiplex')
    parser_demultiplex.add_argument('fasta', metavar='fasta', help='Undemultiplexed FASTA file.')
    parser_demultiplex.add_argument('--threads', type=int, default=8, help='Number of threads')
    parser_demultiplex.add_argument('--prefix', help='Prefix for demultiplexed files')
    parser_demultiplex.add_argument('--no-remove-directory', dest='no_remove_directory', action='store_true')
    parser_demultiplex.set_defaults(func=run_subtool)

    # minion
    parser_minion = subparsers.add_parser('minion', help='Run demultiplex')
    parser_minion.add_argument('scheme', metavar='scheme', help='The name of the scheme.')
    parser_minion.add_argument('sample', metavar='sample', help='The name of the sample.')
    parser_minion.add_argument('--medaka', dest='medaka', type=bool, default=False, help='Use medaka instead of nanopolish for variants')
    parser_minion.add_argument('--medaka-model', dest='medaka_model', default='r941_min_high_g344_model.hdf5', help='Medaka model to use')
    parser_minion.add_argument('--normalise', dest='normalise', type=int, default=100, help='Normalise down to moderate coverage to save runtime.')
    parser_minion.add_argument('--threads', type=int, default=8, help='Number of threads')
    parser_minion.add_argument('--scheme-directory', metavar='scheme_directory', default='/artic/schemes', help='Default scheme directory')
    parser_minion.add_argument('--max-haplotypes', type=int, default=1000000, metavar='max_haplotypes', help='max-haplotypes value for nanopolish')
    parser_minion.add_argument('--read-file', metavar='read_file', help='Use alternative FASTA/FASTQ file to <sample>.fasta')
    parser_minion.add_argument('--nanopolish-read-file', metavar='nanopolish_read_file', help='Use alternative read file (previously indexed)')
    parser_minion.add_argument('--skip-nanopolish', action='store_true')
    parser_minion.set_defaults(func=run_subtool)

    # gather
    parser_gather = subparsers.add_parser('gather', help='Gather up demultiplexed files')
    parser_gather.add_argument('--directory', nargs='+', metavar='directory', help='Basecalled (guppy) results directory or directories.')
    parser_gather.add_argument('--max-length', type=int, metavar='max_length', help='remove reads greater than read length')
    parser_gather.add_argument('--min-length', type=int, metavar='min_length', help='remove reads less than read length')
    parser_gather.add_argument('--prefix', help='Prefix for gathered files')
    parser_gather.add_argument('--run-directory', metavar='run_directory', help='The run directory', default='/var/lib/MinKNOW/data')
    parser_gather.set_defaults(func=run_subtool)

    # filter
    parser_filter = subparsers.add_parser('filter', help='Filter FASTQ files by length')
    parser_filter.add_argument('filename', metavar='filename', help='FASTQ file.')
    parser_filter.add_argument('--max-length', type=int, metavar='max_length', help='remove reads greater than read length')
    parser_filter.add_argument('--min-length', type=int, metavar='min_length', help='remove reads less than read length')
    parser_filter.set_defaults(func=run_subtool)

    # rampart
    parser_rampart = subparsers.add_parser('rampart', help='Interactive prompts to start RAMPART')
    parser_rampart.add_argument('--protocol-directory', metavar='protocol_directory', help='The RAMPART protocols directory.', default='/home/artic/artic/artic-ebov/rampart')
    parser_rampart.add_argument('--run-directory', metavar='run_directory', help='The run directory', default='/var/lib/MinKNOW/data')
    parser_rampart.set_defaults(func=run_subtool)

    # run
    parser_run = subparsers.add_parser('run', help='Process an entire run folder interactively')
    parser_run.set_defaults(func=run_subtool)

    args = parser.parse_args()

    #if args.quiet:
    #    logger.setLevel(logging.ERROR)

    if args.command:
        args.func(parser, args)
    else:
        parser.print_usage()

if __name__ == "__main__":
    main()
