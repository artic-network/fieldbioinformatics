import sys
from Bio import SeqIO
import os
import gzip
import fnmatch
import pandas as pd
from mimetypes import guess_type
from functools import partial
from math import log10
from random import random

# the ambition of this module is to merge some of the functionality from the "gather" and "demultiplex" tasks;
# we are assuming the input is already guppy-demultiplexed data (one-pot barcoding - like) - this workflow will
# also assume that Medaka rather than Nanopolish will be used and will skip the sequencing_summary step - basic
# QC metrics will be derived during the fastq parsing step ...

# This method will also allow for the analysis of gzipped fastq files = makes sense of space

# Have tried to be minimally invasive to the existing code - maintain FieldBioinformatics style


def get_read_mean_quality(record):
    return -10 * log10(
        (10 ** (pd.Series(record.letter_annotations["phred_quality"]) / -10)).mean()
    )


def run(parser, args):
    files = os.listdir(args.directory)
    fastq_files = [
        os.path.join(args.directory, f)
        for f in files
        if fnmatch.fnmatch(f, "*.fastq*") and not f.endswith(".temp")
    ]

    if fastq_files:
        if not args.output:
            fastq_outfn = "%s_%s.fastq" % (
                args.prefix,
                os.path.basename(args.directory),
            )
        else:
            fastq_outfn = args.output

        outfh = open(fastq_outfn, "w")
        print(
            "Processing %s files in %s" % (len(fastq_files), args.directory),
            file=sys.stderr,
        )

        dups = set()

        for fn in fastq_files:
            encoding = guess_type(fn)[1]
            _open = open
            # only accommodating gzip compression at present
            if encoding == "gzip":
                _open = partial(gzip.open, mode="rt")
            with _open(fn) as f:
                try:
                    for rec in SeqIO.parse(f, "fastq"):
                        if args.max_length and len(rec) > args.max_length:
                            continue
                        if args.min_length and len(rec) < args.min_length:
                            continue
                        if (
                            not args.skip_quality_check
                            and get_read_mean_quality(rec) < args.quality
                        ):
                            continue
                        if args.sample < 1:
                            r = random()
                            if r >= args.sample:
                                continue

                        if rec.id not in dups:
                            SeqIO.write([rec], outfh, "fastq")
                            dups.add(rec.id)
                except ValueError:
                    pass

        outfh.close()
        print(f"{fastq_outfn}\t{len(dups)}")
