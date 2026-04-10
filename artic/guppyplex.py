import sys
from Bio import SeqIO
import os
import gzip
import fnmatch
import concurrent.futures
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


def _process_file(args_tuple):
    """Module-level worker: filter reads from a single FASTQ file.

    Returns a list of SeqRecord objects that pass all filters.
    Deduplication is handled by the caller across all files.
    """
    fn, min_length, max_length, quality, skip_quality_check, sample = args_tuple
    encoding = guess_type(fn)[1]
    _open = open
    if encoding == "gzip":
        _open = partial(gzip.open, mode="rt")

    records = []
    with _open(fn) as f:
        try:
            for rec in SeqIO.parse(f, "fastq"):
                if max_length and len(rec) > max_length:
                    continue
                if min_length and len(rec) < min_length:
                    continue
                if not skip_quality_check and get_read_mean_quality(rec) < quality:
                    continue
                if sample < 1:
                    r = random()
                    if r >= sample:
                        continue
                records.append(rec)
        except ValueError:
            pass
    return records


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

        if fastq_outfn.lower().endswith(".gz"):
            outfh = gzip.open(fastq_outfn, "wt")
        else:
            outfh = open(fastq_outfn, "w")

        print(
            "Processing %s files in %s" % (len(fastq_files), args.directory),
            file=sys.stderr,
        )

        dups = set()
        worker_args = [
            (
                fn,
                args.min_length,
                args.max_length,
                args.quality,
                args.skip_quality_check,
                args.sample,
            )
            for fn in fastq_files
        ]

        threads = getattr(args, "threads", 1)
        with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
            for file_records in executor.map(_process_file, worker_args):
                for rec in file_records:
                    if rec.id not in dups:
                        SeqIO.write([rec], outfh, "fastq")
                        dups.add(rec.id)

        outfh.close()
        print(f"{fastq_outfn}\t{len(dups)}")
