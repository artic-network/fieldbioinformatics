import sys
from Bio import SeqIO
import tempfile
import os
import glob
import gzip
import fnmatch
import shutil
import pandas as pd
from collections import defaultdict
from mimetypes import guess_type
from functools import partial
from math import log10
from random import random
from multiprocessing import Pool

# the ambition of this module is to merge some of the functionality from the "gather" and "demultiplex" tasks;
# we are assuming the input is already guppy-demultiplexed data (one-pot barcoding - like) - this workflow will
# also assume that Medaka rather than Nanopolish will be used and will skip the sequencing_summary step - basic
# QC metrics will be derived during the fastq parsing step ...

# This method will also allow for the analysis of gzipped fastq files = makes sense of space

# Have tried to be minimally invasive to the existing code - maintain FieldBioinformatics style

# I have added multiprocessing functionality to accelerate analysis of large datasets. The subprocesses of the
# worker pool each process one fastq file at a time and write their results in temporary fastq files
# named .tmp_[worker_id] to reduce memory load. These are then merged with sensitivity to duplicates, 
# as in the original implementation. I set the default number of processes to 4.


def get_read_mean_quality(record):
    return -10 * log10((10 ** (pd.Series(record.letter_annotations["phred_quality"]) / -10)).mean())


def worker(infn, outfn, args):
    dups = set()

    _open = open
    outfh = open(outfn, "w")
    # only accommodating gzip compression at present
    encoding = guess_type(infn)[1]
    if encoding == "gzip":
        _open = partial(gzip.open, mode="rt")
    with _open(infn) as f:
        try:
            for rec in SeqIO.parse(f, "fastq"):
                if args.max_length and len(rec) > args.max_length:
                    continue
                if args.min_length and len(rec) < args.min_length:
                    continue
                if not args.skip_quality_check and get_read_mean_quality(rec) < args.quality:
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
    return outfn, dups


def run(parser, args):
    files = os.listdir(args.directory)
    fastq_files = [os.path.join(args.directory, f) for f in files if fnmatch.fnmatch(f, '*.fastq*') and not f.endswith('.temp')]

    if fastq_files:
        pool = Pool(args.processes)

        if not args.output:
            fastq_outfn = "%s_%s.fastq" % (args.prefix, os.path.basename(barcode_directory))
        else:
            fastq_outfn = args.output

        print("Processing %s files in %s" % (len(fastq_files), args.directory), file=sys.stderr)

        dups = set()

        results_array = []
        for i,fn in enumerate(fastq_files):
            batch_fn = os.path.join(os.path.dirname(fastq_outfn), ".tmp-{}_{}".format(i, os.path.basename(fastq_outfn)))
            results_array.append(pool.apply_async(worker, args=(fn, batch_fn, args)))

        pool.close()
        pool.join()

        with open(fastq_outfn, "w") as outfh:
            for result_set in results_array:
                batch_fn, _dups = result_set.get()
                duplicates = dups & _dups
                if duplicates:
                    for rec in SeqIO.parse(batch_fn, "fastq"):
                        if rec.id not in duplicates:
                            SeqIO.write([rec], outfh, "fastq")
                else:
                    with open(batch_fn, "r") as f:
                        for line in f.readlines():
                            outfh.write(line)
                os.remove(batch_fn)
                dups = dups | _dups

        print(f"{fastq_outfn}\t{len(dups)}")
