import sys
import os
import gzip
import fnmatch
import concurrent.futures
import numpy as np
from Bio.SeqIO.QualityIO import FastqGeneralIterator
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
    quals = np.asarray(record.letter_annotations["phred_quality"], dtype=np.float64)
    return -10 * log10(np.mean(10.0 ** (quals / -10.0)))


def _mean_phred_quality_from_ascii(qual):
    """Compute mean Phred quality directly from a raw FASTQ quality string.

    Avoids constructing a SeqRecord (and its list of int quality scores) just
    to compute a mean - operates on the ASCII string Biopython already parsed.
    """
    scores = np.frombuffer(qual.encode("ascii"), dtype=np.uint8).astype(np.float64) - 33
    return -10 * log10(np.mean(10.0 ** (scores / -10.0)))


def _process_file(args_tuple):
    """Module-level worker: filter reads from a single FASTQ file.

    Returns a list of (read_id, formatted_fastq_record) tuples for reads that
    pass all filters. Records are pre-formatted here (in the worker process)
    so the main process only has to deduplicate and write raw strings -
    Bio.SeqIO's SeqRecord construction/serialisation is the dominant cost of
    this pipeline stage and is avoided entirely by working with the raw
    title/sequence/quality strings that Bio.SeqIO parses FASTQ into internally.

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
            for title, seq, qual in FastqGeneralIterator(f):
                seq_len = len(seq)
                if max_length and seq_len > max_length:
                    continue
                if min_length and seq_len < min_length:
                    continue
                if not skip_quality_check and _mean_phred_quality_from_ascii(qual) < quality:
                    continue
                if sample < 1:
                    r = random()
                    if r >= sample:
                        continue
                read_id = title.split(maxsplit=1)[0]
                records.append((read_id, f"@{title}\n{seq}\n+\n{qual}\n"))
        except (ValueError, gzip.BadGzipFile, EOFError) as e:
            print(f"Warning: skipping {fn}: {e}", file=sys.stderr)
    return records


def run(parser, args):
    files = os.listdir(args.directory)
    fastq_files = [
        os.path.join(args.directory, f)
        for f in files
        if fnmatch.fnmatch(f, "*.fastq*")
        and not f.endswith(".temp")
        and not f.startswith(".")
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
                for read_id, text in file_records:
                    if read_id not in dups:
                        outfh.write(text)
                        dups.add(read_id)

        outfh.close()
        print(f"{fastq_outfn}\t{len(dups)}")
