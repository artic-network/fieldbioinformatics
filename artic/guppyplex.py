import sys
import os
import gzip
import fnmatch
import tempfile
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


def _split_file(args_tuple):
    """Module-level worker: split one plain-text FASTQ file into n shards.

    Modern basecallers (e.g. Dorado) commonly write a single large FASTQ per
    barcode rather than many smaller chunks, which leaves nothing for
    --threads to parallelise over in the filtering pass below (that pass is
    file-granular). This does one cheap sequential pass - round-robin
    re-emitting each record with no quality computation - so the expensive
    filtering pass has n independent, roughly-equal work units to run in
    parallel instead of one. gzip input isn't handled here since it can't be
    split without a full decompress first; callers should only offer plain-
    text files to this function.

    Returns the list of shard file paths (some may end up empty for small
    inputs, which is harmless).
    """
    fn, n = args_tuple
    shard_files = [
        tempfile.NamedTemporaryFile(
            mode="w", suffix=".fastq", prefix="guppyplex_split_", delete=False
        )
        for _ in range(n)
    ]
    try:
        with open(fn) as f:
            try:
                for i, (title, seq, qual) in enumerate(FastqGeneralIterator(f)):
                    shard_files[i % n].write(f"@{title}\n{seq}\n+\n{qual}\n")
            except (ValueError, EOFError) as e:
                print(f"Warning: skipping {fn}: {e}", file=sys.stderr)
    finally:
        for sf in shard_files:
            sf.close()
    return [sf.name for sf in shard_files]


def _process_file(args_tuple):
    """Module-level worker: filter reads from a single FASTQ file.

    Passing reads are streamed straight to a per-worker temp file on disk as
    they're found, rather than accumulated in a list and returned - for a
    large input file that keeps this worker's memory use O(1) instead of
    O(filtered file size), and avoids pickling a potentially huge list of
    strings back through the process pool's IPC pipe.

    Returns the temp file path (or None if nothing passed the filters);
    deduplication across files is handled by the caller.
    """
    fn, min_length, max_length, quality, skip_quality_check, sample = args_tuple
    encoding = guess_type(fn)[1]
    _open = open
    if encoding == "gzip":
        _open = partial(gzip.open, mode="rt")

    tmp = tempfile.NamedTemporaryFile(
        mode="w", suffix=".fastq", prefix="guppyplex_", delete=False
    )
    wrote_any = False
    try:
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
                    tmp.write(f"@{title}\n{seq}\n+\n{qual}\n")
                    wrote_any = True
            except (ValueError, gzip.BadGzipFile, EOFError) as e:
                print(f"Warning: skipping {fn}: {e}", file=sys.stderr)
    finally:
        tmp.close()

    if not wrote_any:
        os.unlink(tmp.name)
        return None
    return tmp.name


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

        threads = getattr(args, "threads", 1)

        # If there aren't enough (splittable, plain-text) input files to keep
        # every worker busy - e.g. a single big per-barcode FASTQ, as newer
        # basecallers tend to produce - pre-split the large ones into
        # `threads` shards each so the filtering pass below actually has
        # enough independent work units to parallelise over.
        split_shards = []
        work_files = fastq_files
        if threads > 1 and len(fastq_files) < threads:
            splittable = [fn for fn in fastq_files if guess_type(fn)[1] != "gzip"]
            unsplittable = [fn for fn in fastq_files if guess_type(fn)[1] == "gzip"]
            if splittable:
                with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
                    for shards in executor.map(
                        _split_file, [(fn, threads) for fn in splittable]
                    ):
                        split_shards.extend(shards)
                work_files = split_shards + unsplittable

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
            for fn in work_files
        ]

        with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
            futures = [
                executor.submit(_process_file, wa) for wa in worker_args
            ]
            # as_completed (rather than executor.map) merges each worker's
            # shard as soon as it's ready, instead of buffering results that
            # finish early while waiting on an earlier-submitted file.
            for future in concurrent.futures.as_completed(futures):
                shard_path = future.result()
                if shard_path is None:
                    continue
                try:
                    with open(shard_path) as shard:
                        while True:
                            title_line = shard.readline()
                            if not title_line:
                                break
                            seq_line = shard.readline()
                            plus_line = shard.readline()
                            qual_line = shard.readline()
                            read_id = title_line[1:].split(maxsplit=1)[0]
                            if read_id not in dups:
                                outfh.write(title_line)
                                outfh.write(seq_line)
                                outfh.write(plus_line)
                                outfh.write(qual_line)
                                dups.add(read_id)
                finally:
                    os.unlink(shard_path)

        for sp in split_shards:
            if os.path.exists(sp):
                os.unlink(sp)

        outfh.close()
        print(f"{fastq_outfn}\t{len(dups)}")
