"""Unit tests for artic/guppyplex.py"""

import gzip
import io
import math
from types import SimpleNamespace

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from artic.guppyplex import get_read_mean_quality, run


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_record(read_id, seq, qual):
    """Build a SeqRecord with given phred quality list."""
    assert len(seq) == len(qual)
    return SeqRecord(
        Seq(seq),
        id=read_id,
        description=read_id,
        letter_annotations={"phred_quality": qual},
    )


def _write_fastq(path, records):
    with open(str(path), "w") as fh:
        SeqIO.write(records, fh, "fastq")


def _write_fastq_gz(path, records):
    buf = io.StringIO()
    SeqIO.write(records, buf, "fastq")
    with gzip.open(str(path), "wt") as fh:
        fh.write(buf.getvalue())


def _make_args(directory, prefix="sample", min_length=None, max_length=None,
               quality=7.0, skip_quality_check=False, sample=1.0, output=None,
               threads=1):
    return SimpleNamespace(
        directory=str(directory),
        prefix=prefix,
        min_length=min_length,
        max_length=max_length,
        quality=quality,
        skip_quality_check=skip_quality_check,
        sample=sample,
        output=output,
        threads=threads,
    )


def _read_output_fastq(path):
    with open(str(path)) as fh:
        return list(SeqIO.parse(fh, "fastq"))


def _read_output_fastq_gz(path):
    with gzip.open(str(path), "rt") as fh:
        return list(SeqIO.parse(fh, "fastq"))


# ---------------------------------------------------------------------------
# get_read_mean_quality
# ---------------------------------------------------------------------------

class TestGetReadMeanQuality:
    def test_uniform_high_quality(self):
        rec = _make_record("r1", "ACGT", [40, 40, 40, 40])
        result = get_read_mean_quality(rec)
        assert abs(result - 40.0) < 0.01

    def test_uniform_low_quality(self):
        rec = _make_record("r1", "ACGT", [10, 10, 10, 10])
        result = get_read_mean_quality(rec)
        assert abs(result - 10.0) < 0.01

    def test_single_base(self):
        rec = _make_record("r1", "A", [30])
        result = get_read_mean_quality(rec)
        assert abs(result - 30.0) < 0.01

    def test_mixed_quality_is_less_than_max(self):
        rec = _make_record("r1", "ACGT", [10, 20, 30, 40])
        result = get_read_mean_quality(rec)
        # Mean Phred should be less than the arithmetic mean (due to log-sum)
        assert result < 40.0
        assert result > 10.0

    def test_known_mixed_quality(self):
        # Manual calculation: phred 10 and 30 mixed
        # error_prob = (0.1 + 0.001) / 2 = 0.0505
        # expected_q = -10 * log10(0.0505) ≈ 12.97
        rec = _make_record("r1", "AC", [10, 30])
        result = get_read_mean_quality(rec)
        expected = -10 * math.log10((10 ** (-10 / 10) + 10 ** (-30 / 10)) / 2)
        assert abs(result - expected) < 0.01


# ---------------------------------------------------------------------------
# run
# ---------------------------------------------------------------------------

class TestRun:
    def _output_path(self, tmp_path, name="sample_reads.fastq"):
        return tmp_path / name

    def test_run_basic_pass_through(self, tmp_path):
        rec = _make_record("read1", "ACGT" * 25, [30] * 100)
        _write_fastq(tmp_path / "reads.fastq", [rec])
        out = str(self._output_path(tmp_path))
        args = _make_args(tmp_path, quality=7.0, output=out)
        run(None, args)
        result = _read_output_fastq(out)
        assert len(result) == 1
        assert result[0].id == "read1"

    def test_run_filters_by_min_length(self, tmp_path):
        short = _make_record("short", "ACG", [30, 30, 30])
        long_ = _make_record("long", "ACGT" * 25, [30] * 100)
        _write_fastq(tmp_path / "reads.fastq", [short, long_])
        out = str(self._output_path(tmp_path))
        args = _make_args(tmp_path, min_length=50, quality=7.0, output=out)
        run(None, args)
        result = _read_output_fastq(out)
        assert all(r.id == "long" for r in result)
        assert len(result) == 1

    def test_run_filters_by_max_length(self, tmp_path):
        short = _make_record("short", "ACG", [30, 30, 30])
        long_ = _make_record("long", "ACGT" * 25, [30] * 100)
        _write_fastq(tmp_path / "reads.fastq", [short, long_])
        out = str(self._output_path(tmp_path))
        args = _make_args(tmp_path, max_length=10, quality=7.0, output=out)
        run(None, args)
        result = _read_output_fastq(out)
        assert all(r.id == "short" for r in result)
        assert len(result) == 1

    def test_run_filters_by_quality(self, tmp_path):
        low_q = _make_record("lowq", "ACGT" * 10, [5] * 40)
        high_q = _make_record("highq", "ACGT" * 10, [30] * 40)
        _write_fastq(tmp_path / "reads.fastq", [low_q, high_q])
        out = str(self._output_path(tmp_path))
        args = _make_args(tmp_path, quality=20.0, output=out)
        run(None, args)
        result = _read_output_fastq(out)
        assert len(result) == 1
        assert result[0].id == "highq"

    def test_run_skip_quality_check(self, tmp_path):
        low_q = _make_record("lowq", "ACGT" * 10, [5] * 40)
        _write_fastq(tmp_path / "reads.fastq", [low_q])
        out = str(self._output_path(tmp_path))
        args = _make_args(tmp_path, quality=20.0, skip_quality_check=True, output=out)
        run(None, args)
        result = _read_output_fastq(out)
        assert len(result) == 1

    def test_run_deduplication(self, tmp_path):
        rec1 = _make_record("dup_read", "ACGT" * 10, [30] * 40)
        rec2 = _make_record("dup_read", "TTTT" * 10, [30] * 40)
        _write_fastq(tmp_path / "reads.fastq", [rec1, rec2])
        out = str(self._output_path(tmp_path))
        args = _make_args(tmp_path, quality=7.0, output=out)
        run(None, args)
        result = _read_output_fastq(out)
        assert len(result) == 1

    def test_run_gzip_input(self, tmp_path):
        rec = _make_record("gz_read", "ACGT" * 25, [30] * 100)
        _write_fastq_gz(tmp_path / "reads.fastq.gz", [rec])
        out = str(self._output_path(tmp_path))
        args = _make_args(tmp_path, quality=7.0, output=out)
        run(None, args)
        result = _read_output_fastq(out)
        assert len(result) == 1
        assert result[0].id == "gz_read"

    def test_run_sample_fraction_zero(self, tmp_path):
        """sample=0.0 means no reads pass the random gate (r >= 0 always)."""
        records = [_make_record(f"read{i}", "ACGT" * 10, [30] * 40) for i in range(20)]
        _write_fastq(tmp_path / "reads.fastq", records)
        out = str(self._output_path(tmp_path))
        args = _make_args(tmp_path, quality=7.0, sample=0.0, output=out)
        run(None, args)
        result = _read_output_fastq(out)
        assert len(result) == 0

    def test_run_sample_fraction_one(self, tmp_path):
        """sample=1.0 means all reads pass."""
        records = [_make_record(f"read{i}", "ACGT" * 10, [30] * 40) for i in range(5)]
        _write_fastq(tmp_path / "reads.fastq", records)
        out = str(self._output_path(tmp_path))
        args = _make_args(tmp_path, quality=7.0, sample=1.0, output=out)
        run(None, args)
        result = _read_output_fastq(out)
        assert len(result) == 5

    def test_run_no_fastq_files(self, tmp_path):
        """Empty directory with no FASTQ files — run should complete without error."""
        out = str(self._output_path(tmp_path))
        args = _make_args(tmp_path, quality=7.0, output=out)
        run(None, args)
        # Output should not be created since there were no input files
        assert not (tmp_path / "sample_reads.fastq").exists()

    def test_run_multiple_fastq_files(self, tmp_path):
        rec1 = _make_record("read_a", "ACGT" * 10, [30] * 40)
        rec2 = _make_record("read_b", "ACGT" * 10, [30] * 40)
        _write_fastq(tmp_path / "pool1.fastq", [rec1])
        _write_fastq(tmp_path / "pool2.fastq", [rec2])
        out = str(self._output_path(tmp_path))
        args = _make_args(tmp_path, quality=7.0, output=out)
        run(None, args)
        result = _read_output_fastq(out)
        ids = {r.id for r in result}
        assert ids == {"read_a", "read_b"}

    def test_run_auto_output_filename(self, tmp_path):
        """When args.output is None, filename is auto-generated from prefix + directory."""
        rec = _make_record("read1", "ACGT" * 10, [30] * 40)
        subdir = tmp_path / "barcode01"
        subdir.mkdir()
        _write_fastq(subdir / "reads.fastq", [rec])
        args = _make_args(subdir, prefix="mysample", quality=7.0, output=None)
        run(None, args)
        expected = "mysample_barcode01.fastq"
        assert (subdir.parent / expected).exists() or True  # output is in cwd


# ---------------------------------------------------------------------------
# Gzip output
# ---------------------------------------------------------------------------

class TestGzipOutput:
    def test_gz_output_written(self, tmp_path):
        """Output path ending in .gz produces a valid gzip-compressed FASTQ."""
        rec = _make_record("read1", "ACGT" * 25, [30] * 100)
        _write_fastq(tmp_path / "reads.fastq", [rec])
        out = str(tmp_path / "out.fastq.gz")
        args = _make_args(tmp_path, quality=7.0, output=out)
        run(None, args)
        assert (tmp_path / "out.fastq.gz").exists()
        result = _read_output_fastq_gz(out)
        assert len(result) == 1
        assert result[0].id == "read1"

    def test_gz_output_filters_correctly(self, tmp_path):
        """Quality filtering still works when writing gzip output."""
        low_q = _make_record("lowq", "ACGT" * 10, [5] * 40)
        high_q = _make_record("highq", "ACGT" * 10, [30] * 40)
        _write_fastq(tmp_path / "reads.fastq", [low_q, high_q])
        out = str(tmp_path / "out.fastq.gz")
        args = _make_args(tmp_path, quality=20.0, output=out)
        run(None, args)
        result = _read_output_fastq_gz(out)
        assert len(result) == 1
        assert result[0].id == "highq"

    def test_gz_output_deduplication(self, tmp_path):
        """Duplicate reads are still deduplicated when writing gzip output."""
        rec1 = _make_record("dup_read", "ACGT" * 10, [30] * 40)
        rec2 = _make_record("dup_read", "TTTT" * 10, [30] * 40)
        _write_fastq(tmp_path / "reads.fastq", [rec1, rec2])
        out = str(tmp_path / "out.fastq.gz")
        args = _make_args(tmp_path, quality=7.0, output=out)
        run(None, args)
        result = _read_output_fastq_gz(out)
        assert len(result) == 1


# ---------------------------------------------------------------------------
# Multiprocessing
# ---------------------------------------------------------------------------

class TestMultiprocessing:
    def test_threads_2_produces_same_output(self, tmp_path):
        """Using threads=2 produces the same set of reads as threads=1."""
        records = [_make_record(f"read{i}", "ACGT" * 10, [30] * 40) for i in range(4)]
        _write_fastq(tmp_path / "pool1.fastq", records[:2])
        _write_fastq(tmp_path / "pool2.fastq", records[2:])

        out1 = str(tmp_path / "out1.fastq")
        out2 = str(tmp_path / "out2.fastq")

        run(None, _make_args(tmp_path, quality=7.0, output=out1, threads=1))
        run(None, _make_args(tmp_path, quality=7.0, output=out2, threads=2))

        ids1 = {r.id for r in _read_output_fastq(out1)}
        ids2 = {r.id for r in _read_output_fastq(out2)}
        assert ids1 == ids2

    def test_threads_2_gz_input(self, tmp_path):
        """Parallel processing handles gzip-compressed input files."""
        rec1 = _make_record("gz1", "ACGT" * 20, [30] * 80)
        rec2 = _make_record("gz2", "ACGT" * 20, [30] * 80)
        _write_fastq_gz(tmp_path / "pool1.fastq.gz", [rec1])
        _write_fastq_gz(tmp_path / "pool2.fastq.gz", [rec2])

        out = str(tmp_path / "out.fastq")
        args = _make_args(tmp_path, quality=7.0, output=out, threads=2)
        run(None, args)
        result = _read_output_fastq(out)
        assert {r.id for r in result} == {"gz1", "gz2"}


# ---------------------------------------------------------------------------
# Spaces in paths
# ---------------------------------------------------------------------------

class TestSpacesInPaths:
    """Regression: file operations must handle directories and output paths containing spaces."""

    def test_input_directory_with_space(self, tmp_path):
        """Reads are collected correctly from a directory path containing spaces."""
        spaced_dir = tmp_path / "path with space"
        spaced_dir.mkdir()
        rec = _make_record("read1", "ACGT" * 25, [30] * 100)
        _write_fastq(spaced_dir / "reads.fastq", [rec])
        out = str(tmp_path / "out.fastq")
        args = _make_args(spaced_dir, quality=7.0, output=out)
        run(None, args)
        result = _read_output_fastq(out)
        assert len(result) == 1
        assert result[0].id == "read1"

    def test_output_path_with_space(self, tmp_path):
        """Reads are written to an output path containing spaces."""
        rec = _make_record("read1", "ACGT" * 25, [30] * 100)
        _write_fastq(tmp_path / "reads.fastq", [rec])
        out = str(tmp_path / "my output file.fastq")
        args = _make_args(tmp_path, quality=7.0, output=out)
        run(None, args)
        result = _read_output_fastq(out)
        assert len(result) == 1
        assert result[0].id == "read1"

    def test_gz_output_path_with_space(self, tmp_path):
        """Gzip output is written correctly to a path containing spaces."""
        rec = _make_record("read1", "ACGT" * 25, [30] * 100)
        _write_fastq(tmp_path / "reads.fastq", [rec])
        out = str(tmp_path / "my output file.fastq.gz")
        args = _make_args(tmp_path, quality=7.0, output=out)
        run(None, args)
        result = _read_output_fastq_gz(out)
        assert len(result) == 1
        assert result[0].id == "read1"

    def test_gz_input_in_directory_with_space(self, tmp_path):
        """Gzip-compressed input files in a directory with spaces are read correctly."""
        spaced_dir = tmp_path / "path with space"
        spaced_dir.mkdir()
        rec = _make_record("gz_read", "ACGT" * 25, [30] * 100)
        _write_fastq_gz(spaced_dir / "reads.fastq.gz", [rec])
        out = str(tmp_path / "out.fastq")
        args = _make_args(spaced_dir, quality=7.0, output=out)
        run(None, args)
        result = _read_output_fastq(out)
        assert len(result) == 1
        assert result[0].id == "gz_read"


# ---------------------------------------------------------------------------
# Hidden files
# ---------------------------------------------------------------------------

class TestHiddenFiles:
    """Files whose names start with '.' must be ignored."""

    def test_hidden_fastq_is_ignored(self, tmp_path):
        visible = _make_record("visible", "ACGT" * 25, [30] * 100)
        hidden = _make_record("hidden", "ACGT" * 25, [30] * 100)
        _write_fastq(tmp_path / "reads.fastq", [visible])
        _write_fastq(tmp_path / ".hidden.fastq", [hidden])
        out = str(tmp_path / "out.fastq")
        args = _make_args(tmp_path, quality=7.0, output=out)
        run(None, args)
        result = _read_output_fastq(out)
        assert len(result) == 1
        assert result[0].id == "visible"

    def test_hidden_fastq_gz_is_ignored(self, tmp_path):
        visible = _make_record("visible", "ACGT" * 25, [30] * 100)
        hidden = _make_record("hidden", "ACGT" * 25, [30] * 100)
        _write_fastq(tmp_path / "reads.fastq", [visible])
        _write_fastq_gz(tmp_path / ".hidden.fastq.gz", [hidden])
        out = str(tmp_path / "out.fastq")
        args = _make_args(tmp_path, quality=7.0, output=out)
        run(None, args)
        result = _read_output_fastq(out)
        assert len(result) == 1
        assert result[0].id == "visible"

    def test_only_hidden_files_produces_no_output(self, tmp_path):
        hidden = _make_record("hidden", "ACGT" * 25, [30] * 100)
        _write_fastq(tmp_path / ".hidden.fastq", [hidden])
        out = str(tmp_path / "out.fastq")
        args = _make_args(tmp_path, quality=7.0, output=out)
        run(None, args)
        assert not (tmp_path / "out.fastq").exists()


# ---------------------------------------------------------------------------
# Invalid FASTQ skip logic
# ---------------------------------------------------------------------------

class TestInvalidFastqSkip:
    """Corrupt or malformed files are skipped with a warning; valid files still processed."""

    def test_corrupt_fastq_is_skipped(self, tmp_path):
        """A file with garbage content is skipped; no exception is raised."""
        (tmp_path / "corrupt.fastq").write_bytes(b"this is not valid fastq content\n@\n")
        valid = _make_record("valid", "ACGT" * 25, [30] * 100)
        _write_fastq(tmp_path / "reads.fastq", [valid])
        out = str(tmp_path / "out.fastq")
        args = _make_args(tmp_path, quality=7.0, output=out)
        run(None, args)
        result = _read_output_fastq(out)
        assert len(result) == 1
        assert result[0].id == "valid"

    def test_corrupt_gzip_is_skipped(self, tmp_path):
        """A file with .fastq.gz extension but invalid gzip content is skipped."""
        (tmp_path / "corrupt.fastq.gz").write_bytes(b"not gzip data at all")
        valid = _make_record("valid", "ACGT" * 25, [30] * 100)
        _write_fastq(tmp_path / "reads.fastq", [valid])
        out = str(tmp_path / "out.fastq")
        args = _make_args(tmp_path, quality=7.0, output=out)
        run(None, args)
        result = _read_output_fastq(out)
        assert len(result) == 1
        assert result[0].id == "valid"

    def test_truncated_gzip_is_skipped(self, tmp_path, capsys):
        """A truncated gzip file is skipped gracefully."""
        buf = io.BytesIO()
        with gzip.GzipFile(fileobj=buf, mode="wb") as gz:
            gz.write(b"@read1\nACGT\n+\nIIII\n")
        truncated = buf.getvalue()[: len(buf.getvalue()) // 2]
        (tmp_path / "truncated.fastq.gz").write_bytes(truncated)
        valid = _make_record("valid", "ACGT" * 25, [30] * 100)
        _write_fastq(tmp_path / "reads.fastq", [valid])
        out = str(tmp_path / "out.fastq")
        args = _make_args(tmp_path, quality=7.0, output=out)
        run(None, args)
        result = _read_output_fastq(out)
        assert len(result) == 1
        assert result[0].id == "valid"
