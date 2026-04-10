"""Unit tests for artic/mask.py"""

from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import pytest

from artic.mask import read_3col_bed, go


# ---------------------------------------------------------------------------
# read_3col_bed
# ---------------------------------------------------------------------------

class TestRead3ColBed:
    def test_read_valid_bed(self, tmp_path):
        bed = tmp_path / "mask.bed"
        bed.write_text("chrom1\t10\t50\nchrom1\t100\t200\n")
        result = read_3col_bed(str(bed))
        assert len(result) == 2
        assert result[0] == {"chrom": "chrom1", "start": 10, "end": 50}
        assert result[1] == {"chrom": "chrom1", "start": 100, "end": 200}

    def test_read_positions_are_integers(self, tmp_path):
        bed = tmp_path / "mask.bed"
        bed.write_text("MN908947.3\t30\t60\n")
        result = read_3col_bed(str(bed))
        assert isinstance(result[0]["start"], int)
        assert isinstance(result[0]["end"], int)

    def test_read_nonnumeric_start_raises(self, tmp_path):
        bed = tmp_path / "bad.bed"
        bed.write_text("chrom1\tabc\t50\n")
        with pytest.raises(ValueError, match="malformed"):
            read_3col_bed(str(bed))

    def test_read_nonnumeric_end_raises(self, tmp_path):
        bed = tmp_path / "bad.bed"
        bed.write_text("chrom1\t10\txyz\n")
        with pytest.raises(ValueError, match="malformed"):
            read_3col_bed(str(bed))

    def test_read_empty_file(self, tmp_path):
        bed = tmp_path / "empty.bed"
        bed.write_text("")
        result = read_3col_bed(str(bed))
        assert result == []

    def test_read_multiple_rows(self, tmp_path):
        bed = tmp_path / "multi.bed"
        bed.write_text(
            "chrom1\t0\t10\n"
            "chrom1\t20\t30\n"
            "chrom1\t40\t50\n"
        )
        result = read_3col_bed(str(bed))
        assert len(result) == 3

    def test_read_preserves_chrom_name(self, tmp_path):
        bed = tmp_path / "chrom.bed"
        bed.write_text("MN908947.3\t5\t15\n")
        result = read_3col_bed(str(bed))
        assert result[0]["chrom"] == "MN908947.3"


# ---------------------------------------------------------------------------
# go  (mask application)
# ---------------------------------------------------------------------------

def _write_fasta(path, records):
    """Write {id: seq} dict to FASTA."""
    with open(path, "w") as fh:
        for seq_id, seq in records.items():
            fh.write(f">{seq_id}\n{seq}\n")


def _write_bed(path, lines):
    """Write list-of-tuples (chrom, start, end) as 3-col BED."""
    with open(path, "w") as fh:
        for chrom, start, end in lines:
            fh.write(f"{chrom}\t{start}\t{end}\n")


def _make_vcf_record(chrom, pos, ref):
    """pos is 1-based (genomics convention); stored as 0-based to match pysam."""
    v = MagicMock()
    v.chrom = chrom
    v.pos = pos - 1  # pysam is 0-based
    v.ref = ref
    return v


def _mock_vcf_reader(records):
    reader = MagicMock()
    reader.__iter__ = MagicMock(return_value=iter(records))
    reader.__enter__ = MagicMock(return_value=reader)
    reader.__exit__ = MagicMock(return_value=False)
    return reader


class TestGo:
    def _build_args(self, tmp_path, fasta_content, bed_lines, vcf_records=None):
        """Build a SimpleNamespace args and write input files."""
        fasta = tmp_path / "ref.fasta"
        _write_fasta(str(fasta), fasta_content)

        bed = tmp_path / "mask.bed"
        _write_bed(str(bed), bed_lines)

        output = tmp_path / "masked.fasta"

        args = SimpleNamespace(
            reference=str(fasta),
            maskfile=str(bed),
            maskvcf="dummy.vcf",
            output=str(output),
        )
        return args, str(output), vcf_records or []

    @patch("artic.mask.pysam.VariantFile")
    def test_go_masks_bed_regions(self, mock_vcf_cls, tmp_path):
        seq = "A" * 100
        args, output, vcf_records = self._build_args(
            tmp_path,
            {"MN908947.3": seq},
            [("MN908947.3", 10, 20)],
            vcf_records=[],
        )
        mock_vcf_cls.return_value = _mock_vcf_reader([])

        go(args)

        with open(output) as fh:
            lines = fh.read().splitlines()
        consensus = lines[1]
        # Positions 10–19 (0-based) should be N
        assert consensus[10:20] == "N" * 10
        # Positions outside the mask should be A
        assert consensus[0:10] == "A" * 10
        assert consensus[20:30] == "A" * 10

    @patch("artic.mask.pysam.VariantFile")
    def test_go_masks_vcf_variants(self, mock_vcf_cls, tmp_path):
        seq = "A" * 100
        args, output, _ = self._build_args(
            tmp_path,
            {"MN908947.3": seq},
            [],  # no BED masking
        )
        # VCF record at POS=50 (1-based), REF=AAA → masks positions 49,50,51
        vcf_record = _make_vcf_record("MN908947.3", 50, "AAA")
        mock_vcf_cls.return_value = _mock_vcf_reader([vcf_record])

        go(args)

        with open(output) as fh:
            lines = fh.read().splitlines()
        consensus = lines[1]
        assert consensus[49:52] == "N" * 3
        assert consensus[48] == "A"
        assert consensus[52] == "A"

    @patch("artic.mask.pysam.VariantFile")
    def test_go_combined_masking(self, mock_vcf_cls, tmp_path):
        seq = "A" * 100
        args, output, _ = self._build_args(
            tmp_path,
            {"MN908947.3": seq},
            [("MN908947.3", 0, 5)],
        )
        vcf_record = _make_vcf_record("MN908947.3", 80, "A")
        mock_vcf_cls.return_value = _mock_vcf_reader([vcf_record])

        go(args)

        with open(output) as fh:
            lines = fh.read().splitlines()
        consensus = lines[1]
        assert consensus[0:5] == "N" * 5   # BED masked
        assert consensus[79] == "N"          # VCF masked (POS=80, 0-based index 79)
        assert consensus[5] == "A"           # between masks

    @patch("artic.mask.pysam.VariantFile")
    def test_go_output_written(self, mock_vcf_cls, tmp_path):
        seq = "ACGT" * 25  # 100 bases
        args, output, _ = self._build_args(
            tmp_path,
            {"MN908947.3": seq},
            [],
        )
        mock_vcf_cls.return_value = _mock_vcf_reader([])

        go(args)

        assert (tmp_path / "masked.fasta").exists()
        with open(output) as fh:
            content = fh.read()
        assert ">MN908947.3" in content

    @patch("artic.mask.pysam.VariantFile")
    def test_go_multiple_sequences(self, mock_vcf_cls, tmp_path):
        seq_a = "A" * 50
        seq_b = "C" * 50
        args, output, _ = self._build_args(
            tmp_path,
            {"seqA": seq_a, "seqB": seq_b},
            [("seqA", 10, 20)],
        )
        mock_vcf_cls.return_value = _mock_vcf_reader([])

        go(args)

        with open(output) as fh:
            content = fh.read()
        assert ">seqA" in content
        assert ">seqB" in content

    @patch("artic.mask.pysam.VariantFile")
    def test_go_empty_bed_no_masking(self, mock_vcf_cls, tmp_path):
        seq = "ACGT" * 10
        args, output, _ = self._build_args(
            tmp_path,
            {"MN908947.3": seq},
            [],  # empty bed
        )
        mock_vcf_cls.return_value = _mock_vcf_reader([])

        go(args)

        with open(output) as fh:
            lines = fh.read().splitlines()
        assert lines[1] == seq
