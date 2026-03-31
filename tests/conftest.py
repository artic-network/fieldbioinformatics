"""Shared pytest fixtures for artic unit tests."""

import gzip
import hashlib
from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# ---------------------------------------------------------------------------
# BED file fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def tmp_bed_v1(tmp_path):
    """Minimal 6-column V1 BED file (e.g. nCoV-2019_1_LEFT)."""
    bed = tmp_path / "v1.bed"
    bed.write_text(
        "MN908947.3\t30\t54\tnCoV-2019_1_LEFT\t1\t+\n"
        "MN908947.3\t385\t410\tnCoV-2019_1_RIGHT\t1\t-\n"
    )
    return str(bed)


@pytest.fixture
def tmp_bed_v2(tmp_path):
    """Minimal 7-column V2 BED file (V1 primer names + pool column)."""
    bed = tmp_path / "v2.bed"
    bed.write_text(
        "MN908947.3\t30\t54\tnCoV-2019_1_LEFT\t1\t+\t1\n"
        "MN908947.3\t385\t410\tnCoV-2019_1_RIGHT\t1\t-\t1\n"
    )
    return str(bed)


@pytest.fixture
def tmp_bed_v3(tmp_path):
    """Minimal V3 BED file (new primer ID format: name_amplicon_direction_index)."""
    bed = tmp_path / "v3.bed"
    bed.write_text(
        "MN908947.3\t30\t54\tnCoV-2019_1_LEFT_0\t1\t+\t1\n"
        "MN908947.3\t385\t410\tnCoV-2019_1_RIGHT_0\t1\t-\t1\n"
    )
    return str(bed)


# ---------------------------------------------------------------------------
# FASTA fixture
# ---------------------------------------------------------------------------

@pytest.fixture
def tmp_fasta(tmp_path):
    """Two-record FASTA file with 100-base sequences."""
    fasta = tmp_path / "ref.fasta"
    seq = "A" * 100
    fasta.write_text(f">MN908947.3\n{seq}\n>seq2\n{seq}\n")
    return str(fasta)


# ---------------------------------------------------------------------------
# FASTQ fixtures
# ---------------------------------------------------------------------------

def _make_fastq_record(read_id: str, description: str, seq: str = "ACGT" * 10, qual: int = 30):
    """Return a BioPython SeqRecord formatted as a FASTQ record string."""
    record = SeqRecord(
        Seq(seq),
        id=read_id,
        description=description,
        letter_annotations={"phred_quality": [qual] * len(seq)},
    )
    return record


def write_fastq(path, records):
    with open(path, "w") as fh:
        SeqIO.write(records, fh, "fastq")


@pytest.fixture
def tmp_fastq(tmp_path):
    """FASTQ file with one read whose header carries basecall_model_version_id."""
    fastq = tmp_path / "reads.fastq"
    rec = _make_fastq_record(
        "read1",
        "read1 basecall_model_version_id=dna_r10.4.1_e8.2_400bps_hac@v4.3.0",
    )
    write_fastq(str(fastq), [rec])
    return str(fastq)


@pytest.fixture
def tmp_fastq_gz(tmp_path):
    """Gzipped FASTQ with the same content as tmp_fastq."""
    fastq_gz = tmp_path / "reads.fastq.gz"
    rec = _make_fastq_record(
        "read1",
        "read1 basecall_model_version_id=dna_r10.4.1_e8.2_400bps_hac@v4.3.0",
    )
    import io
    buf = io.StringIO()
    SeqIO.write([rec], buf, "fastq")
    with gzip.open(str(fastq_gz), "wt") as fh:
        fh.write(buf.getvalue())
    return str(fastq_gz)


@pytest.fixture
def tmp_fastq_no_model_tag(tmp_path):
    """FASTQ whose read header has no basecall_model_version_id."""
    fastq = tmp_path / "reads_no_model.fastq"
    rec = _make_fastq_record("read1", "read1 some_other_tag=value")
    write_fastq(str(fastq), [rec])
    return str(fastq)


# ---------------------------------------------------------------------------
# VCF variant mock helper
# ---------------------------------------------------------------------------

def make_mock_variant(
    chrom="MN908947.3",
    pos=100,
    ref="A",
    alt=None,
    qual=50.0,
    af=0.9,
    dp=100,
    af_raises=False,
    dp_raises=False,
):
    """Return a MagicMock that looks like a cyvcf2 variant."""
    if alt is None:
        alt = ["T"]

    v = MagicMock()
    v.CHROM = chrom
    v.POS = pos
    v.REF = ref
    v.ALT = alt
    v.QUAL = qual

    def _format(field):
        if field == "AF":
            if af_raises:
                raise KeyError("AF")
            return [[af]]
        if field == "DP":
            if dp_raises:
                raise KeyError("DP")
            return [[dp]]
        raise KeyError(field)

    v.format = MagicMock(side_effect=_format)
    return v
