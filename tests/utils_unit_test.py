"""Unit tests for artic/utils.py"""

import gzip
import hashlib
import json
import os
from io import StringIO
from unittest import mock
from unittest.mock import MagicMock, patch, call

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from artic.utils import (
    getPrimerDirection,
    merge_sites,
    overlaps,
    check_hash,
    identify_bed_file,
    choose_model,
    get_scheme,
)


# ---------------------------------------------------------------------------
# getPrimerDirection
# ---------------------------------------------------------------------------

class TestGetPrimerDirection:
    def test_left_primer(self):
        assert getPrimerDirection("nCoV-2019_1_LEFT") == "+"

    def test_right_primer(self):
        assert getPrimerDirection("nCoV-2019_1_RIGHT") == "-"

    def test_left_in_compound_name(self):
        assert getPrimerDirection("scheme_5_LEFT_alt1") == "+"

    def test_right_in_compound_name(self):
        assert getPrimerDirection("SARS-CoV-2_10_RIGHT_alt2") == "-"

    def test_invalid_primer_raises(self):
        with pytest.raises(SystemExit) as exc:
            getPrimerDirection("nCoV-2019_1_FORWARD")
        assert exc.value.code == 1

    def test_empty_string_raises(self):
        with pytest.raises(SystemExit) as exc:
            getPrimerDirection("")
        assert exc.value.code == 1


# ---------------------------------------------------------------------------
# merge_sites
# ---------------------------------------------------------------------------

def _make_site(start, end, direction="+"):
    return {"start": start, "end": end, "direction": direction}


class TestMergeSites:
    def test_merge_expands_start(self):
        canonical = _make_site(50, 100)
        alt = _make_site(40, 100)
        result = merge_sites(canonical, alt)
        assert result["start"] == 40

    def test_merge_expands_end(self):
        canonical = _make_site(50, 100)
        alt = _make_site(50, 110)
        result = merge_sites(canonical, alt)
        assert result["end"] == 110

    def test_merge_no_expansion(self):
        canonical = _make_site(30, 120)
        alt = _make_site(50, 100)
        result = merge_sites(canonical, alt)
        assert result["start"] == 30
        assert result["end"] == 120

    def test_merge_both_directions_expand(self):
        canonical = _make_site(50, 100)
        alt = _make_site(40, 110)
        result = merge_sites(canonical, alt)
        assert result["start"] == 40
        assert result["end"] == 110

    def test_merge_direction_mismatch_raises(self):
        canonical = _make_site(50, 100, direction="+")
        alt = _make_site(40, 100, direction="-")
        with pytest.raises(SystemExit):
            merge_sites(canonical, alt)

    def test_merge_preserves_other_fields(self):
        canonical = {"start": 50, "end": 100, "direction": "+", "name": "primer1"}
        alt = _make_site(40, 100)
        result = merge_sites(canonical, alt)
        assert result["name"] == "primer1"


# ---------------------------------------------------------------------------
# overlaps
# ---------------------------------------------------------------------------

class TestOverlaps:
    def _coords(self):
        return [
            {"start": 10, "end": 50, "name": "region1"},
            {"start": 100, "end": 200, "name": "region2"},
        ]

    def test_overlaps_hit_middle(self):
        result = overlaps(self._coords(), 30)
        assert result["name"] == "region1"

    def test_overlaps_miss(self):
        assert overlaps(self._coords(), 60) is False

    def test_overlaps_boundary_start(self):
        result = overlaps(self._coords(), 10)
        assert result["name"] == "region1"

    def test_overlaps_boundary_end(self):
        result = overlaps(self._coords(), 50)
        assert result["name"] == "region1"

    def test_overlaps_second_region(self):
        result = overlaps(self._coords(), 150)
        assert result["name"] == "region2"

    def test_overlaps_empty_coords(self):
        assert overlaps([], 30) is False

    def test_overlaps_before_all_regions(self):
        assert overlaps(self._coords(), 0) is False

    def test_overlaps_after_all_regions(self):
        assert overlaps(self._coords(), 999) is False


# ---------------------------------------------------------------------------
# check_hash
# ---------------------------------------------------------------------------

class TestCheckHash:
    def test_check_hash_match(self, tmp_path):
        content = b"hello world"
        f = tmp_path / "file.txt"
        f.write_bytes(content)
        correct_md5 = hashlib.md5(content).hexdigest()
        # Should not raise
        check_hash(str(f), correct_md5)

    def test_check_hash_mismatch_raises(self, tmp_path):
        content = b"hello world"
        f = tmp_path / "file.txt"
        f.write_bytes(content)
        with pytest.raises(SystemExit) as exc:
            check_hash(str(f), "0" * 32)
        assert exc.value.code == 1

    def test_check_hash_empty_file(self, tmp_path):
        f = tmp_path / "empty.txt"
        f.write_bytes(b"")
        correct_md5 = hashlib.md5(b"").hexdigest()
        check_hash(str(f), correct_md5)


# ---------------------------------------------------------------------------
# identify_bed_file
# ---------------------------------------------------------------------------

class TestIdentifyBedFile:
    def test_identify_v1(self, tmp_bed_v1):
        assert identify_bed_file(tmp_bed_v1) == 1

    def test_identify_v2(self, tmp_bed_v2):
        assert identify_bed_file(tmp_bed_v2) == 2

    def test_identify_v3(self, tmp_bed_v3):
        assert identify_bed_file(tmp_bed_v3) == 3

    def test_identify_too_few_cols_raises(self, tmp_path):
        bed = tmp_path / "bad.bed"
        bed.write_text("MN908947.3\t30\t54\n")
        with pytest.raises(SystemExit) as exc:
            identify_bed_file(str(bed))
        assert exc.value.code == 1

    def test_identify_invalid_primer_name_raises(self, tmp_path):
        bed = tmp_path / "bad_name.bed"
        bed.write_text("MN908947.3\t30\t54\tBADNAME\t1\t+\n")
        with pytest.raises(SystemExit) as exc:
            identify_bed_file(str(bed))
        assert exc.value.code == 1

    def test_identify_only_reads_first_non_comment_line(self, tmp_path):
        """Function returns after parsing the first non-comment line."""
        bed = tmp_path / "two_lines.bed"
        # First line is V1 (6-col); second is V3 — only first line is read
        bed.write_text(
            "MN908947.3\t30\t54\tnCoV-2019_1_LEFT\t1\t+\n"
            "MN908947.3\t385\t410\tnCoV-2019_1_RIGHT_0\t1\t-\t1\n"
        )
        assert identify_bed_file(str(bed)) == 1

    def test_identify_skips_comment_lines(self, tmp_path):
        bed = tmp_path / "commented.bed"
        bed.write_text(
            "# this is a comment\n"
            "MN908947.3\t30\t54\tnCoV-2019_1_LEFT\t1\t+\n"
        )
        # Comments are skipped; first real line is V1
        assert identify_bed_file(str(bed)) == 1


# ---------------------------------------------------------------------------
# choose_model
# ---------------------------------------------------------------------------

def _write_fastq(path, read_id, description, seq="ACGT" * 10, qual=30):
    record = SeqRecord(
        Seq(seq),
        id=read_id,
        description=description,
        letter_annotations={"phred_quality": [qual] * len(seq)},
    )
    with open(path, "w") as fh:
        SeqIO.write([record], fh, "fastq")


class TestChooseModel:
    def test_choose_model_4part_hac(self, tmp_path):
        """4-segment basecall_model_version_id with hac preset → r941_prom_hac_g360+g422"""
        fastq = str(tmp_path / "reads.fastq")
        _write_fastq(
            fastq,
            "read1",
            "read1 basecall_model_version_id=dna_r9.4.1_e8_hac@v3.3",
        )
        result = choose_model(fastq)
        assert result["name"] == "r941_prom_hac_g360+g422"

    def test_choose_model_4part_sup(self, tmp_path):
        """4-segment with sup preset → r941_prom_sup_g5014"""
        fastq = str(tmp_path / "reads.fastq")
        _write_fastq(
            fastq,
            "read1",
            "read1 basecall_model_version_id=dna_r9.4.1_e8_sup@v3.3",
        )
        result = choose_model(fastq)
        assert result["name"] == "r941_prom_sup_g5014"

    def test_choose_model_5part(self, tmp_path):
        """5-segment ID with hac@v4.3.0 → r1041_e82_400bps_hac_v430"""
        fastq = str(tmp_path / "reads.fastq")
        _write_fastq(
            fastq,
            "read1",
            "read1 basecall_model_version_id=dna_r10.4.1_e8.2_400bps_hac@v4.3.0",
        )
        result = choose_model(fastq)
        assert result["name"] == "r1041_e82_400bps_hac_v430"

    def test_choose_model_gz(self, tmp_path):
        """Reads a gzipped FASTQ file correctly."""
        import io
        fastq_gz = str(tmp_path / "reads.fastq.gz")
        record = SeqRecord(
            Seq("ACGT" * 10),
            id="read1",
            description="read1 basecall_model_version_id=dna_r9.4.1_e8_hac@v3.3",
            letter_annotations={"phred_quality": [30] * 40},
        )
        buf = io.StringIO()
        SeqIO.write([record], buf, "fastq")
        with gzip.open(fastq_gz, "wt") as fh:
            fh.write(buf.getvalue())
        result = choose_model(fastq_gz)
        assert result["name"] == "r941_prom_hac_g360+g422"

    def test_choose_model_no_tag_exits(self, tmp_path):
        """Missing basecall_model_version_id in header → sys.exit(6)"""
        fastq = str(tmp_path / "reads.fastq")
        _write_fastq(fastq, "read1", "read1 some_other_tag=value")
        with pytest.raises(SystemExit) as exc:
            choose_model(fastq)
        assert exc.value.code == 6

    def test_choose_model_unknown_pore_exits(self, tmp_path):
        """Unrecognised pore type → sys.exit(6)"""
        fastq = str(tmp_path / "reads.fastq")
        _write_fastq(
            fastq,
            "read1",
            "read1 basecall_model_version_id=dna_rXXXX_e8_hac@v1.0",
        )
        with pytest.raises(SystemExit) as exc:
            choose_model(fastq)
        assert exc.value.code == 6

    def test_choose_model_empty_fastq_exits(self, tmp_path):
        """Empty FASTQ file → sys.exit(4)"""
        fastq = str(tmp_path / "empty.fastq")
        open(fastq, "w").close()
        with pytest.raises(SystemExit) as exc:
            choose_model(fastq)
        assert exc.value.code == 4

    def test_choose_model_empty_fastq_gz_exits(self, tmp_path):
        """Empty gzipped FASTQ file → sys.exit(4)"""
        import gzip as _gzip
        fastq_gz = str(tmp_path / "empty.fastq.gz")
        with _gzip.open(fastq_gz, "wt"):
            pass
        with pytest.raises(SystemExit) as exc:
            choose_model(fastq_gz)
        assert exc.value.code == 4


# ---------------------------------------------------------------------------
# choose_model — all Dorado DNA models from the official list
# https://software-docs.nanoporetech.com/dorado/latest/models/list/
# ---------------------------------------------------------------------------

# Dorado DNA models with an exact-version Clair3 counterpart on HKU.
_DORADO_DNA_EXACT = [
    ("dna_r10.4.1_e8.2_400bps_hac@v4.2.0", "r1041_e82_400bps_hac_v420"),
    ("dna_r10.4.1_e8.2_400bps_sup@v4.2.0", "r1041_e82_400bps_sup_v420"),
    ("dna_r10.4.1_e8.2_400bps_hac@v4.3.0", "r1041_e82_400bps_hac_v430"),
    ("dna_r10.4.1_e8.2_400bps_sup@v4.3.0", "r1041_e82_400bps_sup_v430"),
    ("dna_r10.4.1_e8.2_400bps_hac@v5.0.0", "r1041_e82_400bps_hac_v500"),
    ("dna_r10.4.1_e8.2_400bps_sup@v5.0.0", "r1041_e82_400bps_sup_v500"),
    ("dna_r10.4.1_e8.2_400bps_hac@v5.2.0", "r1041_e82_400bps_hac_v520"),
    ("dna_r10.4.1_e8.2_400bps_sup@v5.2.0", "r1041_e82_400bps_sup_v520"),
]

# Dorado fast models have no versioned (v-prefix) Clair3 counterpart on HKU —
# only Guppy-era (g-prefix) fast models exist, which must not be used silently
# with Dorado data as they can produce severely incorrect variant calls.
_DORADO_DNA_FAST_NO_MODEL = [
    "dna_r10.4.1_e8.2_400bps_fast@v4.2.0",
    "dna_r10.4.1_e8.2_400bps_fast@v4.3.0",
    "dna_r10.4.1_e8.2_400bps_fast@v5.0.0",
    "dna_r10.4.1_e8.2_400bps_fast@v5.2.0",
]


@pytest.mark.parametrize("dorado_model,expected_clair3", _DORADO_DNA_EXACT)
def test_choose_model_dorado_dna_exact(tmp_path, dorado_model, expected_clair3):
    """choose_model maps Dorado hac/sup models to their exact Clair3 counterpart."""
    fastq = str(tmp_path / "reads.fastq")
    _write_fastq(fastq, "read1", f"read1 basecall_model_version_id={dorado_model}")
    result = choose_model(fastq)
    assert result["name"] == expected_clair3, (
        f"Dorado model {dorado_model!r} → expected Clair3 {expected_clair3!r}, got {result['name']!r}"
    )


@pytest.mark.parametrize("dorado_model", _DORADO_DNA_FAST_NO_MODEL)
def test_choose_model_dorado_fast_no_versioned_model_exits(tmp_path, dorado_model):
    """Dorado fast models have no versioned Clair3 counterpart — must exit(6) rather
    than silently falling back to an incompatible Guppy-era model."""
    fastq = str(tmp_path / "reads.fastq")
    _write_fastq(fastq, "read1", f"read1 basecall_model_version_id={dorado_model}")
    with pytest.raises(SystemExit) as exc:
        choose_model(fastq)
    assert exc.value.code == 6


# ---------------------------------------------------------------------------
# get_scheme (mocked network)
# ---------------------------------------------------------------------------

def _make_manifest(scheme_name="test-scheme", length="400", version="v1.0.0"):
    """Return a minimal manifest dict mirroring the index.json structure."""
    return {
        "primerschemes": {
            scheme_name: {
                length: {
                    version: {
                        "primer_bed_url": "https://example.com/primer.bed",
                        "primer_bed_md5": "abc123",
                        "reference_fasta_url": "https://example.com/reference.fasta",
                        "reference_fasta_md5": "def456",
                    }
                }
            }
        }
    }


def _make_aliases(alias, real_name):
    return {alias: real_name}


class TestGetScheme:
    def _mock_response(self, json_data=None, text="file content", status_code=200):
        resp = MagicMock()
        resp.json.return_value = json_data or {}
        resp.text = text
        resp.status_code = status_code
        return resp

    @patch("artic.utils.check_hash")
    @patch("artic.utils.requests.get")
    def test_get_scheme_success(self, mock_get, mock_hash, tmp_path):
        """Happy path: scheme found, files downloaded, paths returned."""
        manifest = _make_manifest("test-scheme", "400", "v1.0.0")
        aliases = {}
        mock_get.side_effect = [
            self._mock_response(manifest),   # index.json
            self._mock_response(aliases),    # aliases.json
            self._mock_response(text="bed"),  # primer bed download
            self._mock_response(text="ref", status_code=200),  # ref fasta download
        ]
        mock_hash.return_value = None

        bed, ref, ver = get_scheme(
            scheme_name="test-scheme",
            scheme_version="v1.0.0",
            scheme_directory=str(tmp_path),
        )

        assert bed.endswith("primer.bed")
        assert ref.endswith("reference.fasta")
        assert ver == "v1.0.0"

    @patch("artic.utils.check_hash")
    @patch("artic.utils.requests.get")
    def test_get_scheme_alias_resolution(self, mock_get, mock_hash, tmp_path):
        """scheme_name found via alias → resolves to real scheme name."""
        manifest = _make_manifest("real-scheme", "400", "v1.0.0")
        aliases = {"alias-scheme": "real-scheme"}
        mock_get.side_effect = [
            self._mock_response(manifest),
            self._mock_response(aliases),
            self._mock_response(text="bed"),
            self._mock_response(text="ref", status_code=200),
        ]
        mock_hash.return_value = None

        bed, ref, ver = get_scheme(
            scheme_name="alias-scheme",
            scheme_version="v1.0.0",
            scheme_directory=str(tmp_path),
        )
        assert ver == "v1.0.0"

    @patch("artic.utils.requests.get")
    def test_get_scheme_unknown_scheme_exits(self, mock_get, tmp_path):
        """Scheme not in manifest or aliases → SystemExit(1)."""
        manifest = _make_manifest("some-other-scheme")
        aliases = {}
        mock_get.side_effect = [
            self._mock_response(manifest),
            self._mock_response(aliases),
        ]
        with pytest.raises(SystemExit) as exc:
            get_scheme(
                scheme_name="nonexistent-scheme",
                scheme_version="v1.0.0",
                scheme_directory=str(tmp_path),
            )
        assert exc.value.code == 1

    @patch("artic.utils.requests.get")
    def test_get_scheme_invalid_version_format_exits(self, mock_get, tmp_path):
        """Version not matching vX.X.X format → SystemExit(1)."""
        manifest = _make_manifest("test-scheme")
        aliases = {}
        mock_get.side_effect = [
            self._mock_response(manifest),
            self._mock_response(aliases),
        ]
        with pytest.raises(SystemExit) as exc:
            get_scheme(
                scheme_name="test-scheme",
                scheme_version="1.0",  # missing 'v' prefix and third segment
                scheme_directory=str(tmp_path),
            )
        assert exc.value.code == 1

    @patch("artic.utils.requests.get")
    def test_get_scheme_version_not_found_exits(self, mock_get, tmp_path):
        """Requested version not in manifest → SystemExit(1)."""
        manifest = _make_manifest("test-scheme", "400", "v1.0.0")
        aliases = {}
        mock_get.side_effect = [
            self._mock_response(manifest),
            self._mock_response(aliases),
        ]
        with pytest.raises(SystemExit) as exc:
            get_scheme(
                scheme_name="test-scheme",
                scheme_version="v9.9.9",
                scheme_directory=str(tmp_path),
            )
        assert exc.value.code == 1

    @patch("artic.utils.requests.get")
    def test_get_scheme_network_failure_no_cache_exits(self, mock_get, tmp_path):
        """Network down, no local cache → SystemExit(1)."""
        mock_get.side_effect = Exception("network error")
        with pytest.raises(SystemExit) as exc:
            get_scheme(
                scheme_name="test-scheme",
                scheme_version="v1.0.0",
                scheme_directory=str(tmp_path),
            )
        assert exc.value.code == 1

    @patch("artic.utils.check_hash")
    @patch("artic.utils.requests.get")
    def test_get_scheme_network_failure_with_local_cache(
        self, mock_get, mock_hash, tmp_path
    ):
        """Network down but local index.json + aliases.json present → uses cache."""
        manifest = _make_manifest("test-scheme", "400", "v1.0.0")
        aliases = {}

        # Write the local cache files
        index_path = tmp_path / "index.json"
        aliases_path = tmp_path / "aliases.json"
        index_path.write_text(json.dumps(manifest))
        aliases_path.write_text(json.dumps(aliases))

        # First two requests (manifest + aliases) fail; rest succeed for downloads
        def _side_effect(url, *args, **kwargs):
            if "index.json" in url or "aliases.json" in url:
                raise Exception("network error")
            return self._mock_response(text="content", status_code=200)

        mock_get.side_effect = _side_effect
        mock_hash.return_value = None

        bed, ref, ver = get_scheme(
            scheme_name="test-scheme",
            scheme_version="v1.0.0",
            scheme_directory=str(tmp_path),
        )
        assert ver == "v1.0.0"

    @patch("artic.utils.requests.get")
    def test_get_scheme_multiple_lengths_no_length_arg_exits(self, mock_get, tmp_path):
        """Scheme has multiple lengths and none specified → SystemExit(1)."""
        manifest = {
            "primerschemes": {
                "multi-len-scheme": {
                    "400": {"v1.0.0": {}},
                    "800": {"v1.0.0": {}},
                }
            }
        }
        aliases = {}
        mock_get.side_effect = [
            self._mock_response(manifest),
            self._mock_response(aliases),
        ]
        with pytest.raises(SystemExit) as exc:
            get_scheme(
                scheme_name="multi-len-scheme",
                scheme_version="v1.0.0",
                scheme_directory=str(tmp_path),
            )
        assert exc.value.code == 1
