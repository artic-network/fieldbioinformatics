# pipeline_unit_test.py contains a test for the pipeline parser
import gzip
import os
import shlex
import tempfile

import pytest
from unittest import TestCase, mock
from types import SimpleNamespace

from artic import pipeline, minion


class TestNonZeroExit(TestCase):
    """Test that the pipeline properly calls an OOM error."""

    def test_pipeline_oom(self):
        mock_scheme = mock.MagicMock()
        mock_scheme.bedlines = [mock.MagicMock(pool="pool1")]

        with mock.patch("artic.minion.subprocess.run") as subprocess_mock, \
             mock.patch("artic.minion.os.path.exists", return_value=True), \
             mock.patch("artic.minion.os.system"), \
             mock.patch("artic.minion.os.remove"), \
             mock.patch("artic.minion.get_scheme", return_value=("fake.bed", "fake.ref", "v1.0.0")), \
             mock.patch("artic.minion.Scheme.from_file", return_value=mock_scheme), \
             mock.patch("builtins.open", mock.mock_open()):

            subprocess_mock.return_value.returncode = 137
            subprocess_mock.return_value.stderr = b""

            parser = False

            args = SimpleNamespace(
                model="r941_prom_hac_g360+g422",
                read_file="test-data/MT007544/MT007544.fastq",
                scheme_name="artic-pan-dengue",
                scheme_version="v1.0.0",
                sample="some-prefix",
                threads=1,
                min_depth=20,
                allow_mismatched_primers=False,
                primer_match_threshold=4,
                normalise=200,
                model_dir="/mock/models",
                bed=False,
                ref=False,
                scheme_length=False,
                scheme_directory=".",
                min_mapq=20,
                no_frameshifts=False,
                no_indels=False,
                linearise_fasta=False,
                align_consensus=False,
                dry_run=False,
                min_variant_quality=10,
                min_allele_frequency=0.6,
                min_frameshift_quality=50,
                min_mask_allele_frequency=0.1,
                min_minor_allele_count=5,
            )

            with pytest.raises(SystemExit) as cm:
                minion.run(parser, args)

            assert cm.value.code == 137, "expected exit code 137 for OOM error"


class TestEmptyReadFile(TestCase):
    """Exit 4 on empty read file, even when --model is supplied (no model-selection code path)."""

    def _args(self, read_file):
        return SimpleNamespace(
            model="r941_prom_hac_g360+g422",
            read_file=read_file,
            scheme_name="artic-pan-dengue",
            scheme_version="v1.0.0",
            sample="some-prefix",
            threads=1,
            min_depth=20,
            allow_mismatched_primers=False,
            primer_match_threshold=4,
            normalise=200,
            model_dir="/mock/models",
            bed=False,
            ref=False,
            scheme_length=False,
            scheme_directory=".",
            min_mapq=20,
            no_frameshifts=False,
            no_indels=False,
            linearise_fasta=False,
            align_consensus=False,
            dry_run=False,
        )

    def test_empty_fastq_exits_4(self):
        with tempfile.NamedTemporaryFile(suffix=".fastq", delete=False) as f:
            empty_fastq = f.name
        try:
            with mock.patch("artic.minion.get_scheme", return_value=("fake.bed", "fake.ref", "v1.0.0")), \
                 mock.patch("artic.minion.os.path.exists", return_value=True), \
                 mock.patch("artic.minion.os.path.getsize", return_value=0):
                with pytest.raises(SystemExit) as cm:
                    minion.run(False, self._args(empty_fastq))
            assert cm.value.code == 4
        finally:
            os.unlink(empty_fastq)

    def test_empty_fastq_gz_exits_4(self):
        with tempfile.NamedTemporaryFile(suffix=".fastq.gz", delete=False) as f:
            empty_gz = f.name
        with gzip.open(empty_gz, "wt"):
            pass
        try:
            with mock.patch("artic.minion.get_scheme", return_value=("fake.bed", "fake.ref", "v1.0.0")), \
                 mock.patch("artic.minion.os.path.exists", return_value=True):
                with pytest.raises(SystemExit) as cm:
                    minion.run(False, self._args(empty_gz))
            assert cm.value.code == 4
        finally:
            os.unlink(empty_gz)


class TestAlignTrimExitCode(TestCase):
    """Test that align_trim exit code 1 is remapped to exit code 2."""

    def test_align_trim_no_reads_exits_2(self):
        mock_scheme = mock.MagicMock()
        mock_scheme.bedlines = [mock.MagicMock(pool="pool1")]

        def subprocess_side_effect(cmd, **__):
            result = mock.MagicMock()
            result.returncode = 1 if cmd.startswith("align_trim") else 0
            result.stderr = b""
            return result

        with mock.patch("artic.minion.subprocess.run", side_effect=subprocess_side_effect), \
             mock.patch("artic.minion.os.path.exists", return_value=True), \
             mock.patch("artic.minion.os.system"), \
             mock.patch("artic.minion.os.remove"), \
             mock.patch("artic.minion.get_scheme", return_value=("fake.bed", "fake.ref", "v1.0.0")), \
             mock.patch("artic.minion.Scheme.from_file", return_value=mock_scheme), \
             mock.patch("builtins.open", mock.mock_open()):

            parser = False

            args = SimpleNamespace(
                model="r941_prom_hac_g360+g422",
                read_file="test-data/MT007544/MT007544.fastq",
                scheme_name="artic-pan-dengue",
                scheme_version="v1.0.0",
                sample="some-prefix",
                threads=1,
                min_depth=20,
                allow_mismatched_primers=False,
                primer_match_threshold=4,
                normalise=200,
                model_dir="/mock/models",
                bed=False,
                ref=False,
                scheme_length=False,
                scheme_directory=".",
                min_mapq=20,
                no_frameshifts=False,
                no_indels=False,
                linearise_fasta=False,
                align_consensus=False,
                dry_run=False,
                min_variant_quality=10,
                min_allele_frequency=0.6,
                min_frameshift_quality=50,
                min_mask_allele_frequency=0.1,
                min_minor_allele_count=5,
            )

            with pytest.raises(SystemExit) as cm:
                minion.run(parser, args)

            assert cm.value.code == 2, "expected exit code 2 when align_trim returns 1 (no reads aligned)"


def test_pipeline_parser():
    """basic test for the pipeline parser"""
    # setup a parser
    parser = pipeline.init_pipeline_parser()

    # set up a valid command
    dummyCLI = [
        "minion",
        "--model",
        "r941_prom_hac_g360+g422",
        "--read-file",
        "some_reads.fastq",
        "some-prefix",
    ]

    # try with required arguments missing
    with pytest.raises(SystemExit):
        _ = parser.parse_args(dummyCLI[0:2])

    # now check the valid command passes
    try:
        args = parser.parse_args(dummyCLI)
    except SystemExit:
        print("failed to parse valid command")
        assert False

    assert args.command == dummyCLI[0], "incorrect subcommand registered"

    # for arg, val in vars(args).items():
    #    print(arg, val)


# ---------------------------------------------------------------------------
# Spaces in paths
# ---------------------------------------------------------------------------

class TestSpacesInPaths:
    """Regression: shell commands must handle paths/sample names containing spaces.

    Uses dry_run=True so no real tools are executed; the generated log is
    inspected to verify every path argument is properly shell-quoted.
    """

    MODEL = "r941_prom_hac_g360+g422"

    def _args(self, sample, read_file, bed, ref, model_dir):
        return SimpleNamespace(
            model=self.MODEL,
            read_file=read_file,
            sample=sample,
            scheme_name=None,
            scheme_version=None,
            scheme_length=None,
            scheme_directory=".",
            bed=bed,
            ref=ref,
            threads=1,
            min_depth=20,
            allow_mismatched_primers=False,
            primer_match_threshold=4,
            normalise=200,
            model_dir=model_dir,
            min_mapq=20,
            no_frameshifts=False,
            no_indels=False,
            linearise_fasta=False,
            align_consensus=False,
            dry_run=True,
            min_variant_quality=10,
            min_allele_frequency=0.6,
            min_frameshift_quality=50,
            min_mask_allele_frequency=0.1,
            min_minor_allele_count=5,
        )

    def _dry_run(self, tmp_path, sample_name="my sample"):
        """Run minion in dry-run mode with space-containing paths; return log lines."""
        spaced_dir = tmp_path / "path with space"
        spaced_dir.mkdir()

        read_file = str(spaced_dir / "my reads.fastq")
        bed = str(spaced_dir / "primer scheme.bed")
        ref = str(spaced_dir / "ref sequence.fasta")
        model_dir = str(spaced_dir / "clair3 models")
        sample = str(tmp_path / sample_name)

        mock_scheme = mock.MagicMock()
        mock_scheme.bedlines = [mock.MagicMock(pool="pool1")]

        with mock.patch("artic.minion.os.path.exists", return_value=True), \
             mock.patch("artic.minion.os.path.getsize", return_value=100), \
             mock.patch("artic.minion.os.remove"), \
             mock.patch("artic.minion.Scheme.from_file", return_value=mock_scheme):
            minion.run(None, self._args(sample, read_file, bed, ref, model_dir))

        log_file = sample + ".minion.log.txt"
        with open(log_file) as f:
            lines = f.readlines()
        return lines, read_file, bed, ref, sample, model_dir

    def test_dry_run_completes_with_spaced_paths(self, tmp_path):
        """Pipeline must complete without error when all paths contain spaces."""
        lines, *_ = self._dry_run(tmp_path)
        assert lines, "Expected non-empty dry-run log"

    def test_minimap2_quotes_read_file_and_ref(self, tmp_path):
        lines, read_file, _, ref, _, _ = self._dry_run(tmp_path)
        cmd = next(line for line in lines if line.startswith("minimap2"))
        assert shlex.quote(read_file) in cmd, f"read_file not quoted in: {cmd!r}"
        assert shlex.quote(ref) in cmd, f"ref not quoted in: {cmd!r}"

    def test_align_trim_quotes_bed(self, tmp_path):
        lines, _, bed, _, _, _ = self._dry_run(tmp_path)
        cmd = next(line for line in lines if line.startswith("align_trim"))
        assert shlex.quote(bed) in cmd, f"bed not quoted in: {cmd!r}"

    def test_samtools_index_quotes_sample_bam(self, tmp_path):
        lines, _, _, _, sample, _ = self._dry_run(tmp_path)
        cmd = next(line for line in lines if line.startswith("samtools index"))
        assert shlex.quote(sample + ".sorted.bam") in cmd, f"sample bam not quoted in: {cmd!r}"

    def test_clair3_quotes_model_path(self, tmp_path):
        lines, _, _, _, _, model_dir = self._dry_run(tmp_path)
        cmd = next(line for line in lines if line.startswith("run_clair3.sh"))
        full_model = f"{model_dir}/{self.MODEL}"
        assert shlex.quote(full_model) in cmd, f"model path not quoted in: {cmd!r}"

    def test_merge_vcf_quotes_sample_and_bed(self, tmp_path):
        lines, _, bed, _, sample, _ = self._dry_run(tmp_path)
        cmd = next(line for line in lines if line.startswith("artic_vcf_merge"))
        assert shlex.quote(sample) in cmd, f"sample not quoted in: {cmd!r}"
        assert shlex.quote(bed) in cmd, f"bed not quoted in: {cmd!r}"
