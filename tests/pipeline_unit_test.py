# pipeline_unit_test.py contains a test for the pipeline parser
import gzip
import os
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
