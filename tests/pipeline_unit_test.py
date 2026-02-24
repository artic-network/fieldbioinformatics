# pipeline_unit_test.py contains a test for the pipeline parser
import pathlib
import tarfile
import pytest
from unittest import TestCase, mock
import os
import requests
import tqdm
import sys
from io import StringIO
import contextlib
from types import SimpleNamespace

from artic import pipeline, minion


# download will download and untar a test dataset
def download(url, dataDir):
    filename = url.rsplit("/", 1)[1]
    with open(f"{dataDir}/{filename}", "wb+") as f:
        response = requests.get(url, stream=True)
        total = int(response.headers.get("content-length"))
        if total is None:
            f.write(response.content)
        else:
            with tqdm(total=total, unit="B", unit_scale=True, desc=filename) as pbar:
                for data in tqdm(response.iter_content(chunk_size=1024)):
                    f.write(data)
                    pbar.update(1024)

    tar = tarfile.open(dataDir + "/" + filename, "r:gz")
    tar.extractall(dataDir)
    tar.close()
    os.remove(dataDir + "/" + filename)


class TestNonZeroExit(TestCase):
    """Test that the pipeline properly calls an OOM error."""

    def setUp(self):
        testfile_url = "https://raw.githubusercontent.com/artic-network/fieldbioinformatics/master/test-data/MT007544/MT007544.fastq"
        dataDir = str(pathlib.Path(__file__).parent.parent) + "/test-data/"

        targetPath = dataDir + "MT007544"
        if not os.path.exists(targetPath):
            print("\tno data for {}".format("MT007544"))
            print("\tmaking dir at {}".format(targetPath))
            os.mkdir(targetPath)
            print("\tdownloading from {}".format(testfile_url))
            try:
                download(testfile_url, dataDir)
            except Exception as e:
                print("download failed: ", e)
                sys.exit(1)
        else:
            print("\tfound data dir for {}".format("MT007544"))

    print("validation datasets ready\n")

    def test_pipeline_oom(self):
        patcher = mock.patch("artic.minion.subprocess.run")
        subprocess_mock = patcher.start()

        subprocess_mock.return_value.returncode = 137

        parser = False

        args = SimpleNamespace(
            model="r1041_e82_400bps_sup_v420",
            read_file="test-data/MT007544/MT007544.fastq",
            scheme_name="artic-pan-dengue",
            scheme_version="v1.0.0",
            sample="some-prefix",
            threads=1,
            min_depth=20,
            allow_mismatched_primers=False,
            primer_match_threshold=4,
            normalise=200,
            model_dir=None,
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

        stderr = StringIO()
        with pytest.raises(SystemExit) as cm, contextlib.redirect_stderr(stderr):
            minion.run(parser, args)
            print(stderr, file=sys.stderr)

        assert cm.value.code == 137, "expected exit code 137 for OOM error"
        assert (
            "Process ran out of memory, please try running on a machine with more memory available."
            in stderr.getvalue()
        ), "expected OOM error message not found in stderr"


def test_pipeline_parser():
    """basic test for the pipeline parser"""
    # setup a parser
    parser = pipeline.init_pipeline_parser()

    # set up a valid command
    dummyCLI = [
        "minion",
        "--model",
        "r1041_e82_400bps_sup_v420",
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
