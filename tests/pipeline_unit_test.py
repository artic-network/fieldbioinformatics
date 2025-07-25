# pipeline_unit_test.py contains a test for the pipeline parser
import argparse
import pytest

from artic import pipeline


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
