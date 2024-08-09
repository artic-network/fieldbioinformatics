# vcftagprimersites_unit_test.py contains unit tests for the vcf primer site tagging
import os
import pytest

from artic import vcftagprimersites


# help pytest resolve where test data is kept
TEST_DIR = os.path.dirname(os.path.abspath(__file__))


def test_read_bed_file():

    # process the nCoV-2019 V3 primer scheme
    primerScheme = vcftagprimersites.read_bed_file(
        TEST_DIR + "/../test-data/primer-schemes/nCoV-2019/V3/nCoV-2019.scheme.bed"
    )

    # check the the alts have been collapsed into a canonical primer site
    assert len(primerScheme) == 196, "alts were not collapsed"

    # check that alts are merging by union (not intersection)
    for entry in primerScheme:
        if entry["Primer_ID"] == "nCoV-2019_45_RIGHT":
            assert (
                entry["start"] == 13660
            ), "failed to merge nCov-2019_45_RIGHT alt, bad start"
            assert (
                entry["end"] == 13699
            ), "failed to merge nCov-2019_45_RIGHT alt, bad end"
        if entry["Primer_ID"] == "nCoV-2019_89_LEFT":
            assert (
                entry["start"] == 26835
            ), "failed to merge nCoV-2019_89_LEFT alt, bad start"
            assert (
                entry["end"] == 26860
            ), "failed to merge nCoV-2019_89_LEFT alt, bad end"

    # check access to primer scheme fields
    for row in primerScheme:
        assert "Primer_ID" in row, "failed to parse primer scheme for Primer_ID"
        assert "PoolName" in row, "failed to parse primer scheme for PoolName"

    # process the nCoV-2019 V2 primer scheme
    # this scheme has a single alt that has replaced the original primer so no alts should be collapsed
    primerScheme2 = vcftagprimersites.read_bed_file(
        TEST_DIR + "/../test-data/primer-schemes/nCoV-2019/V2/nCoV-2019.bed"
    )
    assert len(primerScheme2) == 196, "no alts should have been collapsed"

    # TODO: this is a starting point for the unit tests - more to come...
