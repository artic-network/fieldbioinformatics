# align_trim_unit_test.py contains unit tests for alignment trimming
import os
import pysam
import pytest

from ..artic import align_trim
from ..artic import vcftagprimersites


# help pytest resolve where test data is kept
TEST_DIR = os.path.dirname(os.path.abspath(__file__))

# dummy primers (using min required fields)
p1 = {"start": 0, "end": 10, "direction": "+", "primerID": "primer1_LEFT"}
p2 = {"start": 30, "end": 40, "direction": "-", "primerID": "primer1_RIGHT"}
p3 = {"start": 10, "end": 20, "direction": "+", "primerID": "primer2_LEFT"}
p4 = {"start": 40, "end": 50, "direction": "-", "primerID": "primer2_RIGHT"}

# primer scheme to hold dummy primers
dummyPrimerScheme = [p1, p2, p3, p4]

# actual the primer scheme for nCov
primerScheme = vcftagprimersites.read_bed_file(
    TEST_DIR + "/../test-data/primer-schemes/nCoV-2019/V1/nCoV-2019.scheme.bed"
)

# nCov alignment segment (derived from a real nCov read)
seg1 = pysam.AlignedSegment()
seg1.query_name = "0be29940-97ae-440e-b02c-07748edeceec"
seg1.flag = 0
seg1.reference_id = 0
seg1.reference_start = 4294
seg1.mapping_quality = 60
seg1.cigarstring = "40S9M1D8M2D55M1D4M2I12M1I14M1D101M1D53M2D78M1I60M52S"
seg1.query_sequence = "CAGGTTAACACAAAGACACCGACAACTTTCTTCAGCACCTACAGTGCTTAAAAGTGTAAGTGCCTTTTACATTCTACCATCTATTATCTCTAATGAGAAGCAAGAAATTCTTGAACCTTCATACTTGGAATTTTGCGAGAAATGCTGCACATGCAGAAGAAACACGCAAATTAATGCCTGTCTGTGTGGAAACTAAAGCCATAGTTTCAACTATACAGCGTAAATATAAGGGTATTAAAATACAAGGGGTGTGGTTGATTATGGTGCTAGATTTTACTTTTACACCAGTAAAACAACTGGCGTCACTTATCAACACACTTAACGATCTAAATGAAACTCTTGTTACAATGCCACTTGGCTATGTAACACATGGCTTAGAATTTGGAAGAAGCTGCTCGGTATATGAGATCTCTCAAAGTGCCAGCTACAGTTTCTGTTGCGATTGCTGAAAGTTGTCGGTGTCTTTGTGTTAACCTTAGCAATACCCATG"
seg1.query_qualities = [30] * 490

# nCov alignment segment (derived from a real nCov read) - will result in a leading CIGAR deletion during softmasking
seg2 = pysam.AlignedSegment()
seg2.query_name = "15c86d34-a527-4506-9b0c-f62827d01555"
seg2.flag = 0
seg2.reference_id = 0
seg2.reference_start = 4294
seg2.mapping_quality = 60
seg2.cigarstring = (
    "41S9M1D17M1D69M5D1M1D40M2I12M1D41M1D117M1I3M1D4M2D11M2I2M1D18M1D18M1I25M56S"
)
seg2.query_sequence = "CCAGGTTAACACAAAGACACCGACAACTTTCTTCAGCACCTACAGTGCTTAAAAGTGTAAAAGTGCCTTTACATTCTACCATCTATTATCTCTAATGAGAAGCAAGAAATTCTTGGAACTGTTTCTTGGAATTTGCAGCTTGCACATGCAGAAGAAACACGCAAATTAATGCCTGTCTGTGTGTGGAAACTAAGCCATAGTTTCAACTATACAGCGTAAATATAAGGGTATTAAATACAAGAGGGTGTGGTTGATTATGGTGCTAGATTTTACTTTTACACCAGTAAAACAACTGTAGCGTCACTTATCAACACGCTTAACGATCTAAATGAAACTCTTGTTACAATGCACACTGGCTGTAACACATGAAACTAAATTTGGAAGAAGCTGTCGGTATATGAGATCTCTCCAAAGTGCCAGCTACGGTTTCTGTTAGGTGCTGAAAAGAAAGTTGTCGGTGTCTTTGTGTGAACCTTAGCAATACGTAACC"
seg2.query_qualities = [30] * 490

# expected softmasked CIGARs
seg1expectedCIGAR = "64S48M1D4M2I12M1I14M1D101M1D53M2D78M1I38M74S"
seg2expectedCIGAR = "67S69M5D1M1D40M2I12M1D41M1D117M1I3M1D4M2D11M2I2M1D18M1D18M1I3M78S"


def test_find_primer():
    """test for the find primer function"""

    # test the primer finder on the primers themselves
    for primer in dummyPrimerScheme:
        if primer["direction"] == "+":
            result = align_trim.find_primer(
                dummyPrimerScheme, primer["start"], primer["direction"]
            )
        else:
            result = align_trim.find_primer(
                dummyPrimerScheme, primer["end"], primer["direction"]
            )
        assert (
            result[2]["primerID"] == primer["primerID"]
        ), "find_primer did not produce the query primer, which should be nearest"

    # test against other ref positions
    result = align_trim.find_primer(dummyPrimerScheme, 8, "+")
    assert (
        result[2]["primerID"] == "primer2_LEFT"
    ), "find_primer returned incorrect primer"
    result = align_trim.find_primer(dummyPrimerScheme, 25, "-")
    assert (
        result[2]["primerID"] == "primer1_RIGHT"
    ), "find_primer returned incorrect primer"


def test_trim():
    """test for the trim function"""

    def testRunner(seg, expectedCIGAR):

        # get the nearest primers to the alignment segment
        p1 = align_trim.find_primer(primerScheme, seg.reference_start, "+")
        p2 = align_trim.find_primer(primerScheme, seg.reference_end, "-")

        # get the primer positions
        p1_position = p1[2]["end"]
        p2_position = p2[2]["start"]

        # this segment should need forward and reverse softmasking
        assert seg.reference_start < p1_position, (
            "missed a forward soft masking opportunity (read: %s)" % seg.query_name
        )
        assert seg.reference_end > p2_position, (
            "missed a reverse soft masking opportunity (read: %s)" % seg.query_name
        )

        # before masking, get the query_alignment_length and the CIGAR to use for testing later
        originalCigar = seg.cigarstring
        originalQueryAlnLength = seg.query_alignment_length

        # trim the forward primer
        try:
            align_trim.trim(seg, p1_position, False, False)
        except Exception as e:
            raise Exception(
                "problem soft masking left primer in {} (error: {})".format(
                    seg.query_name, e
                )
            )

        # check the CIGAR and query alignment length is updated
        assert seg.cigarstring != originalCigar, (
            "cigar was not updated with a softmask (read: %s)" % seg.query_name
        )
        assert seg.query_alignment_length != originalQueryAlnLength, (
            "query alignment was not updated after softmask (read: %s)" % seg.query_name
        )

        # trim the reverse primer
        try:
            align_trim.trim(seg, p2_position, True, False)
        except Exception as e:
            raise Exception(
                "problem soft masking right primer in {} (error: {})".format(
                    seg.query_name, e
                )
            )

        # check the CIGAR and query alignment length is updated
        assert seg.cigarstring != originalCigar, (
            "cigar was not updated with a softmask (read: %s)" % seg.query_name
        )
        assert seg.query_alignment_length != originalQueryAlnLength, (
            "query alignment was not updated after softmask (read: %s)" % seg.query_name
        )

        # check we have the right CIGAR
        assert seg.cigarstring == expectedCIGAR, (
            "cigar does not match expected cigar string (read: %s)" % seg.query_name
        )

        # check the query alignment now matches the expected primer product
        assert seg.reference_start >= p1_position, (
            "left primer not masked corrrectly (read: %s)" % seg.query_name
        )
        assert seg.reference_end <= p2_position, (
            "right primer not masked correctly (read: %s)" % seg.query_name
        )

    # run the test with the two alignment segments
    testRunner(seg1, seg1expectedCIGAR)
    testRunner(seg2, seg2expectedCIGAR)
