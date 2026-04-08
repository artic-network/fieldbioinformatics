"""
Unit tests for artic.make_depth_mask.collect_depths.

Synthetic BAM files are generated programmatically with pysam into pytest's
tmp_path — no large test files are committed to the repo.
"""
import numpy as np
import pytest
import pysam

from artic.make_depth_mask import collect_depths


REF_NAME = "ref"
REF_LEN = 200


def _make_segment(header, qname, rg, ref_start, cigartuples):
    """Build a minimal mapped AlignedSegment."""
    # Only CIGAR ops that consume query bases contribute to query_sequence length
    qlen = sum(length for op, length in cigartuples if op in (0, 1, 4, 7, 8))
    seg = pysam.AlignedSegment(header)
    seg.query_name = qname
    seg.query_sequence = "A" * qlen
    seg.flag = 0
    seg.reference_id = 0
    seg.reference_start = ref_start
    seg.mapping_quality = 60
    seg.cigartuples = cigartuples
    seg.query_qualities = pysam.qualitystring_to_array("I" * qlen)
    seg.set_tag("RG", rg)
    return seg


@pytest.fixture
def make_bam(tmp_path):
    """
    Fixture factory.  Call as:
        path = make_bam(read_groups, reads, indexed=True)

    read_groups : list of RG ID strings
    reads       : list of (query_name, rg_id, ref_start, cigartuples)
    indexed     : if True, sort + index the BAM
    """
    def _factory(read_groups, reads, indexed=True):
        header = pysam.AlignmentHeader.from_dict({
            "HD": {"VN": "1.6"},
            "SQ": [{"SN": REF_NAME, "LN": REF_LEN}],
            "RG": [{"ID": rg} for rg in read_groups],
        })
        unsorted = str(tmp_path / "unsorted.bam")
        path = str(tmp_path / "test.bam")
        with pysam.AlignmentFile(unsorted, "wb", header=header) as bam:
            for qname, rg, ref_start, cigar in reads:
                bam.write(_make_segment(header, qname, rg, ref_start, cigar))
        pysam.sort("-o", path, unsorted)
        if indexed:
            pysam.index(path)
        return path

    return _factory


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_basic_depth_single_rg(make_bam):
    """A single 50M read in one RG should give depth=1 for positions 0-49."""
    path = make_bam(
        read_groups=["pool1"],
        reads=[("r1", "pool1", 0, [(0, 50)])],  # 50M
    )
    depths, rg_depths = collect_depths(path, REF_NAME, minDepth=1, ignoreDeletions=False, warnRGcov=False)

    assert np.all(depths[0:50] == 1), "covered positions should have depth 1"
    assert np.all(depths[50:] == 0), "uncovered positions should have depth 0"
    assert np.array_equal(depths, rg_depths["pool1"]), "single RG should match combined depths"


def test_min_depth_masking(make_bam):
    """Positions with depth < minDepth should be zeroed out."""
    path = make_bam(
        read_groups=["pool1"],
        reads=[("r1", "pool1", 0, [(0, 50)])],  # depth=1 everywhere
    )
    depths, _ = collect_depths(path, REF_NAME, minDepth=2, ignoreDeletions=False, warnRGcov=False)

    assert np.all(depths == 0), "depth=1 < minDepth=2 so all positions should be masked"


def test_two_rgs_non_overlapping(make_bam):
    """Two RGs with non-overlapping reads; each RG's depth array should only cover its own region."""
    path = make_bam(
        read_groups=["pool1", "pool2"],
        reads=[
            ("r1", "pool1", 0,  [(0, 50)]),   # pool1 covers 0-49
            ("r2", "pool2", 50, [(0, 50)]),   # pool2 covers 50-99
        ],
    )
    depths, rg_depths = collect_depths(path, REF_NAME, minDepth=1, ignoreDeletions=False, warnRGcov=False)

    assert np.all(depths[0:50] == 1)
    assert np.all(depths[50:100] == 1)
    assert np.all(depths[100:] == 0)
    assert np.all(rg_depths["pool1"][0:50] == 1)
    assert np.all(rg_depths["pool1"][50:] == 0)
    assert np.all(rg_depths["pool2"][0:50] == 0)
    assert np.all(rg_depths["pool2"][50:100] == 1)


def test_low_rg_coverage_masking(make_bam):
    """
    Combined depth >= minDepth but NO individual RG reaches minDepth —
    positions should be masked to 0.
    pool1: 1 read (depth=1), pool2: 1 read (depth=1) → combined=2, minDepth=2
    Each RG has depth 1 < 2, so the position should be masked.
    """
    path = make_bam(
        read_groups=["pool1", "pool2"],
        reads=[
            ("r1", "pool1", 0, [(0, 50)]),
            ("r2", "pool2", 0, [(0, 50)]),
        ],
    )
    depths, _ = collect_depths(path, REF_NAME, minDepth=2, ignoreDeletions=False, warnRGcov=False)

    assert np.all(depths[0:50] == 0), (
        "combined depth passes minDepth but neither RG individually does — should be masked"
    )


def test_deletion_counted_by_default(make_bam):
    """With ignoreDeletions=False, D-op positions contribute to depth."""
    # cigar: 25M 5D 25M  →  ref positions 0-54, deletion at 25-29
    path = make_bam(
        read_groups=["pool1"],
        reads=[("r1", "pool1", 0, [(0, 25), (2, 5), (0, 25)])],
    )
    depths, _ = collect_depths(path, REF_NAME, minDepth=1, ignoreDeletions=False, warnRGcov=False)

    assert np.all(depths[0:55] == 1), "deletion positions should be counted when ignoreDeletions=False"
    assert np.all(depths[55:] == 0)


def test_deletion_ignored(make_bam):
    """With ignoreDeletions=True, D-op positions should have depth 0."""
    path = make_bam(
        read_groups=["pool1"],
        reads=[("r1", "pool1", 0, [(0, 25), (2, 5), (0, 25)])],
    )
    depths, _ = collect_depths(path, REF_NAME, minDepth=1, ignoreDeletions=True, warnRGcov=False)

    assert np.all(depths[0:25] == 1)
    assert np.all(depths[25:30] == 0), "deletion positions should be 0 when ignoreDeletions=True"
    assert np.all(depths[30:55] == 1)
    assert np.all(depths[55:] == 0)


def test_reference_skip_not_counted(make_bam):
    """N-op (reference skip) positions should not contribute to depth."""
    # cigar: 25M 10N 25M  →  ref positions 0-59, skip at 25-34
    path = make_bam(
        read_groups=["pool1"],
        reads=[("r1", "pool1", 0, [(0, 25), (3, 10), (0, 25)])],
    )
    depths, _ = collect_depths(path, REF_NAME, minDepth=1, ignoreDeletions=False, warnRGcov=False)

    assert np.all(depths[0:25] == 1)
    assert np.all(depths[25:35] == 0), "N-op positions should not be counted"
    assert np.all(depths[35:60] == 1)
    assert np.all(depths[60:] == 0)


def test_indexed_and_unindexed_give_same_result(make_bam, tmp_path):
    """fetch(refName) and sequential scan fallback should produce identical output."""
    reads = [
        ("r1", "pool1", 0,  [(0, 50)]),
        ("r2", "pool2", 25, [(0, 50)]),
    ]
    rgs = ["pool1", "pool2"]

    indexed_bam = make_bam(rgs, reads, indexed=True)

    # build an unindexed copy by sorting without calling pysam.index
    unindexed_bam = str(tmp_path / "unindexed.bam")
    header = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6"},
        "SQ": [{"SN": REF_NAME, "LN": REF_LEN}],
        "RG": [{"ID": rg} for rg in rgs],
    })
    unsorted = str(tmp_path / "unsorted2.bam")
    with pysam.AlignmentFile(unsorted, "wb", header=header) as bam:
        for qname, rg, ref_start, cigar in reads:
            bam.write(_make_segment(header, qname, rg, ref_start, cigar))
    pysam.sort("-o", unindexed_bam, unsorted)
    # deliberately NOT indexing

    depths_idx, rg_idx = collect_depths(indexed_bam,   REF_NAME, minDepth=1, ignoreDeletions=False, warnRGcov=False)
    depths_raw, rg_raw = collect_depths(unindexed_bam, REF_NAME, minDepth=1, ignoreDeletions=False, warnRGcov=False)

    assert np.array_equal(depths_idx, depths_raw), "indexed and unindexed should give the same combined depths"
    for rg in rgs:
        assert np.array_equal(rg_idx[rg], rg_raw[rg]), f"RG {rg} depths differ between indexed and unindexed"
