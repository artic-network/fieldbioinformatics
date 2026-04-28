"""Unit tests for artic/vcf_filter.py"""

from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest

from artic.vcf_filter import in_frame, Clair3Filter


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _variant(ref, alt):
    """Minimal variant-like object for in_frame tests (pysam-style attributes)."""
    return SimpleNamespace(ref=ref, alts=alt)


def _mock_sample(af, dp, ad=None, af_raises=False, dp_raises=False):
    """Build a mock pysam sample object supporting dict-style access."""
    sample = MagicMock()

    def _getitem(key):
        if key == "AF":
            if af_raises:
                raise KeyError("AF")
            return af
        if key == "DP":
            if dp_raises:
                raise KeyError("DP")
            return dp
        if key == "AD":
            if ad is None:
                raise KeyError("AD")
            return ad
        raise KeyError(key)

    sample.__getitem__ = MagicMock(side_effect=_getitem)
    return sample


def _mock_variant(ref="A", alt=None, qual=50.0, af=0.9, dp=100, ad=None,
                  af_raises=False, dp_raises=False):
    """Full mock variant for Clair3Filter.check_filter tests (pysam-style)."""
    if alt is None:
        alt = ("T",)
    v = MagicMock()
    v.ref = ref
    v.alts = alt
    v.qual = qual
    sample = _mock_sample(af, dp, ad=ad, af_raises=af_raises, dp_raises=dp_raises)
    v.samples = {"sample1": sample}
    return v


# ---------------------------------------------------------------------------
# in_frame
# ---------------------------------------------------------------------------

class TestInFrame:
    def test_snp_is_in_frame(self):
        assert in_frame(_variant("A", ["T"])) is True

    def test_mnp_same_length_is_in_frame(self):
        assert in_frame(_variant("ACG", ["TCA"])) is True

    def test_inframe_deletion_no_alt(self):
        # REF length divisible by 3, no ALT alleles
        assert in_frame(_variant("AAA", [])) is True

    def test_outofframe_deletion_no_alt(self):
        # REF length 2 — not divisible by 3
        assert in_frame(_variant("AA", [])) is False

    def test_inframe_insertion(self):
        # ALT is 3 bases longer than REF → in-frame
        assert in_frame(_variant("A", ["AAAA"])) is True

    def test_frameshift_insertion_plus1(self):
        assert in_frame(_variant("A", ["AA"])) is False

    def test_frameshift_insertion_plus2(self):
        assert in_frame(_variant("A", ["AAA"])) is False

    def test_inframe_deletion_alt_shorter(self):
        # REF 4 bases, ALT 1 base → diff = -3 → in-frame
        assert in_frame(_variant("AAAA", ["A"])) is True

    def test_frameshift_deletion_alt_shorter(self):
        # REF 3 bases, ALT 2 bases → diff = -1 → frameshift
        assert in_frame(_variant("AAA", ["AA"])) is False

    def test_multi_allele_raises(self):
        v = _variant("A", ["T", "C"])
        with pytest.raises(SystemExit):
            in_frame(v)


# ---------------------------------------------------------------------------
# Clair3Filter.check_filter
# ---------------------------------------------------------------------------

class TestClair3Filter:
    def _filter(self, no_frameshifts=False, min_depth=20, min_minor_allele_count=5):
        return Clair3Filter(
            no_frameshifts=no_frameshifts,
            min_depth=min_depth,
            min_minor_allele_count=min_minor_allele_count,
        )

    def test_passes_all_filters(self):
        f = self._filter()
        v = _mock_variant(ref="A", alt=["T"], qual=50.0, af=0.9, dp=100, ad=(10, 90))
        assert f.check_filter(v) == "pass"

    def test_fails_low_qual(self):
        f = self._filter()
        v = _mock_variant(qual=5.0)
        assert f.check_filter(v) == "mask"

    def test_passes_qual_at_threshold(self):
        f = self._filter()
        v = _mock_variant(qual=10.0, af=0.9, dp=100)
        assert f.check_filter(v) == "pass"

    def test_fails_below_qual_threshold(self):
        f = self._filter()
        v = _mock_variant(qual=9.9, af=0.9, dp=100)
        assert f.check_filter(v) == "mask"

    def test_fails_low_allele_frequency(self):
        f = self._filter()
        v = _mock_variant(qual=50.0, af=0.3, dp=100)
        assert f.check_filter(v) == "mask"

    def test_passes_af_at_threshold(self):
        f = self._filter()
        v = _mock_variant(qual=50.0, af=0.6, dp=100)
        assert f.check_filter(v) == "pass"

    def test_fails_af_just_below_threshold(self):
        f = self._filter()
        v = _mock_variant(qual=50.0, af=0.599, dp=100)
        assert f.check_filter(v) == "mask"

    def test_fails_frameshift_low_qual(self):
        # Frameshift indel with QUAL below min_frameshift_quality (50)
        f = self._filter()
        v = _mock_variant(ref="A", alt=["AA"], qual=30.0, af=0.9, dp=100)
        assert f.check_filter(v) == "mask"

    def test_passes_frameshift_high_qual(self):
        f = self._filter()
        v = _mock_variant(ref="A", alt=["AA"], qual=60.0, af=0.9, dp=100)
        assert f.check_filter(v) == "pass"

    def test_no_frameshifts_flag_rejects_frameshift(self):
        f = self._filter(no_frameshifts=True)
        v = _mock_variant(ref="A", alt=["AA"], qual=60.0, af=0.9, dp=100)
        assert f.check_filter(v) == "mask"

    def test_no_frameshifts_flag_allows_snp(self):
        f = self._filter(no_frameshifts=True)
        v = _mock_variant(ref="A", alt=["T"], qual=50.0, af=0.9, dp=100)
        assert f.check_filter(v) == "pass"

    def test_fails_low_depth(self):
        f = self._filter(min_depth=20)
        v = _mock_variant(qual=50.0, af=0.9, dp=5)
        assert f.check_filter(v) == "mask"

    def test_passes_depth_at_threshold(self):
        f = self._filter(min_depth=20)
        v = _mock_variant(qual=50.0, af=0.9, dp=20)
        assert f.check_filter(v) == "pass"

    def test_missing_af_raises_systemexit(self):
        f = self._filter()
        v = _mock_variant(qual=50.0, af_raises=True)
        with pytest.raises(SystemExit) as exc:
            f.check_filter(v)
        assert exc.value.code == 1

    def test_missing_dp_still_passes(self):
        """DP field absent → depth check is skipped, variant can still pass."""
        f = self._filter(min_depth=20)
        v = _mock_variant(qual=50.0, af=0.9, dp_raises=True)
        assert f.check_filter(v) == "pass"

    def test_defaults(self):
        """Verify the filter defaults are as expected."""
        f = self._filter()
        assert f.min_variant_quality == 10
        assert f.min_frameshift_quality == 50
        assert f.min_allele_frequency == 0.6
        assert f.min_minor_allele_count == 5

    def test_fails_low_minor_allele_count(self):
        """High quality, high AF, but too few alt reads → discard (not mask)."""
        f = self._filter(min_minor_allele_count=5)
        v = _mock_variant(qual=50.0, af=0.9, dp=100, ad=(97, 3))
        assert f.check_filter(v) == "discard"

    def test_passes_minor_allele_count_at_threshold(self):
        f = self._filter(min_minor_allele_count=5)
        v = _mock_variant(qual=50.0, af=0.9, dp=100, ad=(95, 5))
        assert f.check_filter(v) == "pass"

    def test_missing_ad_skips_mac_check(self):
        """Absent AD field → minor allele count check is skipped."""
        f = self._filter(min_minor_allele_count=5)
        v = _mock_variant(qual=50.0, af=0.9, dp=100, ad=None)
        assert f.check_filter(v) == "pass"

    def test_low_qual_masks_not_discards(self):
        """Low qual + low depth → mask (not discard), qual check runs first."""
        f = self._filter(min_depth=20)
        v = _mock_variant(qual=1.0, af=0.9, dp=5)
        assert f.check_filter(v) == "mask"

    def test_very_low_af_discards(self):
        """AF below min_mask_allele_frequency → discard (too low to be ambiguous)."""
        f = Clair3Filter(no_frameshifts=False, min_depth=20, min_mask_allele_frequency=0.1)
        v = _mock_variant(qual=50.0, af=0.05, dp=100)
        assert f.check_filter(v) == "discard"

    def test_af_between_thresholds_masks(self):
        """AF between min_mask_allele_frequency and min_allele_frequency → mask."""
        f = Clair3Filter(no_frameshifts=False, min_depth=20, min_mask_allele_frequency=0.1)
        v = _mock_variant(qual=50.0, af=0.3, dp=100)
        assert f.check_filter(v) == "mask"

    def test_af_at_mask_threshold_masks(self):
        """AF exactly at min_mask_allele_frequency → mask (not discard)."""
        f = Clair3Filter(no_frameshifts=False, min_depth=20, min_mask_allele_frequency=0.1)
        v = _mock_variant(qual=50.0, af=0.1, dp=100)
        assert f.check_filter(v) == "mask"

    def test_refcall_is_discarded(self):
        """Clair3 RefCall entries must be discarded, not passed or masked.

        AF for RefCalls is the REF allele frequency — it looks high but is the wrong
        field. The position is fine; keep the reference base.
        """
        f = self._filter()
        v = _mock_variant(ref="T", alt=(".",), qual=22.0, af=0.88, dp=200, ad=(177,))
        v.filter = {"RefCall"}
        assert f.check_filter(v) == "discard"

    def test_non_refcall_not_discarded_by_filter_field(self):
        """A passing variant without RefCall in FILTER is not caught by the RefCall guard."""
        f = self._filter()
        v = _mock_variant(ref="A", alt=("T",), qual=50.0, af=0.9, dp=100, ad=(10, 90))
        v.filter = {"PASS"}
        assert f.check_filter(v) == "pass"
