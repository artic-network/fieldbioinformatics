"""Unit tests for artic/vcf_filter.py"""

from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest

from artic.vcf_filter import in_frame, Clair3Filter


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _variant(ref, alt):
    """Minimal variant-like object for in_frame tests."""
    return SimpleNamespace(REF=ref, ALT=alt)


def _mock_variant(ref="A", alt=None, qual=50.0, af=0.9, dp=100,
                  af_raises=False, dp_raises=False):
    """Full mock variant for Clair3Filter.check_filter tests."""
    if alt is None:
        alt = ["T"]
    v = MagicMock()
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
    def _filter(self, no_frameshifts=False, min_depth=20):
        return Clair3Filter(no_frameshifts=no_frameshifts, min_depth=min_depth)

    def test_passes_all_filters(self):
        f = self._filter()
        v = _mock_variant(ref="A", alt=["T"], qual=50.0, af=0.9, dp=100)
        assert f.check_filter(v) is True

    def test_fails_low_qual(self):
        f = self._filter()
        v = _mock_variant(qual=5.0)
        assert f.check_filter(v) is False

    def test_passes_qual_at_threshold(self):
        f = self._filter()
        v = _mock_variant(qual=10.0, af=0.9, dp=100)
        assert f.check_filter(v) is True

    def test_fails_below_qual_threshold(self):
        f = self._filter()
        v = _mock_variant(qual=9.9, af=0.9, dp=100)
        assert f.check_filter(v) is False

    def test_fails_low_allele_frequency(self):
        f = self._filter()
        v = _mock_variant(qual=50.0, af=0.3, dp=100)
        assert f.check_filter(v) is False

    def test_passes_af_at_threshold(self):
        f = self._filter()
        v = _mock_variant(qual=50.0, af=0.6, dp=100)
        assert f.check_filter(v) is True

    def test_fails_af_just_below_threshold(self):
        f = self._filter()
        v = _mock_variant(qual=50.0, af=0.599, dp=100)
        assert f.check_filter(v) is False

    def test_fails_frameshift_low_qual(self):
        # Frameshift indel with QUAL below min_frameshift_quality (50)
        f = self._filter()
        v = _mock_variant(ref="A", alt=["AA"], qual=30.0, af=0.9, dp=100)
        assert f.check_filter(v) is False

    def test_passes_frameshift_high_qual(self):
        f = self._filter()
        v = _mock_variant(ref="A", alt=["AA"], qual=60.0, af=0.9, dp=100)
        assert f.check_filter(v) is True

    def test_no_frameshifts_flag_rejects_frameshift(self):
        f = self._filter(no_frameshifts=True)
        v = _mock_variant(ref="A", alt=["AA"], qual=60.0, af=0.9, dp=100)
        assert f.check_filter(v) is False

    def test_no_frameshifts_flag_allows_snp(self):
        f = self._filter(no_frameshifts=True)
        v = _mock_variant(ref="A", alt=["T"], qual=50.0, af=0.9, dp=100)
        assert f.check_filter(v) is True

    def test_fails_low_depth(self):
        f = self._filter(min_depth=20)
        v = _mock_variant(qual=50.0, af=0.9, dp=5)
        assert f.check_filter(v) is False

    def test_passes_depth_at_threshold(self):
        f = self._filter(min_depth=20)
        v = _mock_variant(qual=50.0, af=0.9, dp=20)
        assert f.check_filter(v) is True

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
        assert f.check_filter(v) is True

    def test_defaults(self):
        """Verify the filter defaults are as expected."""
        f = self._filter()
        assert f.min_variant_quality == 10
        assert f.min_frameshift_quality == 50
        assert f.min_allele_frequency == 0.6
