import unittest
import os

from artic.utils import get_scheme


class get_scheme_unit_test(unittest.TestCase):
    def test_get_scheme_no_length(self):

        bed, ref, scheme_version = get_scheme(
            scheme_name="nCoV-2019",
            scheme_version="v1.0.0",
            scheme_directory=f"{os.getcwd()}/primer-schemes",
        )

        self.assertEqual(
            bed, f"{os.getcwd()}/primer-schemes/artic-sars-cov-2/400/v1.0.0/primer.bed"
        )

        self.assertEqual(
            ref,
            f"{os.getcwd()}/primer-schemes/artic-sars-cov-2/400/v1.0.0/reference.fasta",
        )

    def test_get_scheme_unclear_length(self):

        with self.assertRaises(SystemExit) as cm:
            bed, ref, scheme_version = get_scheme(
                scheme_name="hbv",
                scheme_version="v1.0.0",
                scheme_directory=f"{os.getcwd()}/primer-schemes",
            )

        self.assertEqual(cm.exception.code, 1)
