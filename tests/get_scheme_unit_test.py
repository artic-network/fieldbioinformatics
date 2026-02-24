import unittest
import os

from artic.utils import get_scheme


class get_scheme_unit_test(unittest.TestCase):
    def test_get_scheme_no_length(self):

        bed, ref, scheme_version = get_scheme(
            scheme_name="artic-pan-dengue",
            scheme_version="v1.0.0",
            scheme_directory=f"{os.getcwd()}/primer-schemes",
        )

        self.assertEqual(
            bed, f"{os.getcwd()}/primer-schemes/artic-pan-dengue/400/v1.0.0/primer.bed"
        )

        self.assertEqual(
            ref,
            f"{os.getcwd()}/primer-schemes/artic-pan-dengue/400/v1.0.0/reference.fasta",
        )

    def test_get_scheme_unclear_length(self):

        with self.assertRaises(SystemExit) as cm:
            bed, ref, scheme_version = get_scheme(
                scheme_name="hav",
                scheme_version="v1.0.0",
                scheme_directory=f"{os.getcwd()}/primer-schemes",
            )

        self.assertEqual(cm.exception.code, 1)
