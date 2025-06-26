import unittest

import os
import pandas as pd
import glob

from stcrpy.tcr_methods.tcr_methods import fetch_mhcs, load_mhcs
from stcrpy.tcr_processing import MHC


class TestTCRMethods(unittest.TestCase):
    def test_fetch_tcr(self):
        import stcrpy
        from stcrpy import fetch_TCRs

        tcr, = fetch_TCRs("6eqa")
        self.assertIsInstance(tcr, stcrpy.tcr_processing.abTCR)

        with self.assertWarns(UserWarning):
            non_tcr = fetch_TCRs("8zt4")
        self.assertEqual(non_tcr, [])


class TestFetchMHCs(unittest.TestCase):
    def test_tcr_class_i_in_stcrdab(self):
        mhc, = fetch_mhcs('8d5q')
        self.assertIsInstance(mhc, MHC)

    def test_tcr_class_i_not_in_stcrdab(self):
        mhc, = fetch_mhcs('7rZd')
        self.assertIsInstance(mhc, MHC)


class TestLoadMHCs(unittest.TestCase):
    def test_class_i(self):
        mhc, = load_mhcs('./test_files/8d5q.cif')
        self.assertIsInstance(mhc, MHC)
