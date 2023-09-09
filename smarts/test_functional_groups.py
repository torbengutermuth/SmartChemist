import unittest
from pathlib import Path

import pandas as pd

from smarts.utils import check_expected_atom_matches


class FunctionGroupPatternTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.df_patterns = pd.read_csv(Path(__file__).parent / "functional_groups.csv", skiprows=1)

    def test_amide(self):
        df = FunctionGroupPatternTest.df_patterns.query("trivialname == 'Amide'")
        self.assertEqual(df.shape[0], 1)
        pattern_str = df.iloc[0]["SMARTS"]

        # example 3 from https://github.com/rdkit/rdkit/blob/9249ca5cc840fc72ea3bb73c2ff1d71a1fbd3f47/Contrib/IFG/ifg.py
        self.assertTrue(
            check_expected_atom_matches(
                "CC(=O)Nc1nnc(s1)S(=O)(=O)N", pattern_str, [(1, 2, 3)]
            )
        )
        # example 6 from https://github.com/rdkit/rdkit/blob/9249ca5cc840fc72ea3bb73c2ff1d71a1fbd3f47/Contrib/IFG/ifg.py
        self.assertTrue(
            check_expected_atom_matches(
                "Cc1onc(c1C(=O)NC2C3SC(C)(C)C(N3C2=O)C(=O)O)c4ccccc4",
                pattern_str,
                [(6, 7, 8)],
            )
        )
        # di-peptide: lysine, histidine
        self.assertTrue(
            check_expected_atom_matches(
                "N[C@@]([H])(CC(C)C)C(=O)N[C@@]([H])(CC1=CN=C-N1)C(=O)O",
                pattern_str,
                [(6, 7, 8)],
            )
        )
        # tri-peptide: polyalanine
        self.assertTrue(
            check_expected_atom_matches(
                "N[C@@]([H])(C)C(=O)N[C@@]([H])(C)C(=O)N[C@@]([H])(C)C(=O)O",
                pattern_str,
                [(3, 4, 5), (8, 9, 10)],
            )
        )

    def test_carbonyl(self):
        df = FunctionGroupPatternTest.df_patterns.query("trivialname == 'Carbonyl'")
        self.assertEqual(df.shape[0], 1)
        pattern_str = df.iloc[0]["SMARTS"]

        # example 15 from https://github.com/rdkit/rdkit/blob/9249ca5cc840fc72ea3bb73c2ff1d71a1fbd3f47/Contrib/IFG/ifg.py
        # match (25, 26) is part of a carboxylic acid
        self.assertTrue(
            check_expected_atom_matches(
                "CC1CN(CC(C)N1)c2c(F)c(N)c3c(=O)c(cn(C4CC4)c3c2F)C(=O)O",
                pattern_str,
                [(14, 15), (25, 26)],
            )
        )

    def test_primary_amine(self):
        df = FunctionGroupPatternTest.df_patterns.query(
            "trivialname == 'Primary Amine'"
        )
        self.assertEqual(df.shape[0], 1)
        pattern_str = df.iloc[0]["SMARTS"]

        # example 1 from https://github.com/rdkit/rdkit/blob/9249ca5cc840fc72ea3bb73c2ff1d71a1fbd3f47/Contrib/IFG/ifg.py
        self.assertTrue(
            check_expected_atom_matches(
                "Cc1nc(NS(=O)(=O)c2ccc(N)cc2)nc(C)c1",
                pattern_str,
                [(12,)],
            )
        )
        # example 2 from https://github.com/rdkit/rdkit/blob/9249ca5cc840fc72ea3bb73c2ff1d71a1fbd3f47/Contrib/IFG/ifg.py
        self.assertTrue(
            check_expected_atom_matches(
                "NC(=N)c1ccc(C=Cc2ccc(cc2O)C(=N)N)cc1",
                pattern_str,
                [(0,), (18,)],
            )
        )

    def test_secondary_amine(self):
        df = FunctionGroupPatternTest.df_patterns.query(
            "trivialname == 'Secondary Amine'"
        )
        self.assertEqual(df.shape[0], 1)
        pattern_str = df.iloc[0]["SMARTS"]

        # example 14 from https://github.com/rdkit/rdkit/blob/9249ca5cc840fc72ea3bb73c2ff1d71a1fbd3f47/Contrib/IFG/ifg.py
        self.assertTrue(
            check_expected_atom_matches(
                "CC(O)C1C2C(C)C(=C(N2C1=O)C(=O)O)SC3CNC(C3)C(=O)N(C)C",
                pattern_str,
                [(18,)],
            )
        )

    def test_secondary_amine(self):
        df = FunctionGroupPatternTest.df_patterns.query(
            "trivialname == 'Tertiary Amine'"
        )
        self.assertEqual(df.shape[0], 1)
        pattern_str = df.iloc[0]["SMARTS"]

        # example 14 from https://github.com/rdkit/rdkit/blob/9249ca5cc840fc72ea3bb73c2ff1d71a1fbd3f47/Contrib/IFG/ifg.py
        self.assertTrue(
            check_expected_atom_matches(
                "CC(O)C1C2C(C)C(=C(N2C1=O)C(=O)O)SC3CNC(C3)C(=O)N(C)C",
                pattern_str,
                [(9,), (23,)],
            )
        )


if __name__ == "__main__":
    unittest.main()
