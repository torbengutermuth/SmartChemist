import ast
import unittest
from pathlib import Path

import pandas as pd

from smarts.utils import check_expected_atom_matches, check_smiles_smarts_matching


class FunctionGroupPatternTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.df_patterns = pd.read_csv(Path(__file__).parent / "functional_groups.csv", skiprows=1)
        cls.df_molecules = pd.read_csv(Path(__file__).parent / "test_molecules_file.csv")

    def test_molecule_file(self):
        current_pattern_to_test = None
        pattern_str = None
        for index,row in FunctionGroupPatternTest.df_molecules.iterrows():
            smiles, pattern, matching, comment, atomindices = row.values
            if pattern != current_pattern_to_test:
                current_pattern_to_test = pattern
                df = FunctionGroupPatternTest.df_patterns.query(f"trivialname == '{pattern}'")
                self.assertEqual(df.shape[0], 1)
                pattern_str = df.iloc[0]["SMARTS"]
            if not pd.isna(atomindices):
                self.assertTrue(check_expected_atom_matches(smiles, pattern_str, ast.literal_eval(atomindices)), msg=f"Problem with atomindex matching for the following case: {smiles, pattern_str, comment}")
            else:
                self.assertEqual(check_smiles_smarts_matching(smiles, pattern_str), matching, msg=f"Problem with patternmatching for the following case: {smiles, pattern_str, pattern , atomindices, comment}")


if __name__ == "__main__":
    unittest.main()
