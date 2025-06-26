import ast
import unittest
from pathlib import Path

import pandas as pd

from smarts.utils import check_expected_atom_matches, check_smiles_smarts_matching


class FunctionGroupPatternTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.functionals = pd.read_csv(Path(__file__).parent / "functional_groups.csv", skiprows=1)
        cls.biologicals = pd.read_csv(Path(__file__).parent / "biologicals.csv", skiprows=1)
        cls.cyclics = pd.read_csv(Path(__file__).parent / "cyclic.csv", skiprows=1)
        cls.df_molecules = pd.read_csv(Path(__file__).parent / "test_molecules_file.csv")
        cls.df_all_patterns = pd.read_csv(Path(__file__).parent / "smarts_with_hierarchy.csv", skiprows=1)

    def get_index_of_pattern(self,pattern_name):
        list_of_indexes = self.df_all_patterns.index[self.df_all_patterns['trivialname'] == pattern_name].tolist()
        self.assertEqual(len(list_of_indexes), 1)
        return list_of_indexes[0]

    def assure_hierarchy_exists(self, index_more_specific, index_less_specific):
        row_of_interest = self.df_all_patterns.iloc[[index_more_specific]]
        hierarchy_list = ast.literal_eval(row_of_interest['Hierarchy'].to_list()[0])
        return int(index_less_specific) in hierarchy_list


    def test_molecule_file(self):
        current_pattern_to_test = None
        pattern_str = None
        for index,row in FunctionGroupPatternTest.df_molecules.iterrows():
            smiles, pattern, matching, comment, atomindices = row.values
            if pattern != current_pattern_to_test:
                current_pattern_to_test = pattern
                df = FunctionGroupPatternTest.functionals.query(f"trivialname == '{pattern}'")
                self.assertEqual(df.shape[0], 1)
                pattern_str = df.iloc[0]["SMARTS"]
            if not pd.isna(atomindices):
                self.assertTrue(check_expected_atom_matches(smiles, pattern_str, ast.literal_eval(atomindices)), msg=f"Problem with atomindex matching for the following case: {smiles, pattern_str, comment}")
            else:
                self.assertEqual(check_smiles_smarts_matching(smiles, pattern_str), matching, msg=f"Problem with patternmatching for the following case: {smiles, pattern_str, pattern , atomindices, comment}")

    def test_necessary_hierarchy(self):
        index_carbonyl = self.get_index_of_pattern("Carbonyl")
        index_acyl = self.get_index_of_pattern("Acyl group")
        index_acyl_halide = self.get_index_of_pattern("Acyl halide")
        index_aldehyde = self.get_index_of_pattern("Aldehyde")
        index_carboxyl = self.get_index_of_pattern("Carboxylic acid")
        index_ketone = self.get_index_of_pattern("Ketone")
        index_imine = self.get_index_of_pattern("Imine")
        index_iminium = self.get_index_of_pattern("Iminium")
        index_quinoneimine = self.get_index_of_pattern("Quinoneimine")
        self.assertTrue(self.assure_hierarchy_exists(index_ketone, index_carbonyl))
        self.assertTrue(self.assure_hierarchy_exists(index_carboxyl, index_carbonyl))
        self.assertTrue(self.assure_hierarchy_exists(index_aldehyde, index_carbonyl))
        self.assertTrue(self.assure_hierarchy_exists(index_ketone, index_acyl))
        self.assertTrue(self.assure_hierarchy_exists(index_aldehyde, index_acyl))
        self.assertTrue(self.assure_hierarchy_exists(index_imine, index_acyl))
        self.assertTrue(self.assure_hierarchy_exists(index_iminium, index_acyl))
        self.assertTrue(self.assure_hierarchy_exists(index_acyl_halide, index_acyl))
        self.assertTrue(self.assure_hierarchy_exists(index_quinoneimine, index_imine))

    def test_no_differences_between_pattern_files(self):
        pattern_dict = {}
        single_dfs = [self.functionals, self.cyclics, self.biologicals]
        for df in single_dfs:
            for index,row in df.iterrows():
                pattern_name = row["trivialname"]
                pattern = row["SMARTS"]
                if pattern_name in pattern_dict:
                    pattern_dict[pattern_name].append(pattern)
                else:
                    pattern_dict[pattern_name] = [pattern]
        for index,row in self.df_all_patterns.iterrows():
            pattern_name = row["trivialname"]
            pattern_smarts = row["SMARTS"]
            self.assertTrue(pattern_name in pattern_dict, msg=f"Pattern {pattern_name} not found in pattern_dict.")
            self.assertTrue(pattern_smarts in pattern_dict[pattern_name], msg=f"Pattern smarts {pattern_smarts} not found in pattern_dict for pattern {pattern_name}.\n {pattern_dict[pattern_name]}")





if __name__ == "__main__":
    unittest.main()
