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

    def get_indexes_of_pattern(self,pattern_name):
        list_of_indexes = self.df_all_patterns.index[self.df_all_patterns['trivialname'] == pattern_name].tolist()
        return list_of_indexes

    def assure_hierarchy_exists(self, index_more_specific, index_less_specific):
        row_of_interest = self.df_all_patterns.iloc[[index_more_specific]]
        hierarchy_list = ast.literal_eval(row_of_interest['Hierarchy'].to_list()[0])
        return int(index_less_specific) in hierarchy_list

    def test_all_patterns_are_covered(self):
        tested_pattern_names = set()
        for index,row in self.df_molecules.iterrows():
            tested_pattern_names.add(row["pattern"])
        patterns_not_tested = []
        for index,row in self.functionals.iterrows():
            tested_pattern_name = row["trivialname"]
            if tested_pattern_name not in tested_pattern_names:
                patterns_not_tested.append(tested_pattern_name)
        self.assertTrue(len(patterns_not_tested) == 0)

    def test_molecule_file(self):
        current_pattern_to_test = None
        pattern_str = None
        for index,row in FunctionGroupPatternTest.df_molecules.iterrows():
            currently_matching = False
            smiles, pattern, matching, comment, atomindices = row.values
            if pattern != current_pattern_to_test:
                current_pattern_to_test = pattern
                df = FunctionGroupPatternTest.functionals.query(f"trivialname == '{pattern}'")
                if df.shape[0] == 0:
                    df = FunctionGroupPatternTest.biologicals.query(f"trivialname == '{pattern}'")
                self.assertGreaterEqual(df.shape[0], 1)
                pattern_str = [x for x in df["SMARTS"]]
            for pattern in pattern_str:
                if not pd.isna(atomindices):
                    if check_expected_atom_matches(smiles, pattern, ast.literal_eval(atomindices)):
                        currently_matching = True
                else:
                    if check_smiles_smarts_matching(smiles, pattern) == matching:
                        currently_matching = True
            if not pd.isna(atomindices) and not currently_matching:
                print(f"Problem with atomindex matching for the following case: {smiles, pattern_str, comment}")
            elif not currently_matching:
                print(f"Problem with patternmatching for the following case: {smiles, pattern_str, pattern, atomindices, comment}")
            self.assertTrue(currently_matching)


    def test_necessary_hierarchy(self):
        index_carbonyl = self.get_index_of_pattern("Carbonyl")
        index_acyl = self.get_index_of_pattern("Acyl group")
        index_acyl_halide = self.get_index_of_pattern("Acyl halide")
        index_aldehyde = self.get_index_of_pattern("Aldehyde")
        index_carboxyl = self.get_index_of_pattern("Carboxylic acid")
        index_ketone = self.get_index_of_pattern("Ketone")
        index_imine = self.get_index_of_pattern("Imine")
        index_aldimine = self.get_index_of_pattern("Aldimine")
        index_ketimine = self.get_index_of_pattern("Ketimine")
        index_iminium = self.get_index_of_pattern("Iminium")
        index_quinoneimine1, index_quinoneimine2 = self.get_indexes_of_pattern("Quinoneimine")
        index_aldoxime = self.get_index_of_pattern("Aldoxime")
        index_ketoxime = self.get_index_of_pattern("Ketoxime")
        index_oxime = self.get_index_of_pattern("Oxime")
        index_amideoxime = self.get_index_of_pattern("Amide oxime")
        index_amidine = self.get_index_of_pattern("Amidine")
        index_carboxamidine = self.get_index_of_pattern("Carboxamidine")
        index_diacylamine = self.get_index_of_pattern("Diacylamine")
        self.assertTrue(self.assure_hierarchy_exists(index_carboxamidine, index_amidine))
        self.assertTrue(self.assure_hierarchy_exists(index_amideoxime, index_oxime))
        self.assertTrue(self.assure_hierarchy_exists(index_aldoxime, index_oxime))
        self.assertTrue(self.assure_hierarchy_exists(index_ketoxime, index_oxime))
        self.assertTrue(self.assure_hierarchy_exists(index_ketone, index_carbonyl))
        self.assertTrue(self.assure_hierarchy_exists(index_carboxyl, index_carbonyl))
        self.assertTrue(self.assure_hierarchy_exists(index_aldehyde, index_carbonyl))
        self.assertTrue(self.assure_hierarchy_exists(index_ketone, index_acyl))
        self.assertTrue(self.assure_hierarchy_exists(index_diacylamine, index_acyl))
        self.assertTrue(self.assure_hierarchy_exists(index_aldehyde, index_acyl))
        self.assertTrue(self.assure_hierarchy_exists(index_imine, index_acyl))
        self.assertTrue(self.assure_hierarchy_exists(index_iminium, index_acyl))
        self.assertTrue(self.assure_hierarchy_exists(index_acyl_halide, index_acyl))
        self.assertTrue(self.assure_hierarchy_exists(index_quinoneimine1, index_imine))
        self.assertTrue(self.assure_hierarchy_exists(index_quinoneimine2, index_imine))
        self.assertTrue(self.assure_hierarchy_exists(index_aldimine, index_imine))
        self.assertTrue(self.assure_hierarchy_exists(index_ketimine, index_imine))


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
        patterns_are_found = {x:False for x in pattern_dict}
        for index,row in self.df_all_patterns.iterrows():
            pattern_name = row["trivialname"]
            pattern_smarts = row["SMARTS"]
            patterns_are_found[pattern_name] = True
            self.assertTrue(pattern_name in pattern_dict, msg=f"Pattern {pattern_name} not found in pattern_dict.")
            self.assertTrue(pattern_smarts in pattern_dict[pattern_name], msg=f"Pattern smarts {pattern_smarts} not found in pattern_dict for pattern {pattern_name}.\n {pattern_dict[pattern_name]}")
        for pattern in patterns_are_found:
            self.assertTrue(patterns_are_found[pattern], msg=f"Pattern {pattern} not found in all smarts.")




if __name__ == "__main__":
    unittest.main()
