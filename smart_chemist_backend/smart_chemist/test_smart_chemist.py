from django.test import TestCase
from rdkit import Chem

from .models import AnnotatedPattern
from .smart_chemist import SmartChemist


class SmartChemistTests(TestCase):


    def test__match_smarts_patterns(self):
        mol = Chem.MolFromSmiles('CN1C=NC2=C1C(=O)N(C(=O)N2C)C')

        anno_pattern = AnnotatedPattern(smarts='c=O', trivial_name='my_name')
        anno_pattern.save()

        matches = SmartChemist._match_smarts_patterns(mol)
        self.assertEqual(len(matches), 2)
        self.assertEqual(len(matches[0]['atom_indices']), 2)
        self.assertEqual(matches[0]['trivial_name']['name'], 'my_name')
        self.assertEqual(matches[0]['trivial_name']['smarts'], 'c=O')
        self.assertEqual(len(matches[1]['atom_indices']), 2)
        self.assertEqual(matches[1]['trivial_name']['name'], 'my_name')
        self.assertEqual(matches[1]['trivial_name']['smarts'], 'c=O')

    def test__mol_to_image_str(self):
        mol = Chem.MolFromSmiles('CN1C=NC2=C1C(=O)N(C(=O)N2C)C')
        image_str = SmartChemist._mol_to_image_str(mol, 400, 400)
        self.assertEqual(type(image_str), str)
        self.assertTrue(len(image_str) > 0)
        self.assertTrue(image_str.strip().endswith('svg>'))  # assert svg

