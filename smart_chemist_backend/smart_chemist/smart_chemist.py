import pathlib
import string
from io import StringIO

import os
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import rdMolDraw2D
import random

from .models import AnnotatedPattern
from .mol_utils import read_sdf

def random_string_generator(str_size, allowed_chars):
    return ''.join(random.choice(allowed_chars) for x in range(str_size))

class SmartChemist:

    @staticmethod
    def _check_overshadowed_patterns(matches:list):
        for match in matches:
            for submatch in matches:
                if match["trivial_name"]["name"] == submatch["trivial_name"]["name"]:
                    continue
                atoms_match = match["atom_indices"]
                atoms_submatch = submatch["atom_indices"]
                n = len(atoms_match)
                if any(atoms_match == atoms_submatch[i:i + n] for i in range(len(atoms_submatch)-n + 1)):
                    match["trivial_name"]["group"] = "overshadowed"
                    break


    @staticmethod
    def _match_smarts_patterns(mol: rdkit.Chem.rdchem.Mol):
        matches = []
        # iterate all annotated SMARTS patterns from our database
        for db_row in AnnotatedPattern.objects.all():
            pattern = Chem.MolFromSmarts(db_row.smarts)
            if mol.HasSubstructMatch(pattern):
                hit_atom_indices_list = mol.GetSubstructMatches(pattern)
                for hit_atom_indices in hit_atom_indices_list:
                    matches.append({'atom_indices': hit_atom_indices,
                                    'trivial_name': {'name': db_row.trivial_name,
                                                     'smarts': db_row.smarts,
                                                     'group': db_row.group}
                                    })
        SmartChemist._check_overshadowed_patterns(matches)
        return matches

    @staticmethod
    def _mol_to_image_str(mol: rdkit.Chem.rdchem.Mol, width: int, height: int) -> str:
        AllChem.Compute2DCoords(mol)
        if not mol.GetConformer(0):
            raise ValueError('Failed to generate 2D conformation for mol')
        drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        result = drawer.GetDrawingText()
        return result

    @staticmethod
    def mol_to_annotation_json(mol: rdkit.Chem.rdchem.Mol) -> dict:

        # search the database for annotations of substructures
        db_matches = SmartChemist._match_smarts_patterns(mol)
        if mol.HasProp("_Name"):
            name = mol.GetProp("_Name")
        else:
            name = random_string_generator(10,string.ascii_letters)
        print(name)
        return {
            'name': name,
            'svg': SmartChemist._mol_to_image_str(mol, 400, 400),
            'matches': db_matches
        }

    @staticmethod
    def handle_string_input(string_input:str):
        smiles_list = string_input.split(",")
        final_json = []
        for x in smiles_list:
            final_json.append(SmartChemist.mol_to_annotation_json(Chem.MolFromSmiles(x)))
        return final_json

    @staticmethod
    def handle_file_input(file_input:str):
        if not os.path.exists(file_input):
            raise RuntimeError("File does not exist")
        if file_input.endswith(".smi"):
            suppl = Chem.SmilesMolSupplier(file_input, titleLine=0)
        elif file_input.endswith(".sdf"):
            suppl = Chem.SDMolSupplier(file_input)
        final_json = []
        for mol in suppl:
            print("One molecule...")
            if mol is None:
                print("didnt work....")
                continue
            final_json.append(SmartChemist.mol_to_annotation_json(mol))
        return final_json
