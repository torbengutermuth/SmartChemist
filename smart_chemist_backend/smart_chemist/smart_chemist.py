import io
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
    return "".join(random.choice(allowed_chars) for x in range(str_size))


def check_subset_of_lists(list1, list2):
    if len(list1) >= len(list2):
        return False
    for element in list1:
        found = False
        for element2 in list2:
            if element == element2:
                found = True
                break
        if not found:
            return False
    return True

class SmartChemist:
    @staticmethod
    def _check_overshadowed_patterns(matches: list):
        for match in matches:
            for submatch in matches:
                atoms_match = match["atom_indices"]
                atoms_submatch = submatch["atom_indices"]
                if check_subset_of_lists(atoms_match, atoms_submatch):
                    match["trivial_name"]["group"] = "overshadowed"
                    break

    @staticmethod
    def _match_smarts_patterns(mol: rdkit.Chem.rdchem.Mol, remove_overshadowed_patterns:bool=False):
        matches = []
        # iterate all annotated SMARTS patterns from our database
        for db_row in AnnotatedPattern.objects.all():
            pattern = Chem.MolFromSmarts(db_row.smarts)
            if mol.HasSubstructMatch(pattern):
                hit_atom_indices_list = mol.GetSubstructMatches(pattern)
                for hit_atom_indices in hit_atom_indices_list:
                    matches.append(
                        {
                            "atom_indices": hit_atom_indices,
                            "trivial_name": {
                                "name": db_row.trivial_name,
                                "smarts": db_row.smarts,
                                "group": db_row.group,
                            },
                        }
                    )
        SmartChemist._check_overshadowed_patterns(matches)
        if remove_overshadowed_patterns:
            matches = [x for x in matches if x["trivial_name"]["group"] != "overshadowed"]
        return matches

    @staticmethod
    def _mol_to_image_str(mol: rdkit.Chem.rdchem.Mol, width: int, height: int) -> str:
        AllChem.Compute2DCoords(mol)
        if not mol.GetConformer(0):
            raise ValueError("Failed to generate 2D conformation for mol")
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
            name = "No Name"
        return {
            "name": name,
            "svg": SmartChemist._mol_to_image_str(mol, 400, 400),
            "matches": db_matches,
        }

    @staticmethod
    def handle_string_input(string_input: str):
        smiles_list = string_input.split(",")
        final_json = []
        for x in smiles_list:
            mol = Chem.MolFromSmiles(x.strip())
            if mol is None:
                continue
            final_json.append(
                SmartChemist.mol_to_annotation_json(mol)
            )
        return final_json

    @staticmethod
    def handle_file_input(file_input: io.BytesIO, file_name: str):
        if file_name.endswith(".smi"):
            file_string = file_input.read().decode("ascii")
            suppl = Chem.SmilesMolSupplierFromText(file_string)
        elif file_name.endswith(".sdf"):
            suppl = Chem.ForwardSDMolSupplier(file_input)
        else:
            raise RuntimeError("Wrong filetype")
        final_json = []
        for mol in suppl:
            print("One molecule...")
            if mol is None:
                print("didnt work....")
                continue
            final_json.append(SmartChemist.mol_to_annotation_json(mol))
        return final_json
