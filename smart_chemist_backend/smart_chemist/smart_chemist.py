import io
from typing import Iterable

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
import random

from .models import AnnotatedPattern, PatternMatchingJob, PatternMatchingOutputModel


def random_string_generator(str_size, allowed_chars):
    return "".join(random.choice(allowed_chars) for x in range(str_size))


def check_subset_of_lists(list1, list2):
    if len(list1) >= len(list2):
        return False
    return check_list_in_list(list1, list2)


def check_lists_equal(list1, list2):
    if len(list1) != len(list2):
        return False
    return check_list_in_list(list1, list2)


def check_list_in_list(list1, list2):
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
    def _check_hierarchy(matches: list):
        for match in matches:
            hierarchy = match["trivial_name"]["hierarchy"]
            atoms_match = match["atom_indices"]
            if hierarchy != "[]":
                hierarchy_pattern_indexes = [int(x) + 1 for x in hierarchy.replace("[", "").replace("]", "").split(",")]
                for submatch in matches:
                    atoms_submatch = submatch["atom_indices"]
                    submatch_id = submatch["trivial_name"]["index"]
                    if submatch_id == "":
                        continue
                    if submatch_id in hierarchy_pattern_indexes and atoms_submatch[0] in atoms_match:
                        submatch["trivial_name"]["group"] = "overshadowed"

    @staticmethod
    def _check_overshadowed_patterns(matches: list):
        SmartChemist._check_hierarchy(matches)
        for match in matches:
            if match["trivial_name"]["group"] == "overshadowed":
                continue
            for submatch in matches:
                atoms_match = match["atom_indices"]
                atoms_submatch = submatch["atom_indices"]
                if check_subset_of_lists(atoms_match, atoms_submatch):
                    match["trivial_name"]["group"] = "overshadowed"
                    break
                elif check_lists_equal(atoms_match, atoms_submatch):
                    bonds1 = match["trivial_name"]["bonds"]
                    bonds2 = submatch["trivial_name"]["bonds"]
                    if bonds1 > bonds2:
                        submatch["trivial_name"]["group"] = "overshadowed"
                    elif bonds2 < bonds1:
                        match["trivial_name"]["group"] = "overshadowed"

    @staticmethod
    def _match_smarts_patterns(mol: rdkit.Chem.rdchem.Mol, remove_overshadowed_patterns: bool = False):
        matches = []
        # iterate all annotated SMARTS patterns from our database
        for db_row in AnnotatedPattern.objects.all():
            pattern = Chem.MolFromSmarts(db_row.smarts)
            if mol.HasSubstructMatch(pattern, useChirality=True):
                hit_atom_indices_list = mol.GetSubstructMatches(pattern, useChirality=True)
                for hit_atom_indices in hit_atom_indices_list:
                    matches.append(
                        {
                            "atom_indices": hit_atom_indices,
                            "trivial_name": {
                                "name": db_row.trivial_name,
                                "smarts": db_row.smarts,
                                "group": db_row.group,
                                "bonds": pattern.GetNumBonds(),
                                "hierarchy": db_row.hierarchy,
                                "index": db_row.index_file,
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
            "smiles": Chem.MolToSmiles(mol),
        }

    @staticmethod
    def handle_string_input(string_input: str, nof_molecules_allowed: int):
        smiles_list = string_input.split(",")
        final_json = []
        n = 0
        number_worked = 0
        number_skipped = 0
        number_problems = 0
        for x in smiles_list:
            if n > nof_molecules_allowed:
                number_skipped += 1
                continue
            mol = Chem.MolFromSmiles(x.strip())
            if mol is None:
                number_problems += 1
                continue
            n += 1
            number_worked += 1
            final_json.append(
                SmartChemist.mol_to_annotation_json(mol)
            )
        final_json.append(
            {"number_worked": number_worked, "number_skipped": number_skipped, "number_problems": number_problems})
        return final_json

    @staticmethod
    def handle_file_input(mol_supplier: Iterable, nof_molecules_allowed: int):
        n = 0
        number_worked = 0
        number_skipped = 0
        number_problems = 0
        final_json = []
        for mol in mol_supplier:
            if n > nof_molecules_allowed:
                number_skipped += 1
                continue
            if mol is None:
                print("didnt work....")
                number_problems += 1
                continue
            n += 1
            number_worked += 1
            final_json.append(SmartChemist.mol_to_annotation_json(mol))
        final_json.append(
            {"number_worked": number_worked, "number_skipped": number_skipped, "number_problems": number_problems})
        return final_json

    @staticmethod
    def handle_pattern_matching(job: PatternMatchingJob) -> None:
        """Handles pattern matching using the PatternMatchingJob object.

        :param job: The job object.
        """
        # read input data and call according function for the input format/type
        if job.input_info.input_format == 'smiles_list':
            output_json = SmartChemist.handle_string_input(
                job.input_info.input_string,
                job.input_info.input_max_nof_molecules_allowed)
        elif job.input_info.input_format == '.smi':
            suppl = Chem.SmilesMolSupplierFromText(job.input_info.input_string, titleLine=False)
            output_json = SmartChemist.handle_file_input(
                suppl,
                job.input_info.input_max_nof_molecules_allowed)
        elif job.input_info.input_format == '.sdf':
            suppl = Chem.ForwardSDMolSupplier(io.BytesIO(str.encode(job.input_info.input_string)))
            output_json = SmartChemist.handle_file_input(
                suppl,
                job.input_info.input_max_nof_molecules_allowed)
        else:
            raise RuntimeError("Wrong filetype")

        print('output_json', output_json)

        # save results to database
        job.output_info = PatternMatchingOutputModel(parent_pattern_matching_job=job, output_json=output_json)
        job.output_info.save()
        job.save()
