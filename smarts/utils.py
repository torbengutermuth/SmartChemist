from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import MolFromSmiles


def check_expected_atom_matches(
    smiles: str, pattern_str: str, expected_atom_indices_list: list[tuple[int, ...]]
) -> bool:
    mol = MolFromSmiles(smiles)
    if not mol:
        return False
    pattern = Chem.MolFromSmarts(pattern_str)
    if not pattern:
        return False
    hit_atom_indices_list = mol.GetSubstructMatches(pattern, useChirality=True)
    if len(hit_atom_indices_list) != len(expected_atom_indices_list):
        return False
    for i, exp_atom_indices in enumerate(expected_atom_indices_list):
        if not hit_atom_indices_list[i] == exp_atom_indices:
            return False
    return True


def check_smiles_smarts_matching(smiles: str, pattern_str: str) -> bool:
    mol = MolFromSmiles(smiles)
    if not mol:
        return False
    pattern = Chem.MolFromSmarts(pattern_str)
    if not pattern:
        return False
    return mol.HasSubstructMatch(pattern, useChirality=True)