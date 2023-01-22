from io import StringIO

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import rdMolDraw2D

from .models import AnnotatedPattern
from .mol_utils import read_sdf


class SmartChemist:

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
                                                     'smarts': db_row.smarts}
                                    })
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

        return {
            # 'name': mol.GetProp("_Name"),
            'svg': SmartChemist._mol_to_image_str(mol, 400, 400),
            'matches': db_matches
        }

    @staticmethod
    def mol_str_to_annotation_json(mol_str: str, file_type: str) -> dict:
        """Compute SmartChemist annotations for given molstring

        :param mol_str: file string of the molecule.
        :param file_type: file type of the molecule, i.e. smi for smiles or sdf.
        """

        mol = None
        if file_type == 'smi':
            mol = Chem.MolFromSmiles(mol_str)
        elif file_type == 'sdf':
            with StringIO(mol_str) as file_handle:
                mol = next(read_sdf(file_handle), None)
        else:
            raise ValueError('Unknown file format')

        if not mol:
            raise ValueError('RDKit failed to read molecule string')

        return SmartChemist.mol_to_annotation_json(mol)
