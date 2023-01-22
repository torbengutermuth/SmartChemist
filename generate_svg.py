from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D

mol = Chem.MolFromSmiles('CN1C=NC2=C1C(=O)N(C(=O)N2C)C')
AllChem.Compute2DCoords(mol)
print(mol.GetConformer(0))
for atom in mol.GetAtoms():
    print(atom.GetIdx(), atom.GetSymbol())


print(Draw.MolToFile(mol,"caffeine.svg", imageType="svg"))