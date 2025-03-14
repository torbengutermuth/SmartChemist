import django
import pandas as pd
import os
from rdkit import Chem

os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'smart_chemist_backend.settings')
django.setup()
from smart_chemist.models import AnnotatedPattern

AnnotatedPattern.objects.all().delete()

data = pd.read_csv("./../smarts/smarts_with_hierarchy.csv", skiprows=1)
for index, row in data.iterrows():
    if row["group"] == "cyclic":
        mol = Chem.MolFromSmiles(row["SMARTS"])
    else:
        mol = Chem.MolFromSmarts(row["SMARTS"])
    n_nitrogen = n_sulfur = n_oxygen = n_carbon = n_halogens = n_phospor =  0
    for atom in mol.GetAtoms():
        number = atom.GetAtomicNum()
        if number == 6:
            n_carbon += 1
        elif number == 7:
            n_nitrogen += 1
        elif number == 8:
            n_oxygen += 1
        elif number == 15:
            n_phospor += 1
        elif number == 16:
            n_sulfur += 1
        elif number == 9 or number == 17 or number == 35 or number == 53:
            n_halogens += 1
    n_heavy_atoms = mol.GetNumHeavyAtoms()
    n_other = n_heavy_atoms - n_nitrogen - n_sulfur - n_oxygen - n_carbon - n_halogens - n_phospor
    number_rings = len(Chem.GetSymmSSSR(mol))
    if n_other > 0:
        print(row['trivialname'], n_heavy_atoms, number_rings,n_nitrogen, n_halogens, n_sulfur, n_carbon, n_phospor, row["SMARTS"])
    test = AnnotatedPattern.objects.create(smarts=row['SMARTS'],
                                           trivial_name=row['trivialname'],
                                           group=row['group'],
                                           hierarchy=row["Hierarchy"],
                                           index_file=index+1,
                                           heavy_atoms=n_heavy_atoms,
                                           num_rings=number_rings,
                                           n_nitrogens=n_nitrogen,
                                           n_sulfur=n_sulfur,
                                           n_oxygen=n_oxygen,
                                           n_carbon=n_carbon,
                                           n_phosphor=n_phospor,
                                           n_halogens=n_halogens,
                                           n_other_atom=n_other)
    test.save()
