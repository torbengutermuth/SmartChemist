import django
import pandas as pd
import os
from pathlib import Path

HERE = Path(__file__).parent.resolve()
SMARTS_HIERARCHY_PATH = HERE.parent.joinpath("smarts", "smarts_with_hierarchy.csv").resolve()

os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'smart_chemist_backend.settings')
django.setup()
from smart_chemist.models import AnnotatedPattern

AnnotatedPattern.objects.all().delete()

data = pd.read_csv(SMARTS_HIERARCHY_PATH, skiprows=1)
for index, row in data.iterrows():
    test = AnnotatedPattern.objects.create(smarts=row['SMARTS'],
                                           trivial_name=row['trivialname'],
                                           group=row['group'],
                                           hierarchy=row["Hierarchy"],
                                           index_file=index+1)
    test.save()
