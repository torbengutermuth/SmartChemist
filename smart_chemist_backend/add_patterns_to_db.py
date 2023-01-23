import django
import pandas as pd
import os, math


os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'smart_chemist_backend.settings')
django.setup()
from smart_chemist.models import AnnotatedPattern

data = pd.read_csv("smarts_smart_chemist.csv", index_col=0)
for index, row in data.iterrows():
    print(row)
    test = AnnotatedPattern.objects.create(smarts=row['SMARTS'],
                                           trivial_name=row['trivialname'],
                                           group=row['group'])
    test.save()