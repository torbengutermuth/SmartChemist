import django
import pandas as pd
import os, math


os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'smart_chemist_backend.settings')
django.setup()
from smart_chemist.models import AnnotatedPattern


for chunk in pd.read_csv("smarts_smart_chemist.csv", chunksize=200):
    for index, row in chunk.iterrows():
        print(row)
        test = AnnotatedPattern.objects.create(smarts=row['SMARTS'],
                                               trivial_name=row['trivialname'],
                                               group=row['group'])

        # try:
        #     test = SmartChemist.objects.create(**values)
        #     test.save()
        # except:
        #     print(f"Something wrong with row {index}")
        #     print(values)
        #     print(row)
