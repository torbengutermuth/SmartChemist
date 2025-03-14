import pandas as pd
from pathlib import Path

smarts_path = Path("/work/gutermuth/smart_chemist_backend/smarts")
functional_group_path = smarts_path / "functional_groups.csv"
biological_path = smarts_path / "biologicals.csv"
functional_groups = pd.read_csv(functional_group_path, skiprows=1)
biologicals = pd.read_csv(biological_path, skiprows=1)
total = pd.concat([functional_groups, biologicals])
new_names = []
for index, row in total.iterrows():
    name = row["trivialname"]
    new_name = name.replace(" ", "_")
    new_names.append(new_name)
print(total.columns.values)
total.drop("group",axis=1, inplace=True)
total.drop("trivialname",axis=1, inplace=True)
total["trivialname"] = new_names
total.to_csv(smarts_path / "for_upload.smarts", sep=" ", index=None, header=None)
