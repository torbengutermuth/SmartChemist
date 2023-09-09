import pandas as pd
import subprocess
from pathlib import Path
import sys


smartscompare = Path("/local/gutermuth/naomi/bin/SMARTScompare_release")

citation_string = "# Usage of the SMARTS in these files is prohibited without proper citation/reference to Github and/or Publication\n"

biologicals = pd.read_csv("biologicals.csv", skiprows=1)
biologicals.sort_values(by=["trivialname"], inplace=True)
with open("biologicals.csv", "w") as f:
    f.write(citation_string)
    biologicals.to_csv(f, index=None)

cyclic = pd.read_csv("cyclic.csv", skiprows=1)
cyclic.sort_values(by=["trivialname"], inplace=True)
with open("cyclic.csv", "w") as f:
    f.write(citation_string)
    cyclic.to_csv(f, index=None)

functional = pd.read_csv("functional_groups.csv", skiprows=1)
functional.sort_values(by=["trivialname"], inplace=True)
with open("functional_groups.csv", "w") as f:
    f.write(citation_string)
    functional.to_csv(f, index=None)

all_df = [functional, cyclic, biologicals]
full_df = pd.concat(all_df)
#smarts_data = pd.read_csv("smarts_smart_chemist.csv")
print(full_df.head())
new_data = full_df["SMARTS"]
print(new_data.head())
new_data.to_csv("all_smarts_raw_automatic", index=None, header=None, sep="\t")

todo = [smartscompare.as_posix(), "-m", "subsetoffirst", "-f", "all_smarts_raw_automatic"]
process = subprocess.run(todo, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
print(process.stdout.decode())


all_patterns = set()
subset_data = []
pattern_to_index = {}

for line in process.stdout.decode().split("\n"):
    line = line.strip()
    if not line.startswith("["):
        continue
    first_pattern = line.split(" ")[1].split("\t")[0][1:-1]
    second_pattern = line.split(" ")[3][1:-1]
    print(first_pattern, second_pattern)
    all_patterns.add(first_pattern)
    all_patterns.add(second_pattern)
    subset_data.append([first_pattern, second_pattern])

for index,row in full_df.iterrows():
    current_smarts = row["SMARTS"]
    for pattern in all_patterns:
        if pattern == current_smarts:
            pattern_to_index[pattern] = index

subset_column = []
subset_indexes = []
print(subset_data)
for subset in subset_data:
    print(subset)
    first_index = pattern_to_index[subset[0]]
    second_index = pattern_to_index[subset[1]]
    subset_indexes.append([first_index, second_index])

j = 0
for i in range(full_df.shape[0]):
    current_data = []
    for index_pair in subset_indexes:
        if i == index_pair[1]:
            current_data.append(index_pair[0])
            j += 1
    subset_column.append(current_data)

print(subset_column)
print(len(subset_data), len(subset_indexes), j)

print(subset_column)
print(len([x for x in subset_column if not len(x) == 0]))
print(len(subset_data))
full_df["Hierarchy"] = subset_column

with open("smarts_with_hierarchy.csv", "w") as f:
    f.write(citation_string)
    full_df.to_csv(f, index=None)



