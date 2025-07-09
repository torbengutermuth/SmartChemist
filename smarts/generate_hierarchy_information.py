import pandas as pd
import subprocess
from pathlib import Path
import sys
from PIL import Image


smartscompare = Path("/home/torben/arbeit/SMARTScompareViewer_1.2.0/SMARTScompare")

citation_string = "# Usage of the SMARTS in these files is allowed only adhering to the CC-BY-ND 4.0 license in this directory from the University of Hamburg, ZBH - Center for Bioinformatics, Albert-Einstein-Ring 8-10, 22761 Germany Website https://uhh.de/amd.\n"


smartsview_path = "/local/gutermuth/naomi/bin/SmartsViewer_release"

def create_picture(queryString, imagePath):
    # start reaction viewer
    executable_reactionViewer = smartsview_path
    args_reactionViewer = ["-s ", queryString, "-o ", imagePath, "-p ","0 ", "0 ", "0 ", "0 ", "1", "0 ", "1 ", "0 "]
    subprocess.run([executable_reactionViewer] + args_reactionViewer)
    # Load the image
    image = Image.open(imagePath)

    # Convert to grayscale and get pixel data
    gray_image = image.convert('L')
    width, height = gray_image.size

    # Find the top boundary
    top = 0
    for y in range(height):
        if not all(gray_image.getpixel((x, y)) == 255 for x in range(width)):
            top = y
            break

    # Find the bottom boundary
    bottom = height
    for y in range(height - 1, -1, -1):
        if not all(gray_image.getpixel((x, y)) == 255 for x in range(width)):
            bottom = y
            break

            # Crop the image
    cropped_image = image.crop((0, top - 20, width, bottom - 20))

    return cropped_image

def uppercase_trivialnames(pandas_dataframe):
    for index, row in pandas_dataframe.iterrows():
        old_trivialname = row["trivialname"]
        new_trivialname = old_trivialname[0].upper() + old_trivialname[1:]
        print(old_trivialname, new_trivialname)
        row["trivialname"] = new_trivialname

biologicals = pd.read_csv("biologicals.csv", skiprows=1)
biologicals.sort_values(by=["trivialname", "SMARTS"], inplace=True)
with open("biologicals.csv", "w") as f:
    f.write(citation_string)
    biologicals.to_csv(f, index=None)

cyclic = pd.read_csv("cyclic.csv", skiprows=1)
cyclic.sort_values(by=["trivialname", "SMARTS"], inplace=True)
with open("cyclic.csv", "w") as f:
    f.write(citation_string)
    cyclic.to_csv(f, index=None)

functional = pd.read_csv("functional_groups.csv", skiprows=1)
functional.sort_values(by=["trivialname", "SMARTS"], inplace=True)
with open("functional_groups.csv", "w") as f:
    f.write(citation_string)
    functional.to_csv(f, index=None)

hierarchy_df = pd.concat([functional, biologicals])
hierarchy_df.reset_index(inplace=True,drop=True)
full_df = pd.concat([functional, cyclic, biologicals])
full_df.reset_index(inplace=True,drop=True)
#smarts_data = pd.read_csv("smarts_smart_chemist.csv")
#print(full_df.head())
new_data = hierarchy_df["SMARTS"]
#print(new_data.head())
new_data.to_csv("hierarchy_smarts_automatic", index=None, header=None, sep="\t")
todo = [smartscompare.as_posix(), "-m", "subsetoffirst", "-f", "hierarchy_smarts_automatic", "-M", "-1"]
process = subprocess.run(todo, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
#print(process.stdout.decode())


all_patterns = set()
subset_data = []
pattern_to_index = {}

for line in process.stdout.decode().split("\n"):
    line = line.strip()
    if not line.startswith("["):
        continue
    first_pattern = line.split(" ")[1].split("\t")[0][1:-1]
    second_pattern = line.split(" ")[3][1:-1]
    #print(first_pattern, second_pattern)
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
#print(subset_data)
for subset in subset_data:
    first_index = pattern_to_index[subset[0]]
    second_index = pattern_to_index[subset[1]]
    print(subset,first_index, second_index)
    subset_indexes.append([first_index, second_index])

j = 0
for i in range(full_df.shape[0]):
    current_data = []
    for index_pair in subset_indexes:
        if i == index_pair[1]:
            current_data.append(index_pair[0])
            print(index_pair)
            j += 1
    subset_column.append(current_data)

full_df["Hierarchy"] = subset_column

with open("smarts_with_hierarchy.csv", "w") as f:
    f.write(citation_string)
    full_df.to_csv(f)

sys.exit()

functional = pd.read_csv("functional_groups.csv", skiprows=1)
for index,row in functional.iterrows():
    name = row["trivialname"]
    smarts = row["SMARTS"]
    output_path = Path("anki_pics") / (name + ".png")
    #print(name,smarts)
    #create_picture(smarts, output_path)

biologicals = pd.read_csv("biologicals.csv", skiprows=1)
for index,row in biologicals.iterrows():
    name = row["trivialname"]
    smarts = row["SMARTS"]
    output_path = Path("anki_pics") / (name + ".png")
    #print(name,smarts)
    #create_picture(smarts, output_path)

cyclic = pd.read_csv("cyclic.csv", skiprows=1)

existing_names = []
i = 0
for index,row in cyclic.iterrows():
    name = row["trivialname"]
    smarts = row["SMARTS"]
    if name in existing_names:
        continue
    if name.__contains__(",") or name.__contains__("[") or smarts.__contains__("+") or smarts.__contains__("-"):
        continue
    existing_names.append(name)
    output_path = Path("anki_pics") / (name + ".png")
    print(name,smarts)
    i += 1
    create_picture(smarts, output_path)
print(i)
