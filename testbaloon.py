#from smart_chemist_backend.smart_chemist.smart_chemist import SmartChemist
from rdkit import Chem
import ast

def check_subset_of_lists(list1, list2):
    if len(list1) >= len(list2):
        return False
    return check_list_in_list(list1, list2)


def check_lists_equal(list1, list2):
    if len(list1) != len(list2):
        return False
    return check_list_in_list(list1, list2)


def check_list_in_list(list1, list2):
    for element in list1:
        found = False
        for element2 in list2:
            if element == element2:
                found = True
                break
        if not found:
            return False
    return True

def check_hierarchy(matches: list):
    for match in matches:
        hierarchy = match["trivial_name"]["hierarchy"]
        atoms_match = match["atom_indices"]
        if hierarchy != "[]":
            hierarchy_pattern_indexes = [int(x) + 1 for x in hierarchy.replace("[", "").replace("]", "").split(",")]
            for submatch in matches:
                atoms_submatch = submatch["atom_indices"]
                submatch_id = submatch["trivial_name"]["index"]
                if submatch_id is "":
                    continue
                if submatch_id in hierarchy_pattern_indexes and atoms_submatch[0] in atoms_match:
                    print(f"overshadow hierarchy {match['trivial_name']['name']} {submatch['trivial_name']['name']}")
                    print(submatch_id, hierarchy_pattern_indexes)
                    submatch["trivial_name"]["group"] = "overshadowed"


def check_overshadowed_patterns(matches: list):
    check_hierarchy(matches)
    for match in matches:
        if match["trivial_name"]["group"] == "overshadowed":
            continue
        for submatch in matches:
            atoms_match = match["atom_indices"]
            atoms_submatch = submatch["atom_indices"]
            if check_subset_of_lists(atoms_match, atoms_submatch):
                match["trivial_name"]["group"] = "overshadowed"
                break
            elif check_lists_equal(atoms_match, atoms_submatch):
                bonds1 = match["trivial_name"]["bonds"]
                bonds2 = submatch["trivial_name"]["bonds"]
                if bonds1 > bonds2:
                    submatch["trivial_name"]["group"] = "overshadowed"
                elif bonds2 < bonds1:
                    match["trivial_name"]["group"] = "overshadowed"

current_debug_output = ast.literal_eval("[{\"atom_indices\":[10],\"trivial_name\":{\"name\":\"Methyl\",\"bonds\":0,\"group\":\"functional_group\",\"index\":91,\"smarts\":\"[#6H3X4]\",\"hierarchy\":\"[]\"}},{\"atom_indices\":[13],\"trivial_name\":{\"name\":\"Methyl\",\"bonds\":0,\"group\":\"functional_group\",\"index\":91,\"smarts\":\"[#6H3X4]\",\"hierarchy\":\"[]\"}},{\"atom_indices\":[0,1,2,3],\"trivial_name\":{\"name\":\"Sulfonamide\",\"bonds\":3,\"group\":\"functional_group\",\"index\":132,\"smarts\":\"[SX4!R](=O)(=O)[NX3]\",\"hierarchy\":\"[21]\"}},{\"atom_indices\":[11,12,4,5,6],\"trivial_name\":{\"name\":\"1,3,4-Thiadiazole\",\"bonds\":5,\"group\":\"cyclic\",\"index\":2900,\"smarts\":\"C1=NN=CS1\",\"hierarchy\":\"[]\"}},{\"atom_indices\":[6,7],\"trivial_name\":{\"name\":\"Acyl group\",\"bonds\":1,\"group\":\"overshadowed\",\"index\":15,\"smarts\":\"[#6,#15,#16;X3](=[#8,#7,#15,#16])\",\"hierarchy\":\"[14, 35, 6, 23]\"}},{\"atom_indices\":[8,9],\"trivial_name\":{\"name\":\"Acyl group\",\"bonds\":1,\"group\":\"overshadowed\",\"index\":15,\"smarts\":\"[#6,#15,#16;X3](=[#8,#7,#15,#16])\",\"hierarchy\":\"[14, 35, 6, 23]\"}},{\"atom_indices\":[0,1,3],\"trivial_name\":{\"name\":\"Amide\",\"bonds\":2,\"group\":\"overshadowed\",\"index\":22,\"smarts\":\"[*!R](=O)[#7X3]\",\"hierarchy\":\"[11, 14, 35, 6, 23]\"}},{\"atom_indices\":[0,2,3],\"trivial_name\":{\"name\":\"Amide\",\"bonds\":2,\"group\":\"overshadowed\",\"index\":22,\"smarts\":\"[*!R](=O)[#7X3]\",\"hierarchy\":\"[11, 14, 35, 6, 23]\"}},{\"atom_indices\":[8,9],\"trivial_name\":{\"name\":\"Carbonyl\",\"bonds\":1,\"group\":\"overshadowed\",\"index\":36,\"smarts\":\"[#6]=O\",\"hierarchy\":\"[14, 35, 6, 23]\"}},{\"atom_indices\":[0,1,2],\"trivial_name\":{\"name\":\"Sulfone\",\"bonds\":2,\"group\":\"overshadowed\",\"index\":134,\"smarts\":\"[SX4$(*[!1])](=O)(=O)\",\"hierarchy\":\"[]\"}}]")
current_mol = "CC(=O)N=c1sc(S(N)(=O)=O)nn1C"

for output in current_debug_output:
    print(output)
    if output["trivial_name"]["group"] == "overshadowed":
        output["trivial_name"]["group"] = "not longer overshadowed"

check_hierarchy(matches=current_debug_output)