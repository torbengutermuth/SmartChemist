from rdkit import Chem


def read_sdf(file_handle):
    """Generator to read molecules from SDF.
    """
    with Chem.ForwardSDMolSupplier(file_handle) as fsuppl:
        for mol in fsuppl:
            if mol is None:
                continue
            yield mol
            

def read_mol2(mol2_file_handle):
    not_read = 0
    entry = []
    molecules = []
    read_entry = lambda entry: Chem.rdmolfiles.MolFromMol2Block(
        ''.join(entry),
        sanitize=True,
        removeHs=True)
    for line in mol2_file_handle:
        line = line.decode('utf-8')
        if line.startswith('@<TRIPOS>MOLECULE') and len(entry) > 0:
            mol = read_entry(entry)
            if mol:
                yield mol
            else:
                not_read = not_read + 1
            entry = []
        entry.append(line)

    if entry:
        mol = read_entry(entry)
        if mol:
            yield mol
        else:
            not_read = not_read + 1
