from rdkit import Chem
from rdkit.Chem import Descriptors

def calc_props(smiles: str) -> dict:
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        raise ValueError("Error: Invalid SMILES")

    props = {}
    for n, f in Descriptors.descList:
        try:
            props[n] = f(m)
        except Exception:
            props[n] = None
    return props