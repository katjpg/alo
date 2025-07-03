from rdkit import Chem
from typing import Optional, Tuple

def validate_ligand(ligand_data: str, input_type: str) -> Tuple[bool, Optional[str]]:
    """
    Validate a ligand from SMILES string or file content.
    Returns (is_valid, canonical_smiles)
    """
    mol = parse_ligand(ligand_data, input_type)
    
    if mol is None:
        return False, None
    
    # Sanitize the molecule
    try:
        Chem.SanitizeMol(mol)
        canonical_smiles = Chem.MolToSmiles(mol, canonical=True)
        return True, canonical_smiles
    except:
        return False, None

def parse_ligand(ligand_data: str, input_type: str) -> Optional[Chem.Mol]:
    """
    Parse ligand data into RDKit Mol object.
    """
    if input_type == "smiles":
        return Chem.MolFromSmiles(ligand_data)
    elif input_type in ["sdf", "mol"]:
        return Chem.MolFromMolBlock(ligand_data)
    return None