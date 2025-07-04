from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import base64
import io
from typing import List, Tuple

def draw_molecule(mol: Chem.Mol, size: Tuple[int, int] = (300, 300)) -> str:
    """
    Draw a single molecule and return as base64 encoded PNG
    """
    # Generate 2D coordinates
    AllChem.Compute2DCoords(mol)
    
    # Generate image
    img = Draw.MolToImage(mol, size=size)
    
    # Conversion: bytes -> base64
    buffer = io.BytesIO()
    img.save(buffer, format='PNG')
    img_base64 = base64.b64encode(buffer.getvalue()).decode('utf-8')
    
    return img_base64