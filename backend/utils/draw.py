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

def draw_fragments(molecules: List[Chem.Mol], labels: List[str]) -> str:
    """
    Draw multiple fragments in a grid and return as base64 encoded PNG
    """
    # Generate 2D coords for each molecule
    for mol in molecules:
        AllChem.Compute2DCoords(mol)
    
    # Use individual molecule images + combine them
    images = []
    for mol in molecules:
        img = Draw.MolToImage(mol, size=(200, 200))
        images.append(img)
    
    # Create a simple grid via PIL
    from PIL import Image as PILImage
    
    # Calculate grid dimensions
    grid_width = 3 * 200  # 3 columns, 200px each
    grid_height = ((len(images) + 2) // 3) * 200  # rows needed
    
    # Create blank canvas
    grid_img = PILImage.new('RGB', (grid_width, grid_height), 'white')
    
    # Paste images into grid
    for i, img in enumerate(images):
        row = i // 3
        col = i % 3
        x = col * 200
        y = row * 200
        grid_img.paste(img, (x, y))
    
    # Conversion: base64
    buffer = io.BytesIO()
    grid_img.save(buffer, format='PNG')
    img_base64 = base64.b64encode(buffer.getvalue()).decode('utf-8')
    
    return img_base64

def generate_fragment_images(fragments: List[Tuple[str, str]]) -> List[Tuple[str, str]]:
    """
    Generate individual images for each fragment
    Returns list of (label, base64_image) tuples
    """
    fragment_images = []
    
    for smiles, label in fragments:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            img_base64 = draw_molecule(mol, size=(200, 200))
            fragment_images.append((label, img_base64))
    
    return fragment_images