from fastapi import APIRouter, HTTPException
from models.ligand import (
    LigandPropRequest, LigandPropResponse,
    LigandValidateRequest, LigandValidateResponse,
    LigandDrawRequest, LigandDrawResponse
    )

from services.properties import calc_props
from utils.parse_validate import validate_ligand, parse_ligand
from utils.draw import draw_molecule

router = APIRouter(tags=["Ligand"])

@router.post("/ligand/properties", response_model=LigandPropResponse)
def get_props(request: LigandPropRequest):
    try:
        props = calc_props(request.smiles)
        return LigandPropResponse(smiles=request.smiles, properties=props)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    
@router.post("/ligand/validate", response_model=LigandValidateResponse)
def validate(request: LigandValidateRequest):
    is_valid, canonical_smiles = validate_ligand(request.ligand_data, request.input_type)
    return LigandValidateResponse(is_valid=is_valid, canonical_smiles=canonical_smiles)


@router.post("/ligand/draw", response_model=LigandDrawResponse)
def draw(request: LigandDrawRequest):
    # Parse ligand
    mol = parse_ligand(request.ligand_data, request.input_type)
    if mol is None:
        raise HTTPException(status_code=400, detail="Invalid ligand data")
    
    # Draw single molecule
    image = draw_molecule(mol, size=request.image_size)
    
    return LigandDrawResponse(image=image)