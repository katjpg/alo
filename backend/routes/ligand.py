from fastapi import APIRouter, HTTPException
from models.ligand import (
    LigandPropRequest, LigandPropResponse,
    LigandValidateRequest, LigandValidateResponse
    )

from services.properties import calc_props
from utils.parse_validate import validate_ligand

router = APIRouter(tags=["Properties"])

@router.post("/properties/", response_model=LigandPropResponse)
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


# TODO: @router.post("/ligand/draw", response_model=LigandDrawResponse)
# TODO: def draw(request: LigandDrawRequest)
# TODO: parse ligand via parse_ligand -> draw molecule 