from pydantic import BaseModel
from typing import Optional, List, Literal

class LigandPropRequest(BaseModel):
    smiles: str

class LigandPropResponse(BaseModel):
    smiles: str
    properties: dict
    
# Validation models
class LigandValidateRequest(BaseModel):
    ligand_data: str
    input_type: Literal["smiles", "sdf", "mol"]

class LigandValidateResponse(BaseModel):
    is_valid: bool
    canonical_smiles: Optional[str] = None