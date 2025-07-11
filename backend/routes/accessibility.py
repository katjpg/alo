"""
Routes for molecular accessibility validation.
"""

from fastapi import APIRouter, HTTPException
from models.accessibility import AccessibilityRequest, AccessibilityResponse
from services.accessibility_validator import validate_accessibility

router = APIRouter(tags=["Accessibility"])


@router.post("/accessibility/check", response_model=AccessibilityResponse)
async def check_accessibility(request: AccessibilityRequest):
    """
    Check synthetic accessibility of a molecule.
    
    Args:
        request: AccessibilityRequest containing SMILES string
        
    Returns:
        AccessibilityResponse with validation results
    """
    try:
        # Validate accessibility (verbose=False by default)
        result = validate_accessibility(request.smiles)
        
        # Return response model
        return AccessibilityResponse(**result)
        
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))