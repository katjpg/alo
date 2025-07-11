"""
Routes for molecular accessibility validation.
"""

from fastapi import APIRouter, HTTPException
from models.accessibility import AccessibilityRequest, AccessibilityResponse, PropertyAnalysis
from services.accessibility_validator import validate_accessibility

router = APIRouter(tags=["Accessibility"])


@router.post("/accessibility/check", response_model=AccessibilityResponse)
async def check_accessibility(request: AccessibilityRequest):
    """
    Check synthetic accessibility of a molecule.
    
    Args:
        request: AccessibilityRequest containing SMILES string and optional target properties
        
    Returns:
        AccessibilityResponse with validation results
    """
    try:
        # Convert target properties to dict if provided
        target_props = None
        if request.target_properties:
            target_props = request.target_properties.dict(exclude_none=True)
        
        # Validate accessibility with property targets
        result = validate_accessibility(request.smiles, target_properties=target_props)
        
        # Return response model
        return AccessibilityResponse(**result)
        
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))