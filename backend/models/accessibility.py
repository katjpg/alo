"""
Models for accessibility validation endpoints.
"""

from pydantic import BaseModel
from typing import List, Dict, Optional, Any


class PropertyTarget(BaseModel):
    """Target property values."""
    # Physicochemical properties
    mw: Optional[float] = None
    logp: Optional[float] = None
    hba: Optional[int] = None
    hbd: Optional[int] = None
    tpsa: Optional[float] = None
    rotatable_bonds: Optional[int] = None
    aromatic_rings: Optional[int] = None
    fsp3: Optional[float] = None
    qed: Optional[float] = None
    
    # ADMET properties
    bbbp: Optional[float] = None  # Blood-brain barrier permeability (0-1)
    hia: Optional[float] = None   # Human intestinal absorption (0-1)
    mutag: Optional[float] = None # Mutagenicity (0-1)
    drd2: Optional[float] = None  # Dopamine receptor D2 (0-1)
    plogp: Optional[float] = None # Partition coefficient


class PropertyAnalysis(BaseModel):
    """Analysis of a single property."""
    current: float
    target: float
    gap: float
    achievable: bool
    difficulty: Optional[str] = None  # 'easy', 'moderate', 'hard', 'impossible'


class AccessibilityRequest(BaseModel):
    """Request model for accessibility check."""
    smiles: str
    target_properties: Optional[PropertyTarget] = None


class AccessibilityResponse(BaseModel):
    """Response model for accessibility check."""
    smiles: str
    accessible: bool
    score: str  # 'pass', 'warning', or 'fail'
    issues: List[str]
    warnings: List[str]
    details: Dict[str, Any]
    sa_score: Optional[float] = None
    
    # Property feasibility results
    current_properties: Optional[Dict[str, float]] = None
    property_analysis: Optional[Dict[str, PropertyAnalysis]] = None
    required_modifications: Optional[List[Dict[str, Any]]] = None
    property_conflicts: Optional[List[str]] = None