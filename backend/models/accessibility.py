"""
Models for accessibility validation endpoints.
"""

from pydantic import BaseModel
from typing import List, Dict, Optional, Any


class AccessibilityRequest(BaseModel):
    """Request model for accessibility check."""
    smiles: str


class AccessibilityResponse(BaseModel):
    """Response model for accessibility check."""
    smiles: str
    accessible: bool
    score: str  # 'pass', 'warning', or 'fail'
    issues: List[str]
    warnings: List[str]
    details: Dict[str, Any]
    sa_score: Optional[float] = None