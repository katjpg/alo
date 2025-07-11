from pydantic import BaseModel, Field
from typing import Dict, List, Optional
from datetime import datetime


class ChemblDocument(BaseModel):
    document_chembl_id: str
    title: Optional[str] = None
    abstract: Optional[str] = None
    authors: Optional[str] = None
    year: Optional[int] = None
    journal: Optional[str] = None
    journal_full_title: Optional[str] = None
    volume: Optional[str] = None
    issue: Optional[str] = None
    first_page: Optional[str] = None
    last_page: Optional[str] = None
    doi: Optional[str] = None
    doi_chembl: Optional[str] = None
    pubmed_id: Optional[int] = None
    patent_id: Optional[str] = None


class SimilarCompound(BaseModel):
    molecule_chembl_id: str
    pref_name: Optional[str] = None
    similarity: float


class ChemblSearchRequest(BaseModel):
    smiles: str = Field(..., description="Input SMILES string")
    properties: Dict[str, float] = Field(..., description="Target properties and values")
    year_from: int = Field(default=2020, description="Filter papers from this year")
    max_results: int = Field(default=10, ge=1, le=100, description="Maximum results to return")


class SearchMetadata(BaseModel):
    total_documents: int
    properties_searched: List[str]
    year_range: str
    query_used: str
    timestamp: str


class ChemblSearchResponse(BaseModel):
    input_smiles: str
    target_properties: Dict[str, float]
    similar_compounds: List[SimilarCompound]
    documents: List[ChemblDocument]
    search_metadata: SearchMetadata
    error: Optional[str] = None