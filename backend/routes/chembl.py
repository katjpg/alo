from datetime import datetime
from fastapi import APIRouter, HTTPException
from typing import Dict

from models.chembl import ChemblSearchRequest, ChemblSearchResponse, SearchMetadata
from services.chembl_retriever import ChemblRetriever


router = APIRouter(tags=["ChEMBL"])
retriever = ChemblRetriever()


@router.post("/chembl/search", response_model=ChemblSearchResponse)
async def search_optimization_papers(request: ChemblSearchRequest):
    """
    Search ChEMBL for papers about multi-property optimization.
    
    Finds relevant publications based on the input molecule and target properties.
    """
    try:
        # Perform the search
        results = retriever.find_optimization_papers(
            smiles=request.smiles,
            properties=request.properties,
            year_from=request.year_from,
            max_results=request.max_results
        )
        
        # Build search metadata
        search_metadata = SearchMetadata(
            total_documents=len(results['documents']),
            properties_searched=list(request.properties.keys()),
            year_range=f"{request.year_from}-present",
            query_used=results['query_used'],
            timestamp=results['timestamp']
        )
        
        # Build response
        response = ChemblSearchResponse(
            input_smiles=request.smiles,
            target_properties=request.properties,
            similar_compounds=results['similar_compounds'],
            documents=results['documents'],
            search_metadata=search_metadata
        )
        
        return response
        
    except Exception as e:
        # Log the error and return a clean error response
        error_msg = f"ChEMBL search failed: {str(e)}"
        
        # Return response with error field populated
        return ChemblSearchResponse(
            input_smiles=request.smiles,
            target_properties=request.properties,
            similar_compounds=[],
            documents=[],
            search_metadata=SearchMetadata(
                total_documents=0,
                properties_searched=list(request.properties.keys()),
                year_range=f"{request.year_from}-present",
                query_used="",
                timestamp=datetime.now().isoformat()
            ),
            error=error_msg
        )