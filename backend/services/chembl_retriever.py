import logging
from datetime import datetime
from typing import Dict, List, Optional

from chembl_webresource_client.new_client import new_client

from models.chembl import ChemblDocument, SimilarCompound
from services.query_builder import QueryBuilder


logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)


class ChemblRetriever:
    """Handle all ChEMBL API interactions for document and compound retrieval."""
    
    def __init__(self):
        self.document = new_client.document
        self.molecule = new_client.molecule
        self.similarity = new_client.similarity
        self.query_builder = QueryBuilder()
    
    def find_similar_compounds(self, smiles: str, similarity_threshold: int = 85) -> List[SimilarCompound]:
        """
        Find compounds similar to the input SMILES.
        
        Args:
            smiles: Input SMILES string
            similarity_threshold: Minimum similarity percentage (default 85)
            
        Returns:
            List of similar compounds
        """
        similar_compounds = []
        
        try:
            similar = self.similarity.filter(
                smiles=smiles,
                similarity=similarity_threshold
            ).only(['molecule_chembl_id', 'pref_name', 'similarity'])[:5]
            
            for comp in similar:
                similar_compounds.append(SimilarCompound(
                    molecule_chembl_id=comp['molecule_chembl_id'],
                    pref_name=comp.get('pref_name'),
                    similarity=comp['similarity']
                ))
                
        except Exception as e:
            logger.error(f"Similarity search failed: {e}")
        
        return similar_compounds
    
    def search_optimization_papers(
        self, 
        query: str, 
        year_from: int, 
        max_results: int
    ) -> List[ChemblDocument]:
        """
        Search ChEMBL for publications matching the query.
        
        Args:
            query: Search query string
            year_from: Minimum publication year
            max_results: Maximum number of results
            
        Returns:
            List of ChEMBL documents
        """
        documents = []
        
        try:
            # Search with extra buffer to allow for filtering
            docs = self.document.search(query).filter(
                doc_type='PUBLICATION',
                year__gte=year_from
            ).only([
                'document_chembl_id',
                'title',
                'abstract',
                'authors',
                'year',
                'journal',
                'journal_full_title',
                'volume',
                'issue',
                'first_page',
                'last_page',
                'doi',
                'doi_chembl',
                'pubmed_id',
                'patent_id'
            ])[:max_results * 2]  # Get extra for filtering
            
            for doc in docs:
                documents.append(ChemblDocument(**doc))
                
        except Exception as e:
            logger.error(f"Document search failed: {e}")
        
        return documents
    
    def filter_relevant_papers(
        self, 
        documents: List[ChemblDocument], 
        properties: Dict[str, float],
        max_results: int
    ) -> List[ChemblDocument]:
        """
        Filter papers to ensure relevance to optimization goals.
        
        Args:
            documents: List of ChEMBL documents
            properties: Target properties (for context)
            max_results: Maximum number of results to return
            
        Returns:
            Filtered list of relevant documents
        """
        optimization_indicators = self.query_builder.get_optimization_indicators()
        relevant_docs = []
        
        for doc in documents:
            # Calculate relevance score based on title and abstract
            title = (doc.title or '').lower()
            abstract = (doc.abstract or '').lower()
            
            relevance_score = 0
            for indicator in optimization_indicators:
                if indicator in title:
                    relevance_score += 2  # Title matches are weighted higher
                if indicator in abstract:
                    relevance_score += 1
            
            # Include documents with sufficient relevance
            if relevance_score >= 2:
                # Remove abstract to reduce response size
                doc.abstract = None
                relevant_docs.append(doc)
                
                if len(relevant_docs) >= max_results:
                    break
        
        return relevant_docs
    
    def find_optimization_papers(
        self,
        smiles: str,
        properties: Dict[str, float],
        year_from: int = 2020,
        max_results: int = 10
    ) -> Dict:
        """
        Main method to find optimization papers based on input molecule and properties.
        
        Args:
            smiles: Input SMILES string
            properties: Target properties and values
            year_from: Minimum publication year
            max_results: Maximum results to return
            
        Returns:
            Complete search results with metadata
        """
        # Find similar compounds
        similar_compounds = self.find_similar_compounds(smiles)
        
        # Build optimization query
        query = self.query_builder.build_optimization_query(list(properties.keys()))
        
        # Search for documents
        documents = self.search_optimization_papers(query, year_from, max_results)
        
        # Filter for relevance
        filtered_documents = self.filter_relevant_papers(documents, properties, max_results)
        
        return {
            'similar_compounds': similar_compounds,
            'documents': filtered_documents,
            'query_used': query,
            'timestamp': datetime.now().isoformat()
        }