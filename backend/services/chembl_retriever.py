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
    
    def find_similar_compounds(self, smiles: str, similarity_threshold: int = 80) -> List[SimilarCompound]:
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
        Main method to find optimization papers using multi-strategy search.
        
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
        
        # Get compound and target context if similar compounds found
        compound_names = []
        target_names = []
        
        if similar_compounds:
            mol_ids = [c.molecule_chembl_id for c in similar_compounds]
            
            # Get compound names and targets for better search
            compound_names = self.query_builder.get_compound_names(mol_ids)
            target_names = self.query_builder.get_target_names(mol_ids[:3])
        
        # Multi-strategy document search
        all_docs = []
        seen_doc_ids = set()
        
        # Strategy 1: Target-based search
        if target_names:
            query = self.query_builder.build_target_query(
                target_names[0], list(properties.keys())[:2]
            )
            docs = self.search_optimization_papers(query, year_from, 5)
            
            for doc in docs:
                if doc.document_chembl_id not in seen_doc_ids:
                    seen_doc_ids.add(doc.document_chembl_id)
                    doc.search_strategy = 'target_specific'
                    all_docs.append(doc)
        
        # Strategy 2: Compound-specific search
        if compound_names and len(all_docs) < max_results:
            query = self.query_builder.build_compound_query(compound_names[0])
            docs = self.search_optimization_papers(query, year_from, 5)
            
            for doc in docs:
                if doc.document_chembl_id not in seen_doc_ids:
                    seen_doc_ids.add(doc.document_chembl_id)
                    doc.search_strategy = 'compound_specific'
                    all_docs.append(doc)
        
        # Strategy 3: Property-focused search
        if len(all_docs) < max_results:
            molecule_hint = self.query_builder.infer_molecule_type(target_names)
            query = self.query_builder.build_property_query(
                list(properties.keys())[:3], molecule_hint
            )
            remaining = max_results - len(all_docs)
            docs = self.search_optimization_papers(query, year_from, remaining * 2)
            
            # Filter for relevance
            filtered = self.filter_relevant_papers(docs, properties, remaining)
            
            for doc in filtered:
                if doc.document_chembl_id not in seen_doc_ids:
                    seen_doc_ids.add(doc.document_chembl_id)
                    doc.search_strategy = 'property_optimization'
                    all_docs.append(doc)
        
        # Sort by year (newest first) and limit
        all_docs.sort(key=lambda x: x.year or 0, reverse=True)
        final_docs = all_docs[:max_results]
        
        return {
            'similar_compounds': similar_compounds,
            'documents': final_docs,
            'query_used': query if 'query' in locals() else '',
            'compounds_searched': compound_names[:3],
            'targets_found': target_names,
            'timestamp': datetime.now().isoformat()
        }