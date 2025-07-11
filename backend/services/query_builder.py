from typing import List, Dict


class QueryBuilder:
    """Build ChEMBL search queries for multi-property optimization papers."""
    
    def __init__(self):
        self.property_keywords = {
            'BBBP': ['"blood brain barrier"', '"BBB permeability"', '"CNS penetration"'],
            'HIA': ['"intestinal absorption"', '"oral bioavailability"', '"human intestinal absorption"'],
            'Mutag': ['"mutagenicity"', '"Ames test"', '"genotoxicity"'],
            'DRD2': ['"dopamine D2"', '"DRD2 receptor"', '"D2 antagonist"'],
            'plogP': ['"lipophilicity optimization"', '"logP optimization"', '"partition coefficient"'],
            'QED': ['"drug-likeness"', '"QED score"', '"quantitative estimate of drug-likeness"']
        }
        
        self.optimization_terms = [
            '"structure-activity relationship"',
            '"lead optimization"',
            '"ADMET optimization"',
            '"property-based optimization"',
            '"multi-parameter optimization"'
        ]
    
    def build_optimization_query(self, properties: List[str]) -> str:
        """
        Build a query for finding optimization papers based on target properties.
        
        Args:
            properties: List of property names (e.g., ['BBBP', 'HIA'])
            
        Returns:
            ChEMBL search query string
        """
        # Select up to 2 most relevant properties for focused search
        selected_props = properties[:2]
        
        # Build property part of query
        prop_queries = []
        for prop in selected_props:
            if prop in self.property_keywords:
                # Use the first (most specific) keyword for each property
                prop_queries.append(self.property_keywords[prop][0])
        
        # Fallback if no property keywords found
        if not prop_queries:
            return '"ADMET optimization" OR "property optimization" OR "multi-parameter optimization"'
        
        # Combine property keywords
        prop_part = f"({' OR '.join(prop_queries)})"
        
        # Select optimization terms for better relevance
        opt_terms = ['"lead optimization"', '"structure-activity relationship"']
        opt_part = f"({' OR '.join(opt_terms)})"
        
        # Build final query with AND to ensure both aspects are present
        query = f"{prop_part} AND {opt_part}"
        
        return query
    
    def get_optimization_indicators(self) -> List[str]:
        """
        Get list of terms that indicate a paper is about optimization.
        Used for post-filtering results.
        """
        return [
            'optimization', 'optimizing', 'optimized',
            'structure-activity', 'sar',
            'lead', 'scaffold',
            'improve', 'improved', 'enhancement',
            'design', 'discovery'
        ]