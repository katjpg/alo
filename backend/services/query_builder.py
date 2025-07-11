import logging
from typing import Dict, List, Optional

from chembl_webresource_client.new_client import new_client


logger = logging.getLogger(__name__)


class QueryBuilder:
    """Build sophisticated ChEMBL search queries using compound and target context."""
    
    def __init__(self):
        self.molecule = new_client.molecule
        self.activity = new_client.activity
        
        # Simplified keyword mappings for faster searches
        self.property_keywords = {
            'BBBP': ['blood brain barrier', 'CNS penetration', 'BBB'],
            'HIA': ['intestinal absorption', 'oral bioavailability', 'Caco-2'],
            'Mutag': ['mutagenicity', 'Ames test', 'genotoxicity'],
            'DRD2': ['dopamine D2 receptor', 'DRD2', 'D2 receptor'],
            'plogP': ['lipophilicity', 'logP', 'partition coefficient'],
            'QED': ['drug-likeness', 'QED', 'oral drug-like']
        }
        
        self.optimization_terms = [
            '"lead optimization"',
            '"property optimization"', 
            '"ADMET optimization"',
            '"multi-parameter optimization"'
        ]
    
    def get_compound_names(self, chembl_ids: List[str]) -> List[str]:
        """Get compound names with fallback to synonyms."""
        compound_names = []
        
        try:
            molecules = self.molecule.filter(
                molecule_chembl_id__in=chembl_ids
            ).only(['molecule_chembl_id', 'pref_name', 'molecule_synonyms'])
            
            for mol in molecules:
                if mol.get('pref_name'):
                    compound_names.append(mol['pref_name'])
                elif mol.get('molecule_synonyms'):
                    for syn in mol['molecule_synonyms'][:2]:
                        if syn.get('molecule_synonym'):
                            compound_names.append(syn['molecule_synonym'])
                            break
                            
        except Exception as e:
            logger.error(f"Failed to get compound names: {e}")
            
        return compound_names[:5]
    
    def get_target_names(self, chembl_ids: List[str]) -> List[str]:
        """Get target names for compounds to improve search relevance."""
        target_names = []
        
        try:
            activities = self.activity.filter(
                molecule_chembl_id__in=chembl_ids[:3]  # Limit for speed
            ).only(['target_pref_name'])[:20]
            
            seen_targets = set()
            for act in activities:
                target = act.get('target_pref_name')
                if target and target not in seen_targets:
                    seen_targets.add(target)
                    target_names.append(target)
                    
        except Exception as e:
            logger.error(f"Failed to get target names: {e}")
            
        return target_names[:5]
    
    def infer_molecule_type(self, target_names: List[str]) -> Optional[str]:
        """Infer molecule type from target names."""
        if not target_names:
            return None
            
        target_str = ' '.join(target_names).lower()
        
        if 'kinase' in target_str:
            return 'kinase inhibitor'
        elif 'protease' in target_str:
            return 'protease inhibitor'
        elif 'receptor' in target_str:
            if 'gpcr' in target_str or 'g-protein' in target_str:
                return 'GPCR ligand'
            else:
                return 'receptor ligand'
        elif 'enzyme' in target_str:
            return 'enzyme inhibitor'
            
        return None
    
    def build_target_query(self, target_name: str, properties: List[str]) -> str:
        """Build query for target-specific papers with property optimization."""
        clean_target = target_name.replace('(', '').replace(')', '')
        
        prop_term = ""
        if properties and properties[0] in self.property_keywords:
            prop_term = f' AND {self.property_keywords[properties[0]][0]}'
        
        return f'"{clean_target}" AND ("lead optimization" OR "structure-activity"){prop_term}'
    
    def build_compound_query(self, compound_name: str) -> str:
        """Build query for specific compound."""
        return f'"{compound_name}" AND ("optimization" OR "structure-activity")'
    
    def build_property_query(self, properties: List[str], 
                           molecule_hint: Optional[str] = None) -> str:
        """Build query focused on property optimization."""
        prop_terms = []
        for prop in properties[:3]:  # Limit to top 3
            if prop in self.property_keywords:
                prop_terms.append(self.property_keywords[prop][0])
        
        if not prop_terms:
            return ' OR '.join(self.optimization_terms)
        
        prop_part = f'({" OR ".join(prop_terms)})'
        opt_part = f'({" OR ".join(self.optimization_terms[:2])})'
        
        if molecule_hint:
            return f'"{molecule_hint}" {prop_part} AND {opt_part}'
        else:
            return f'{prop_part} AND {opt_part}'
    
    def build_optimization_query(self, properties: List[str]) -> str:
        """Legacy method for compatibility - builds simple property query."""
        return self.build_property_query(properties)
    
    def get_optimization_indicators(self) -> List[str]:
        """Get terms that indicate optimization papers."""
        return [
            'optimization', 'optimizing', 'optimized',
            'structure-activity', 'sar',
            'lead', 'scaffold',
            'improve', 'improved', 'enhancement',
            'design', 'discovery'
        ]