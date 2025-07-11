"""
Property prediction service for molecular accessibility assessment.
"""

from typing import Dict, Optional
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski, rdMolDescriptors
import numpy as np


def calculate_current_properties(smiles: str) -> Optional[Dict[str, float]]:
    """
    Calculate current molecular properties.
    
    Args:
        smiles: SMILES string of the molecule
        
    Returns:
        Dictionary of calculated properties or None if invalid
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None
    
    properties = {
        # Basic physicochemical properties
        'mw': round(Descriptors.MolWt(mol), 2),
        'logp': round(Crippen.MolLogP(mol), 2),
        'hba': Lipinski.NumHAcceptors(mol),
        'hbd': Lipinski.NumHDonors(mol),
        'tpsa': round(Descriptors.TPSA(mol), 2),
        'rotatable_bonds': Descriptors.NumRotatableBonds(mol),
        'aromatic_rings': Descriptors.NumAromaticRings(mol),
        'fsp3': round(rdMolDescriptors.CalcFractionCSP3(mol), 3),
        'qed': round(Descriptors.qed(mol), 3),
        
        # Calculate plogP (same as logp for now, could be refined)
        'plogp': round(Crippen.MolLogP(mol), 2)
    }
    
    # Add ADMET predictions (placeholders for now)
    admet_props = predict_admet_properties(mol)
    properties.update(admet_props)
    
    return properties


def predict_admet_properties(mol: Chem.Mol) -> Dict[str, float]:
    """
    Predict ADMET properties using simple rules or ML models.
    
    Note: These are placeholder predictions based on simple rules.
    In production, these would use trained ML models.
    
    Args:
        mol: RDKit molecule object
        
    Returns:
        Dictionary of ADMET predictions (values between 0-1)
    """
    # Calculate descriptors for rule-based predictions
    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    tpsa = Descriptors.TPSA(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    
    # Blood-Brain Barrier Permeability (BBBP)
    # Simple rule: Low MW, moderate LogP, low TPSA favor BBB penetration
    bbbp_score = 1.0
    if mw > 450:
        bbbp_score *= 0.7
    if logp < 1 or logp > 4:
        bbbp_score *= 0.8
    if tpsa > 90:
        bbbp_score *= 0.6
    if hbd > 5:
        bbbp_score *= 0.7
    
    # Human Intestinal Absorption (HIA)
    # Rule of 5 compliance generally indicates good absorption
    hia_score = 1.0
    if mw > 500:
        hia_score *= 0.8
    if logp > 5:
        hia_score *= 0.8
    if hbd > 5:
        hia_score *= 0.9
    if hba > 10:
        hia_score *= 0.9
    
    # Mutagenicity (Mutag)
    # Check for common mutagenic substructures
    mutag_score = 0.1  # Low baseline
    mutag_patterns = [
        '[N+](=O)[O-]',  # Nitro groups
        'N=N',           # Azo compounds
        'c1ccc2c(c1)ccc3c2cccc3',  # Polycyclic aromatics
    ]
    for pattern in mutag_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            mutag_score = max(mutag_score, 0.7)
            break
    
    # Dopamine Receptor D2 (DRD2) binding
    # Placeholder - would need specific pharmacophore model
    # Basic check: aromatic rings and basic nitrogen
    drd2_score = 0.3  # Baseline
    if Descriptors.NumAromaticRings(mol) >= 2:
        drd2_score += 0.2
    # Check for basic nitrogen
    basic_n_pattern = Chem.MolFromSmarts('[NX3;H2,H1,H0]')
    if mol.HasSubstructMatch(basic_n_pattern):
        drd2_score += 0.2
    
    return {
        'bbbp': round(min(max(bbbp_score, 0), 1), 3),
        'hia': round(min(max(hia_score, 0), 1), 3),
        'mutag': round(min(max(mutag_score, 0), 1), 3),
        'drd2': round(min(max(drd2_score, 0), 1), 3)
    }


def check_property_feasibility(current_properties: Dict[str, float], 
                             target_properties: Dict[str, float]) -> tuple:
    """
    Check if target properties are achievable from current molecule.
    
    Args:
        current_properties: Current molecular properties
        target_properties: Target property values
        
    Returns:
        Tuple of (issues, warnings)
    """
    issues = []
    warnings = []
    
    for prop, target_value in target_properties.items():
        if prop not in current_properties or target_value is None:
            continue
            
        current_value = current_properties[prop]
        delta = target_value - current_value
        
        # Physicochemical property checks
        if prop == 'logp':
            if delta > 3:
                issues.append(f"LogP increase of {delta:.1f} too large - would require extensive hydrophobic additions")
            elif delta > 2:
                warnings.append(f"LogP increase of {delta:.1f} challenging - requires significant modifications")
            elif delta < -3:
                issues.append(f"LogP decrease of {abs(delta):.1f} too large - would require many polar groups")
                
        elif prop == 'mw':
            if target_value > 500 and current_value < 300:
                warnings.append("Large MW increase may reduce oral bioavailability")
            if target_value < current_value:
                issues.append("Cannot reduce molecular weight without removing atoms")
                
        elif prop == 'hba':
            if delta > 5:
                issues.append(f"Adding {delta} H-bond acceptors would require many heteroatoms")
            elif delta > 3:
                warnings.append(f"Adding {delta} H-bond acceptors requires significant modifications")
                
        elif prop == 'hbd':
            if delta > 3:
                issues.append(f"Adding {delta} H-bond donors would require many OH/NH groups")
            elif delta > 2:
                warnings.append(f"Adding {delta} H-bond donors is challenging")
                
        elif prop == 'tpsa':
            if target_value > 140 and current_value < 100:
                warnings.append("High TPSA target may compromise membrane permeability")
            if delta > 60:
                issues.append(f"TPSA increase of {delta:.1f} would require extensive polar modifications")
                
        elif prop == 'rotatable_bonds':
            if delta > 5:
                issues.append(f"Adding {delta} rotatable bonds would make molecule too flexible")
                
        # ADMET property checks
        elif prop == 'bbbp':
            if target_value > 0.8 and current_value < 0.3:
                issues.append("Large increase in BBB permeability requires major structural changes")
            elif abs(delta) > 0.5:
                warnings.append(f"Significant change in BBB permeability ({delta:.2f}) is challenging")
                
        elif prop == 'hia':
            if target_value > 0.9 and current_value < 0.5:
                warnings.append("Achieving high HIA from low baseline requires careful optimization")
                
        elif prop == 'mutag':
            if target_value < 0.2 and current_value > 0.5:
                warnings.append("Reducing mutagenicity requires removing problematic substructures")
    
    return issues, warnings


def check_property_combinations(target_properties: Dict[str, float]) -> list:
    """
    Check if property combinations are realistic.
    
    Args:
        target_properties: Target property values
        
    Returns:
        List of conflicts found
    """
    conflicts = []
    
    # Filter out None values
    targets = {k: v for k, v in target_properties.items() if v is not None}
    
    # Physicochemical conflicts
    if 'logp' in targets and 'tpsa' in targets:
        if targets['logp'] > 5 and targets['tpsa'] > 120:
            conflicts.append("High LogP (>5) with high TPSA (>120) is very rare - conflicting goals")
            
    if 'mw' in targets and 'logp' in targets:
        if targets['mw'] < 300 and targets['logp'] > 4:
            conflicts.append("Low MW (<300) with high LogP (>4) difficult - limited hydrophobic space")
            
    if 'hbd' in targets and 'logp' in targets:
        if targets['hbd'] > 5 and targets['logp'] > 4:
            conflicts.append("Many H-bond donors (>5) with high LogP (>4) are contradictory")
    
    # ADMET conflicts
    if 'bbbp' in targets and 'tpsa' in targets:
        if targets['bbbp'] > 0.8 and targets['tpsa'] > 90:
            conflicts.append("High BBB permeability with high TPSA (>90) is unlikely")
            
    if 'bbbp' in targets and 'mw' in targets:
        if targets['bbbp'] > 0.8 and targets['mw'] > 500:
            conflicts.append("High BBB permeability with MW >500 is very challenging")
    
    # Extreme value checks
    if 'mw' in targets and targets['mw'] > 600:
        conflicts.append("MW > 600 generally not suitable for oral drugs")
        
    if 'rotatable_bonds' in targets and targets['rotatable_bonds'] > 10:
        conflicts.append("Too many rotatable bonds (>10) - poor oral bioavailability")
        
    if 'qed' in targets and targets['qed'] > 0.9:
        conflicts.append("QED > 0.9 is extremely rare and difficult to achieve")
    
    return conflicts


def identify_required_modifications(current_properties: Dict[str, float],
                                  target_properties: Dict[str, float]) -> list:
    """
    Identify what modifications would be needed.
    
    Args:
        current_properties: Current molecular properties
        target_properties: Target property values
        
    Returns:
        List of required modifications
    """
    modifications = []
    
    for prop, target_value in target_properties.items():
        if prop not in current_properties or target_value is None:
            continue
            
        current_value = current_properties[prop]
        delta = target_value - current_value
        
        # LogP modifications
        if prop == 'logp' and abs(delta) > 0.5:
            if delta > 0:
                modifications.append({
                    'property': 'logp',
                    'type': 'increase_lipophilicity',
                    'magnitude': round(delta, 2),
                    'difficulty': 'moderate' if delta < 2 else 'hard',
                    'options': [
                        f'Add alkyl chains (each CH2 adds ~0.5 LogP)',
                        'Add aromatic rings (adds ~1.5-2 LogP)',
                        'Replace polar groups with hydrophobic ones',
                        'Remove H-bond donors/acceptors'
                    ]
                })
            else:
                modifications.append({
                    'property': 'logp',
                    'type': 'increase_polarity',
                    'magnitude': round(abs(delta), 2),
                    'difficulty': 'moderate' if abs(delta) < 2 else 'hard',
                    'options': [
                        'Add hydroxyl groups (each OH reduces ~1.5 LogP)',
                        'Add amide groups',
                        'Replace hydrophobic groups with polar ones',
                        'Add charged functional groups'
                    ]
                })
        
        # MW modifications
        elif prop == 'mw' and delta > 50:
            modifications.append({
                'property': 'mw',
                'type': 'increase_size',
                'magnitude': round(delta, 1),
                'difficulty': 'easy' if delta < 100 else 'moderate',
                'options': [
                    f'Need to add ~{int(delta/14)} heavy atoms',
                    'Consider adding functional groups or ring systems',
                    'Extend existing chains or scaffolds'
                ]
            })
            
        # H-bond modifications
        elif prop in ['hba', 'hbd'] and delta > 1:
            modifications.append({
                'property': prop,
                'type': f'increase_{prop}',
                'magnitude': int(delta),
                'difficulty': 'moderate' if delta < 3 else 'hard',
                'options': [
                    f'Add groups containing {"oxygen/nitrogen" if prop == "hba" else "OH/NH"}',
                    'Modify existing groups to increase H-bonding',
                    'Consider heterocyclic replacements'
                ]
            })
            
        # TPSA modifications
        elif prop == 'tpsa' and abs(delta) > 20:
            modifications.append({
                'property': 'tpsa',
                'type': 'modify_polarity',
                'magnitude': round(abs(delta), 1),
                'difficulty': 'moderate' if abs(delta) < 40 else 'hard',
                'options': [
                    'Add/remove polar functional groups',
                    'Modify heterocyclic content',
                    'Adjust H-bond donor/acceptor balance'
                ]
            })
    
    return modifications


def get_max_achievable_delta(property_name: str) -> float:
    """
    Define realistic modification limits for each property.
    
    Args:
        property_name: Name of the property
        
    Returns:
        Maximum achievable change for the property
    """
    limits = {
        # Physicochemical properties
        'logp': 3.0,              # Hard to change LogP by more than 3 units
        'mw': 200,                # Adding more than 200 Da gets difficult
        'hba': 5,                 # Adding >5 acceptors is challenging
        'hbd': 3,                 # Adding >3 donors is challenging  
        'tpsa': 80,               # TPSA changes >80 require major modifications
        'rotatable_bonds': 5,     # Adding >5 rotatable bonds impacts flexibility too much
        'aromatic_rings': 2,      # Adding >2 aromatic rings is difficult
        'fsp3': 0.5,              # Large Fsp3 changes alter entire structure
        'qed': 0.3,               # QED changes >0.3 are very challenging
        
        # ADMET properties (more conservative as these are complex)
        'bbbp': 0.5,              # BBB permeability hard to change dramatically
        'hia': 0.4,               # HIA somewhat modifiable
        'mutag': 0.6,             # Can reduce mutagenicity by removing groups
        'drd2': 0.5,              # Receptor binding is very specific
        'plogp': 3.0              # Same as logp
    }
    return limits.get(property_name, float('inf'))