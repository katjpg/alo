"""
Accessibility validator for hit-to-lead optimization.
Assesses whether a molecule is synthetically accessible based on SMILES string.
"""

from typing import Dict, List, Optional, Tuple
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, Lipinski
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

# Try to import SA Score - handle if not available
try:
    from rdkit.Contrib.SA_Score import sascorer
    SA_SCORE_AVAILABLE = True
except ImportError:
    SA_SCORE_AVAILABLE = False


def is_valid_molecule(smiles: str) -> Tuple[bool, str]:
    """
    Check if SMILES is valid and suitable for hit-to-lead optimization.
    
    Args:
        smiles: SMILES string to validate
        
    Returns:
        Tuple of (is_valid, message)
    """
    try:
        # Check if SMILES can be parsed
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, "Invalid SMILES string"
        
        # Check for disconnected fragments
        if '.' in smiles:
            return False, "Contains disconnected fragments"
        
        # Check molecular size limits for hit-to-lead
        num_heavy_atoms = mol.GetNumHeavyAtoms()
        if num_heavy_atoms > 50:
            return False, f"Too large for hit-to-lead ({num_heavy_atoms} heavy atoms, max 50)"
        elif num_heavy_atoms < 10:
            return False, f"Too small for hit-to-lead ({num_heavy_atoms} heavy atoms, min 10)"
        
        # Check if molecule has at least one ring (typical for drug-like molecules)
        if rdMolDescriptors.CalcNumRings(mol) == 0:
            return False, "No ring systems present - atypical for hit-to-lead"
        
        return True, "Valid molecule"
        
    except Exception as e:
        return False, f"Failed to parse SMILES: {str(e)}"


def check_reactive_groups(smiles: str) -> List[Dict[str, any]]:
    """
    Identify problematic reactive groups that make synthesis difficult.
    
    Args:
        smiles: SMILES string to check
        
    Returns:
        List of found reactive groups with details
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return []
    
    # Define reactive patterns to avoid (SMARTS patterns)
    reactive_patterns = {
        'acyl_halide': '[C](=[O])[Cl,Br,I]',
        'alkyl_halide_primary': '[CH2][Cl,Br,I]',
        'michael_acceptor': '[C]=[C]-[C]=[O]',
        'peroxide': '[OX2][OX2]',
        'hydrazine': '[NX3][NX3]',
        'epoxide': 'C1OC1',
        'aziridine': 'C1NC1',
        'isocyanate': '[NX2]=C=O',
        'isothiocyanate': '[NX2]=C=S',
        'aldehyde': '[CX3H1](=O)',
        'terminal_alkyne': '[C]#[CH]',
        'anhydride': '[CX3](=[OX1])[OX2][CX3](=[OX1])',
        'azide': '[N-]=[N+]=[NX1]',
        'diazo': '[N+]#[N]',
        'nitroso': '[NX2]=[OX1]',
        'sulfonic_acid': '[SX4](=[OX1])(=[OX1])[OX2H1]'
    }
    
    found_groups = []
    for name, smarts in reactive_patterns.items():
        pattern = Chem.MolFromSmarts(smarts)
        if pattern and mol.HasSubstructMatch(pattern):
            matches = mol.GetSubstructMatches(pattern)
            found_groups.append({
                'group': name,
                'count': len(matches),
                'severity': 'high' if name in ['acyl_halide', 'isocyanate', 'azide', 'diazo'] else 'medium'
            })
    
    return found_groups


def check_molecular_stability(smiles: str) -> Dict[str, List[str]]:
    """
    Check for groups that may be unstable during synthesis or storage.
    
    Args:
        smiles: SMILES string to check
        
    Returns:
        Dictionary with stability issues categorized by type
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return {'chemical': [], 'structural': []}
    
    stability_issues = {
        'chemical': [],
        'structural': []
    }
    
    # Chemical stability patterns
    unstable_patterns = {
        'n_oxide': '[n+][O-]',
        'nitro': '[N+](=O)[O-]',
        'quaternary_nitrogen': '[NX4+]',
        'gem_dihalide': '[CX4]([Cl,Br,I])[Cl,Br,I]',
        'enol_ether': '[CX3]=[CX3][OX2]',
        'hemiacetal': '[CX4]([OH])[OX2]',
        'imine': '[CX3]=[NX2]',
        'hydrazone': '[CX3]=[NX2][NX3]',
        'oxime': '[CX3]=[NX2][OH]'
    }
    
    for name, smarts in unstable_patterns.items():
        pattern = Chem.MolFromSmarts(smarts)
        if pattern and mol.HasSubstructMatch(pattern):
            stability_issues['chemical'].append(name)
    
    # Structural stability checks
    ring_info = mol.GetRingInfo()
    
    # Check for strained rings
    for ring in ring_info.AtomRings():
        ring_size = len(ring)
        if ring_size == 3:
            stability_issues['structural'].append('3-membered ring (high strain)')
        elif ring_size == 4:
            stability_issues['structural'].append('4-membered ring (moderate strain)')
        elif ring_size > 8:
            stability_issues['structural'].append(f'{ring_size}-membered ring (conformational flexibility)')
    
    # Check for bridgehead atoms
    bridgehead_atoms = rdMolDescriptors.CalcNumBridgeheadAtoms(mol)
    if bridgehead_atoms > 2:
        stability_issues['structural'].append(f'Complex bridged system ({bridgehead_atoms} bridgehead atoms)')
    
    # Check for spiro centers
    spiro_atoms = rdMolDescriptors.CalcNumSpiroAtoms(mol)
    if spiro_atoms > 1:
        stability_issues['structural'].append(f'Multiple spiro centers ({spiro_atoms})')
    
    return stability_issues


def check_pains(smiles: str) -> List[str]:
    """
    Check for Pan-Assay Interference Compounds (PAINS).
    
    Args:
        smiles: SMILES string to check
        
    Returns:
        List of PAINS alerts
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return []
    
    # Initialize PAINS filter catalog
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    catalog = FilterCatalog(params)
    
    # Check for PAINS matches
    pains_hits = []
    if catalog.HasMatch(mol):
        # Get all matching entries
        matches = catalog.GetMatches(mol)
        for match in matches:
            pains_hits.append(match.GetDescription())
    
    return pains_hits


def calculate_sa_score(smiles: str) -> Optional[float]:
    """
    Calculate synthetic accessibility score.
    
    Args:
        smiles: SMILES string to score
        
    Returns:
        SA score (1=easy to 10=hard) or None if not available
    """
    if not SA_SCORE_AVAILABLE:
        return None
        
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None
    
    try:
        sa_score = sascorer.calculateScore(mol)
        return round(sa_score, 2)
    except:
        return None


def check_complexity_factors(smiles: str) -> List[str]:
    """
    Check for structural features that increase synthetic complexity.
    
    Args:
        smiles: SMILES string to check
        
    Returns:
        List of complexity factors
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return []
    
    complexity_factors = []
    
    # Check stereocenters
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    if len(chiral_centers) > 2:
        complexity_factors.append(f"Multiple chiral centers ({len(chiral_centers)})")
    elif len(chiral_centers) > 0:
        complexity_factors.append(f"{len(chiral_centers)} chiral center(s)")
    
    # Check ring complexity
    ring_info = mol.GetRingInfo()
    
    # Fused rings
    ring_systems = rdMolDescriptors.CalcNumRings(mol)
    if ring_systems > 3:
        complexity_factors.append(f"Complex ring system ({ring_systems} rings)")
    
    # Check for macrocycles
    if ring_info.AtomRings():
        largest_ring = max(len(ring) for ring in ring_info.AtomRings())
        if largest_ring > 7:
            complexity_factors.append(f"Macrocycle present (size {largest_ring})")
    
    # Check heteroatom density
    num_heteroatoms = rdMolDescriptors.CalcNumHeteroatoms(mol)
    heteroatom_ratio = num_heteroatoms / mol.GetNumHeavyAtoms()
    if heteroatom_ratio > 0.5:
        complexity_factors.append(f"High heteroatom content ({heteroatom_ratio:.1%})")
    
    return complexity_factors


def validate_accessibility(smiles: str, verbose: bool = False) -> Dict[str, any]:
    """
    Main function to check synthetic accessibility for hit-to-lead optimization.
    
    Args:
        smiles: SMILES string to validate
        verbose: Whether to print detailed results
        
    Returns:
        Dictionary with validation results
    """
    results = {
        'smiles': smiles,
        'accessible': True,
        'score': 'pass',  # pass/warning/fail
        'issues': [],
        'warnings': [],
        'details': {}
    }
    
    # Step 1: Basic validation
    is_valid, msg = is_valid_molecule(smiles)
    if not is_valid:
        results['accessible'] = False
        results['score'] = 'fail'
        results['issues'].append(msg)
        return results
    
    # Step 2: Check reactive groups
    reactive_groups = check_reactive_groups(smiles)
    if reactive_groups:
        high_severity = [g for g in reactive_groups if g['severity'] == 'high']
        if high_severity:
            results['accessible'] = False
            results['score'] = 'fail'
            results['issues'].append(f"Contains highly reactive groups: {[g['group'] for g in high_severity]}")
        else:
            results['score'] = 'warning' if results['score'] == 'pass' else results['score']
            results['warnings'].append(f"Contains reactive groups: {[g['group'] for g in reactive_groups]}")
        results['details']['reactive_groups'] = reactive_groups
    
    # Step 3: Check PAINS
    pains_hits = check_pains(smiles)
    if pains_hits:
        results['accessible'] = False
        results['score'] = 'fail'
        results['issues'].append(f"PAINS alert: {pains_hits}")
        results['details']['pains'] = pains_hits
    
    # Step 4: Check stability
    stability_issues = check_molecular_stability(smiles)
    if stability_issues['chemical']:
        results['score'] = 'warning' if results['score'] == 'pass' else results['score']
        results['warnings'].append(f"Chemical stability concerns: {stability_issues['chemical']}")
    if stability_issues['structural']:
        if any('high strain' in issue for issue in stability_issues['structural']):
            results['accessible'] = False
            results['score'] = 'fail'
            results['issues'].append(f"Structural issues: {stability_issues['structural']}")
        else:
            results['score'] = 'warning' if results['score'] == 'pass' else results['score']
            results['warnings'].append(f"Structural considerations: {stability_issues['structural']}")
    results['details']['stability'] = stability_issues
    
    # Step 5: Check complexity
    complexity_factors = check_complexity_factors(smiles)
    if len(complexity_factors) > 3:
        results['accessible'] = False
        results['score'] = 'fail'
        results['issues'].append(f"Too complex: {complexity_factors}")
    elif complexity_factors:
        results['score'] = 'warning' if results['score'] == 'pass' else results['score']
        results['warnings'].append(f"Complexity factors: {complexity_factors}")
    results['details']['complexity'] = complexity_factors
    
    # Step 6: Calculate SA score (if available)
    sa_score = calculate_sa_score(smiles)
    if sa_score is not None:
        results['sa_score'] = sa_score
        if sa_score > 6:
            results['score'] = 'warning' if results['score'] == 'pass' else results['score']
            results['warnings'].append(f"High SA score: {sa_score} (>6 indicates difficult synthesis)")
    
    # Print results if verbose
    if verbose:
        print(f"\nValidation Results for: {smiles}")
        print(f"Accessible: {results['accessible']} (Score: {results['score']})")
        
        if results['issues']:
            print("\nIssues (synthesis blockers):")
            for issue in results['issues']:
                print(f"  • {issue}")
        
        if results['warnings']:
            print("\nWarnings (proceed with caution):")
            for warning in results['warnings']:
                print(f"  • {warning}")
        
        if sa_score is not None:
            print(f"\nSA Score: {sa_score} (1=easy, 10=hard)")
    
    return results