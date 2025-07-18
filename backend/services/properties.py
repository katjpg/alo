# fix rdkit version issue
import sys
import pickle
sys.modules['rdkit.six'] = type(sys)('rdkit.six')
sys.modules['rdkit.six.moves'] = type(sys)('rdkit.six.moves')
sys.modules['rdkit.six.moves'].cPickle = pickle
sys.modules['rdkit.six'].iteritems = lambda d : d.items()

from rdkit import Chem
from rdkit.Chem import Descriptors
from utils.parse_validate import validate_ligand, parse_ligand
from gym_molecule.envs.molecule import reward_penalized_log_p
from admet_ai import ADMETModel
from tdc import Oracle

admet_model = ADMETModel()
drd2_model = Oracle(name='DRD2')

def calc_props(smiles: str) -> dict:
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        raise ValueError("Error: Invalid SMILES")

    props = {}
    for n, f in Descriptors.descList:
        try:
            props[n] = f(m)
        except Exception:
            props[n] = None
    return props

def property_score (smiles: str, properties: str) -> dict : 

    model_map = {
        'drd2': "tdc", 
        'bbbp': "admet", 
        'mutagenicity': "admet", 
        'plogp': "gym", 
        'qed': "admet", 
        'ampa': "admet", 
        'carc': "admet", 
        'erg': "admet", 
        'hia': "admet", 
        'liver': "admet"
    }
    admet_map = {
        'bbbp': "BBB_Martins", 
        'mutagenicity': "AMES", 
        'qed': "QED", 
        'ampa': "PAMPA_NCATS", 
        'carc': "Carcinogens_Lagunin", 
        'erg': "hERG", 
        'hia': "HIA_Hou", 
        'liver': "DILI"
    }

    valid, smiles = validate_ligand(smiles, "smiles")
    if not valid : 
        raise ValueError("Invalid SMILES")
    mol = parse_ligand(smiles, "smiles")

    scores = {}
    admet_preds = admet_model.predict(smiles=smiles)
    for prop in properties.split('+'):
        prop_model = model_map[prop]
        prop_score = None
        if prop_model is None : 
            raise ValueError("Invalid property")
        if prop_model == 'gym' and prop == "plogp": 
            prop_score = reward_penalized_log_p(mol)
        elif prop_model == 'tdc' and prop == 'drd2' : 
            prop_score = drd2_model([smiles])[0]
        elif prop_model == 'admet':
            prop_score = admet_preds[admet_map[prop]]
        scores[prop] = prop_score

    return scores