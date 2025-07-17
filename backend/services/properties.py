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

admet_model = ADMETModel()

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

def property_score (smiles: str, properties: list) -> dict : 
    valid, smiles = validate_ligand(smiles, "smiles")
    if not valid : 
        raise ValueError("Invalid SMILES")
    mol = parse_ligand(smiles, "smiles")

    scores = {}

    if "plogp" in properties: 
        plogp = reward_penalized_log_p(mol)
        scores["plogp"] = plogp

    # DRD2 

    # what if some properties returned are longer sentences?
    # change properties from type dict to type string with '+' separating each prop (same format as molecule_to_json)???
    preds = admet_model.predict(smiles=smiles)
    for key, value in preds.items() : 
        if key in properties : 
            scores[key] = value

    return scores