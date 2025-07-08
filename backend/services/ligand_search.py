import requests
from urllib.parse import quote
from utils.parse_validate import validate_ligand
from openai import OpenAI

def get_chembl_id (smiles : str, data_link : str) -> str:
    encoded_smiles = quote(smiles)
    response = requests.get(f"{data_link}/molecule/{encoded_smiles}?format=json")
    if response.status_code == 200 : 
        response = response.json()
        chembl_id = response.get('molecule_chembl_id')
        return chembl_id
    else : 
        return f"Error: {response.status_code}, {response.text}"

def get_ligand_docs (ligand_id : str, data_link : str) -> dict:
    info = []
    response = requests.get(f"{data_link}/activity.json?molecule_chembl_id={ligand_id}")
    if response.status_code == 200 : 
        response = response.json()
        for activity in response.get('activities') : 
            doc_id = activity.get('document_chembl_id')
            doc_info = requests.get(f"{data_link}/document/{doc_id}.json")
            info.append(doc_info.text)
        return info
    else : 
        return f"Error: {response.status_code}, {response.text}"

def get_ligand_details (ligand_id : str, data_link : str) -> dict : 
    response = requests.get(f"{data_link}/molecule/{ligand_id}.json")
    if response.status_code == 200 : 
        response = response.json()
        return response
    else : 
        return f"Error: {response.status_code}, {response.text}"
    
def ligand_search (smiles : str) -> str:

    chembl_data_link = "https://www.ebi.ac.uk/chembl/api/data"

    # Standardize smiles
    _, smiles = validate_ligand(smiles, "smiles")
    if smiles is None : 
        return "Error : Invalid SMILES string"

    # Retrieve Ligand Information
    chembl_id = get_chembl_id(smiles, chembl_data_link)
    info = get_ligand_details(chembl_id, chembl_data_link)
    docs = get_ligand_docs(chembl_id, chembl_data_link)
    
    # Create summary from information outputted
    client = OpenAI()
    messages = [{
        "role": "system", 
        "content": "Using the documents given, analyze the molecule list and summarize the information into meaningful insights."
    },
    {
        "role": "user", 
        "content" : f"Molecule information : {info}\n Documents : {docs}"
    }]
    
    response = client.chat.completions.create(
        model="gpt-4.1", 
        messages= messages
    )

    return response.choices[0].message.content