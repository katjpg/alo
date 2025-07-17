import sys
import json
import os
from openai import OpenAI
from rdkit import Chem  # ✅ RDKit
from rdkit.Chem import rdFMCS

def extract_smiles_from_response(response):
    """Extract SMILES string from response text."""
    if "<SMILES>" in response and "</SMILES>" in response:
        start = response.find("<SMILES>") + len("<SMILES>")
        end = response.find("</SMILES>")
        return response[start:end].strip()
    else:
        return response.strip()

def canonicalize_smiles(smiles):
    """Convert SMILES to canonical form to ensure uniqueness."""
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        return Chem.MolToSmiles(mol, canonical=True)
    return None

def is_input_preserved(input_smiles, candidate_smiles, threshold=0.6):
    """Check if input molecule is largely preserved in the candidate molecule."""
    try:
        input_mol = Chem.MolFromSmiles(input_smiles)
        candidate_mol = Chem.MolFromSmiles(candidate_smiles)
        if not input_mol or not candidate_mol:
            return False

        mcs = rdFMCS.FindMCS([input_mol, candidate_mol])
        mcs_mol = Chem.MolFromSmarts(mcs.smartsString)
        if not mcs_mol:
            return False

        mcs_atoms = mcs_mol.GetNumAtoms()
        input_atoms = input_mol.GetNumAtoms()
        ratio = mcs_atoms / input_atoms
        return ratio >= threshold
    except Exception as e:
        print(f"⚠️ RDKit error: {e}")
        return False

def ai_check_preservation(input_smiles, candidate_smiles):
    """Use OpenAI to check similarity between input and candidate."""
    prompt = (
        f"Does the molecule '{candidate_smiles}' preserve the core structure of '{input_smiles}'? "
        f"Answer 'yes' if most of the original structure remains, otherwise 'no'."
    )
    try:
        response = client.chat.completions.create(
            model="gpt-3.5-turbo",  # 빠른 모델
            messages=[{"role": "user", "content": prompt}],
            temperature=0
        )
        reply = response.choices[0].message.content.strip().lower()
        return "yes" in reply
    except Exception as e:
        print(f"❌ OpenAI API Error: {e}")
        return False

def check_and_filter_responses(input_smiles, responses, output_path):
    """Filter responses and save only those preserving input structure."""
    filtered = []
    for idx, response in enumerate(responses, start=1):
        print(f"🔸 Checking response {idx}/{len(responses)}")
        candidate_smiles = extract_smiles_from_response(response)
        if not candidate_smiles:
            print("⚠️ Empty response, skipping.")
            continue

        canonical_smiles = canonicalize_smiles(candidate_smiles)
        if not canonical_smiles:
            print("⚠️ Invalid SMILES, skipping.")
            continue

        if is_input_preserved(input_smiles, canonical_smiles):
            filtered.append(canonical_smiles)
            print(f"✅ Preserved (RDKit): {canonical_smiles}")
        elif ai_check_preservation(input_smiles, canonical_smiles):
            filtered.append(canonical_smiles)
            print(f"✅ Preserved (OpenAI): {canonical_smiles}")
        else:
            print(f"❌ Filtered out: {canonical_smiles}")

    # ✅ 중복 제거 (canonical SMILES 기준)
    unique_filtered = sorted(set(filtered))
    duplicates_removed = len(filtered) - len(unique_filtered)
    print(f"🔹 Removed duplicates: {duplicates_removed} duplicates found.")

    result = {
        "filtered_responses": unique_filtered,
        "count": len(unique_filtered)
    }
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2)
    print(f"📦 Saved {len(unique_filtered)} unique responses to {output_path}")

def main(input_json_path, output_json_path, input_smiles):
    """Main function."""
    print(f"🔹 Loading input from: {input_json_path}")
    with open(input_json_path, "r", encoding="utf-8") as f:
        data = json.load(f)

    responses = data[0]["response"]
    print(f"🔹 Checking {len(responses)} responses for input: {input_smiles}")

    check_and_filter_responses(input_smiles, responses, output_json_path)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python filter_with_openai.py <input_json> <output_json> <input_smiles>")
        sys.exit(1)

    input_json_file = sys.argv[1]
    output_json_file = sys.argv[2]
    input_smiles_str = sys.argv[3]

    main(input_json_file, output_json_file, input_smiles_str)