import os
import sys
import json
import random
from datetime import datetime
from datasets import Dataset

def molecule_to_json(source_smiles: str, property_comb: str, output_dir: str = "Input_output/01_input") -> str:
    """
    Save a single molecule into a JSON file formatted for GeLLMO-C.
    The JSON will always be saved to Input_output/01_input folder.
    """
    # List of 10 property keys
    all_properties = ['drd2', 'bbbp', 'mutagenicity', 'plogp', 'qed', 'ampa', 'carc', 'erg', 'hia', 'liver']
    properties_dict = {prop: None for prop in all_properties}
    random_instr_idx = random.randint(0, 5)

    # Build data point
    my_data = [{
        'source_smiles': source_smiles,
        'target_smiles': None,
        'property_comb': property_comb,
        'improved': [],
        'stable': [],
        'swap_needed': None,
        'task': f'C:{property_comb}',
        'properties': properties_dict,
        'split': 'test',
        'instr_idx': random_instr_idx,
        'instr_setting': 'seen'
    }]

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    current_time = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_path = os.path.join(output_dir, f"ALO_INPUT_{current_time}.json")

    # Wrap data with top-level "test" key
    wrapped_data = {"test": my_data}

    # Save JSON
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(wrapped_data, f, indent=2)
    print(f"JSON file saved to: {output_path}")
    return output_path
