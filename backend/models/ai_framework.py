from pydantic import BaseModel

class ConvertInputRequest(BaseModel):
    smiles: str
    properties: str

class ConvertInputResponse(BaseModel):
    output_path: str

class GenerationRequest(BaseModel):
    base_model: str = 'mistralai/Mistral-7B-Instruct-v0.3'
    lora_weights: str = 'NingLab/GeLLMO-P6-Mistral'
    task: str = 'abmp'
    setting: str = 'seen'
    data_path: str = 'Input_output/01_input/ALO_INPUT.json'
    output_dir: str = 'Input_output/02_output'
    num_return_sequences: int = 1

class GenerationResponse(BaseModel):
    output_path: str