from fastapi import APIRouter, HTTPException
from models.ai_framework import (
    ConvertInputRequest, ConvertInputResponse,
    GenerationRequest, GenerationResponse,
    FilterResultRequest, FilterResultResponse
)

from services.molecule_to_json import molecule_to_json
from services.inference_cpu import generate
from services.filter_with_openai import filter_result
from utils.parse_validate import validate_ligand


router = APIRouter(tags=["AI-Framework"])

@router.post("/ai_framework/convert_input", response_model=ConvertInputResponse)
def convert_to_json(request: ConvertInputRequest):
    path = molecule_to_json(request.smiles, request.properties)
    return ConvertInputResponse(output_path=path)

@router.post("/ai_framework/generate", response_model=GenerationResponse)
def generate_molecule(request: GenerationRequest):
    path = generate(
                base_model=request.base_model,
                lora_weights=request.lora_weights,
                task=request.task,
                setting=request.setting, 
                data_path=request.data_path,
                output_dir=request.output_dir,
                num_return_sequences=request.num_return_sequences
            )
    return GenerationResponse(output_path=path)

@router.post("/ai_framework/filter_result", response_model=FilterResultResponse)
def filter(request: FilterResultRequest):
    valid, smiles = validate_ligand(request.smiles, "smiles")
    if not valid : 
        raise HTTPException(status_code=400, detail="Invalid SMILES")
    path = filter_result(request.input_data, request.output_data, request.smiles)
    return FilterResultResponse(output_path=path)
    

