from fastapi import APIRouter, HTTPException
from models.ai_framework import (
    ConvertInputRequest, ConvertInputResponse,
    GenerationRequest, GenerationResponse
)

from services.molecule_to_json import molecule_to_json
from services.inference import generate


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
    

