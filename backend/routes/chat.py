from fastapi import APIRouter, HTTPException
from models.chat import ChatRequest, ChatResponse
from openai import OpenAI
from services.ligand_search import ligand_search
import json

client = OpenAI()
router = APIRouter(tags=["Chat"])

@router.post("/chat", response_model = ChatResponse)
def ai_chat(request: ChatRequest) : 
    messages = [
        {
            "role": "system",
            "content": ("Answer questions only with reputable sources (reliable scientific journal, .gov, and .edu source). Do not guess or provide unverified information.")
        },
        {
            "role": "user", 
            "content": request.input_message
        }
    ]
    tools = [{
        "type": "function",
        "name": "ligand_search",
        "description" : "Given a ligand, output information",
        "parameters":{
            "type": "object",
            "properties": {
                "smiles": {"type" : "string"}
            },
            "required": ["smiles"],
            "additionalProperties": False
        }
    }]
    response = client.responses.create(
        model="gpt-4.1", 
        input= messages, 
        tools=tools
    )
    
    tool_called = response.output
    if tool_called:
        tool_called = tool_called[0]
        if tool_called.type == "function_call" and tool_called.name == "ligand_search":
            args = json.loads(tool_called.arguments)
            result =  ligand_search(**args)
            return ChatResponse(chat_response=result)

    return ChatResponse(chat_response=response.output_text)
