from fastapi import APIRouter, HTTPException
from models.chat import ChatRequest, ChatResponse
from openai import OpenAI

client = OpenAI()
router = APIRouter(tags=["Chat"])

@router.post("/chat", response_model = ChatResponse)
def ai_chat(request: ChatRequest) : 
    messages = [
        {
            "role": "system",
            "content": ("Answer questions only with reputable sources. If you cannot find a reliable scientific journal, .gov, or .edu source, do not guess or provide unverified information.")
        },
        {
            "role": "user", 
            "content": request.input_message
        }
    ]
    response = client.chat.completions.create(
        model="gpt-4.1", 
        messages= messages, 

    )
    
    return ChatResponse(chat_response=response.choices[0].message.content)
