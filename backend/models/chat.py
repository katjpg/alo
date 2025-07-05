from pydantic import BaseModel

class ChatResponse (BaseModel): 
    chat_response : str

class ChatRequest (BaseModel) : 
    input_message : str