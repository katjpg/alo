from fastapi import FastAPI
from routes import ligand, chat, ai_framework

app = FastAPI()
app.include_router(ligand.router)
app.include_router(chat.router)
app.include_router(ai_framework.router)