from fastapi import FastAPI
from routes import ligand, chat

app = FastAPI()
app.include_router(ligand.router)
app.include_router(chat.router)