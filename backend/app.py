from fastapi import FastAPI
from routes import ligand

app = FastAPI()
app.include_router(ligand.router, prefix="/api")