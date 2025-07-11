from fastapi import FastAPI
from routes import ligand, chembl, accessibility

app = FastAPI()
app.include_router(ligand.router, prefix="/api")
app.include_router(chembl.router, prefix="/api")
app.include_router(accessibility.router, prefix="/api")