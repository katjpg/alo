# ALO Backend Integration

## Getting Started

```bash

git clone https://github.com/katjpg/alo.git
cd alo/backend

# Setting up venv
python3.11 -m venv venv
source venv/bin/activate  # Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Run the server
uvicorn app:app --reload
```

## Project Structure

```
backend/
├── app.py              # FastAPI application entry point
├── models/             # Pydantic models for requests/responses
├── routes/             # API route definitions
├── services/           # Application logic (property calculations)
├── utils/              # Utilities (parsing, validation, drawing)
└── requirements.txt    # Python dependencies
```

## API Endpoints

- `POST /properties/` - Calculate molecular properties from SMILES
- `POST /ligand/validate` - Validate ligand data (SMILES/SDF/MOL)
- `POST /ligand/draw` - Generate molecular structure visualization