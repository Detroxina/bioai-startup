from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, QED # <--- IMPORTAMOS QED

app = FastAPI()

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

class SolicitudAnalisis(BaseModel):
    smiles: str

@app.post("/analizar")
def analizar_molecula(datos: SolicitudAnalisis):
    # 1. Crear Molécula
    mol = Chem.MolFromSmiles(datos.smiles)
    if mol is None:
        return {"error": "Estructura inválida"}

    # 2. CALCULAR QED (Tu porcentaje de calidad real)
    # QED devuelve un número entre 0 y 1 (ej: 0.85). Lo multiplicamos por 100.
    score_qed = QED.qed(mol)
    porcentaje_exito = round(score_qed * 100, 1)

    # 3. Reglas de Lipinski (Datos técnicos)
    peso = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    violaciones = 0
    if peso > 500: violaciones += 1
    if logp > 5: violaciones += 1
    if Lipinski.NumHDonors(mol) > 5: violaciones += 1
    if Lipinski.NumHAcceptors(mol) > 10: violaciones += 1

    # 4. Veredicto
    es_viable = violaciones <= 1
    mensaje = f"Calidad QED: {porcentaje_exito}%"
    
    # Elegir color de proteína según el Score
    pdb_id = "1CRN" if porcentaje_exito > 50 else "4HHB"

    return {
        "datos_reales": {
            "peso": round(peso, 2),
            "logp": round(logp, 2),
            "violaciones": violaciones
        },
        "resultado": {
            "es_viable": es_viable,
            "score_qed": porcentaje_exito, # <--- Enviamos el porcentaje real
            "mensaje": mensaje,
            "pdb_id": pdb_id
        }
    }