import streamlit as st
import pandas as pd
import numpy as np
import re
np.product = np.prod
from rdkit.Chem import AllChem, Draw
from rdkit import Chem
from mordred import Calculator, descriptors
from joblib import load

def All_Mordred_descriptors(data):
    calc = Calculator(descriptors, ignore_3D=False)
    mols = [Chem.MolFromSmiles(smi) for smi in data]
    df = calc.pandas(mols)
    return df

def is_allowed_smiles(smiles):
    if pd.isna(smiles) or not isinstance(smiles, str):
        return False
    if '.' in smiles:
        return False
    allowed_elements = {'H', 'B', 'C', 'N', 'O', 'S', 'P', 'F', 'CL', 'BR', 'I'}
    bracketed_symbols = re.findall(r'\[\d*([A-Z][a-z]?)', smiles)
    for sym in bracketed_symbols:
        if sym.upper() not in allowed_elements:
            return False
    outside = re.sub(r'\[[^\]]*\]', '', smiles)
    outside_symbols = re.findall(r'Cl|Br|cl|br|[A-Za-z]', outside)
    for sym in outside_symbols:
        if sym.upper() not in allowed_elements:
            return False
    return True

def get_canonical(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(mol, canonical=True) if mol else None

def predict_mlp(m1, m2, m3, m4, m5, X):
    s1 = m1.predict(X)
    s2 = m2.predict(X)
    s3 = m3.predict(X)
    s4 = m4.predict(X)
    s5 = m5.predict(X)
    arr_pred = np.zeros(s1.shape[0])
    for i in range(arr_pred.shape[0]):
        arr_pred[i] = (s1[i] + s2[i] + s3[i] + s4[i] + s5[i]) / 5
    return arr_pred

best_ens1 = load('model/mlp1_ens_saved.joblib')
best_ens2 = load('model/mlp2_ens_saved.joblib')
best_ens3 = load('model/mlp3_ens_saved.joblib')
best_ens4 = load('model/mlp4_ens_saved.joblib')
best_ens5 = load('model/mlp5_ens_saved.joblib')

st.sidebar.title('Welcome Scientist, Are you ready for magic ??')
st.sidebar.image("science.jpg")
st.sidebar.divider()
st.sidebar.write("CNRS")
st.sidebar.write("LRGP")
st.sidebar.divider()
st.sidebar.write("No reproduction is authorized")

st.title("Tool for prediction of thermodynamics properties using AI models")
st.header("Enter your smiles and wait for results.....")
st.subheader("Prediction of enthalpy, entropy and so on")

smiles = st.text_input('**Enter your smile:**:', value=None, max_chars=None)

col1, col2 = st.columns(2)

if col1.button("Submit"):
    st.write("Process in ...")
    smiles_can = get_canonical(smiles)
    if smiles_can and is_allowed_smiles(smiles_can):
        descripteur = pd.read_excel("liste_descripteurs_retenus.xlsx")
        descripteurs_list = descripteur.iloc[:, 0].tolist()
        df_all_descriptors = All_Mordred_descriptors([smiles_can])  # ← liste
        smile_final = df_all_descriptors[descripteurs_list]
        for col in smile_final.columns:  # ← correction bug col
            smile_final[col] = pd.to_numeric(smile_final[col], errors='coerce').fillna(0)
        result = predict_mlp(best_ens1, best_ens2, best_ens3, best_ens4, best_ens5, smile_final)
        results_kjmol = (result * (8.31446261815324 * 298.15)) / 1000
        mol = Chem.MolFromSmiles(smiles_can)
        img = Draw.MolToImage(mol)  # ← Draw maintenant importé
        st.image(img, caption="Représentation 2D")
        st.write('**Enthalpie de formation : %.4f kJ/mol**' % results_kjmol[0])
    else:
        st.error("SMILES invalide")

if col2.button("Reset"):
    st.rerun()

st.markdown("---")
st.caption("Merci d'avoir utilisé cette plateforme.......")
