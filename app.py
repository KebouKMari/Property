import streamlit as st
import pandas as pd
import numpy as np
np.product = np.prod
from rdkit.Chem import AllChem
from rdkit import Chem
from mordred import Calculator, descriptors
import mordred
from joblib import load

st.markdown(
    """
    <style>
    .stApp { background-color: black; }
    section[data-testid="stSidebar"] { background-color: #111; }
    </style>
    """,
    unsafe_allow_html=True
)

# Transforme les smiles en descripteurs
def All_Mordred_descriptors(data):
    calc = Calculator(descriptors, ignore_3D=False)
    mols = [Chem.MolFromSmiles(smi) for smi in data]
    df = calc.pandas(mols)
    return df

#
# Filtrage des smiles
#
def is_allowed_smiles(smiles):
    """
    Filtre strictement les SMILES pour ne garder que les atomes :
    H, B, C, N, O, S, P, F, Cl, Br, I.
    Gère les isotopes (ex: [3He]), les aromatiques (ex: c, n) et les charges.
    Exclut également les structures déconnectées (contenant un point '.').
    """
    if pd.isna(smiles) or not isinstance(smiles, str):
        return False
    
    # Filtrage des structures déconnectées (sels, mélanges) ---
    if '.' in smiles:
        return False
    
    # Liste officielle des symboles autorisés (insensible à la casse pour l'aromaticité)
    allowed_elements = {'H', 'B', 'C', 'N', 'O', 'S', 'P', 'F', 'CL', 'BR', 'I'}
    
    # 1. Analyse des atomes dans les crochets [ ]
    # Cette regex ignore les chiffres (isotopes), les signes +/-, etc., pour extraire le symbole.
    bracketed_symbols = re.findall(r'\[\d*([A-Z][a-z]?)', smiles)
    for sym in bracketed_symbols:
        if sym.upper() not in allowed_elements:
            return False # Capture 'He' même s'il y a un '3' devant.

    # 2. Analyse des atomes hors crochets
    # On retire d'abord tout ce qui est entre crochets pour ne pas fausser l'analyse
    outside = re.sub(r'\[[^\]]*\]', '', smiles)
    # On cherche les symboles à 2 lettres (Cl, Br) puis les lettres simples (C, N, O, c, n...)
    # On ignore les chiffres de structure (cycles) et les liaisons (=, #)
    outside_symbols = re.findall(r'Cl|Br|cl|br|[A-Za-z]', outside)
    
    for sym in outside_symbols:
        if sym.upper() not in allowed_elements:
            return False
                
    return True

#
# Conversion smile <-> smile canonique - <-> InchiKey
#
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


# Essayons voir si le meilleur model capture le comportement physique des molecules - MLP
best_ens1 = load('model/mlp1_ens_saved.joblib')
best_ens2 = load('model/mlp2_ens_saved.joblib')
best_ens3 = load('model/mlp3_ens_saved.joblib')
best_ens4 = load('model/mlp4_ens_saved.joblib')
best_ens5 = load('model/mlp5_ens_saved.joblib')


# Faire la  barre 
st.sidebar.title('Welcome Scientist, Hope you are doing well')
st.sidebar.image("science.jpg")

#ligne 
st.sidebar.divider()
st.sidebar.write("CNRS")
st.sidebar.write("LRGP")

st.sidebar.divider()

st.sidebar.write("This Work provided from the phd thesis Work, so no reproduction is authorized")

#Main Menu
st.title("Tool for prediction of thermodynamics properties using AI models")
st.header("Enter your smiles and wait for results.....")
st.subheader("Prediction of Enthalpy, entropy, LIE and many others")

# Lets 'Go

smiles = st.text_input('**Insert your smiles:**:', value= None, max_chars = None)

col1, col2 = st.columns(2)

if col1.button("Submit"):
    st.write("Process in ...")
    #Pretraitement 
    smiles_can = get_canonical(smiles)
    if is_allowed_smiles(smiles_can):
        descripteur = pd.read_excel("liste_descripteurs_retenus.xlsx")
        desc
        df_all_descriptors = All_Mordred_descriptors(smiles_can)
        smile_final = df_all_descriptors[descripteur]
        smile_final[col] = pd.to_numeric(smile_final[col], errors='coerce').fillna(0)
        result = predict_mlp(best_ens1, best_ens2, best_ens3, best_ens4, best_ens5, smile_final)
        results = (results * (8.31446261815324*298.15))/1000
        mol = Chem.MolFromSmiles(smiles_can)
        img = Draw.MolToImage(mol) 
        st.image(img, caption="Représentation 2D")
        text = st.text_area(label = result, height = 70)
        st.write('Enthalpie de formation en kJ/mol **%s**' % text)     

    else: 
        st.error("SMILES invalide")
    
if col2.button("Reset"):
    st.rerun()

st.markdown("---")
st.caption("Merci d'avoir utilisé la plateforme")
