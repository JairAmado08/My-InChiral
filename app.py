import streamlit as st
import itertools
import os
import zipfile
import tempfile
import io
import sys
import warnings
warnings.filterwarnings('ignore')

# Importación de RDKit
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    RDKIT_AVAILABLE = True
except ImportError:
    st.error("❌ RDKit no está instalado. Por favor instala RDKit para usar la funcionalidad de conversión a XYZ.")
    st.info("Instala con: pip install rdkit")
    RDKIT_AVAILABLE = False

# ==============================
# INICIO DE ESTILOS PERSONALIZADOS
# ==============================
def aplicar_estilos():
    st.markdown(
        """
        <style>
        /* Fondo general */
        .stApp {
            background: linear-gradient(120deg, #0f2027, #203a43, #2c5364);
            color: #ffffff;
            font-family: 'Segoe UI', sans-serif;
        }

        /* Títulos */
        h1, h2, h3, h4 {
            color: #00e6b8;
            text-align: center;
            font-weight: 700;
        }

        /* Caja de texto */
        .stTextInput > div > div > input {
            background-color: #1c1c1c;
            color: #00ffcc;
            border: 1px solid #00e6b8;
            border-radius: 10px;
        }

        /* Botones */
        .stButton > button {
            background: linear-gradient(90deg, #00e6b8, #0099ff);
            color: white;
            font-weight: bold;
            border-radius: 10px;
            padding: 0.6em 1.2em;
            transition: all 0.3s ease-in-out;
        }
        .stButton > button:hover {
            background: linear-gradient(90deg, #0099ff, #00e6b8);
            transform: scale(1.05);
        }

        /* Tabs */
        .stTabs [data-baseweb="tab-list"] {
            gap: 10px;
        }
        .stTabs [data-baseweb="tab"] {
            background-color: #1c1c1c;
            color: #00e6b8;
            border-radius: 10px;
            padding: 8px 16px;
        }
        .stTabs [aria-selected="true"] {
            background-color: #00e6b8;
            color: black;
        }

        /* Cajas de información */
        .stAlert {
            border-radius: 10px;
        }

        /* Pie de página */
        footer {
            visibility: hidden;
        }
        </style>
        """,
        unsafe_allow_html=True
    )

# ==============================
# APP PRINCIPAL
# ==============================
def main():
    st.set_page_config(
        page_title="Inchiral - Generador de Estereoisómeros",
        page_icon="🧬",
        layout="wide",
        initial_sidebar_state="expanded"
    )

    aplicar_estilos()

    st.title("🧬 Generador de Estereoisómeros")
    st.markdown("<h4>✨ Genera todos los estereoisómeros posibles y convierte a formato XYZ ✨</h4>", unsafe_allow_html=True)

    # ---------------- Sidebar ----------------
    with st.sidebar:
        try:
            st.image("imagenes1/inchiral final.png", width=200)
        except:
            st.markdown("**🧬 Inchiral**")
        st.markdown("---")
        st.markdown("## ℹ️ Información")
        st.info("""
        **Instrucciones:**
        1. Ingresa un código SMILES (con o sin quiralidad especificada)
        2. El sistema detecta automáticamente si la molécula es quiral
        3. Si tiene centros quirales especificados (@ o @@), genera todos los estereoisómeros
        4. Máximo 3 centros quirales para evitar demasiados isómeros
        5. Opcionalmente convierte a formato XYZ para visualización 3D
        """)

    # ---------------- Entrada ----------------
    st.subheader("📝 Entrada de Datos")
    smiles_input = st.text_input(
        "👉 Ingresa el código SMILES:",
        placeholder="Ejemplo: C[C@H](O)[C@@H](N)C"
    )

    # TODO: Mantener aquí toda la lógica que ya tienes
    # (detectar_quiralidad, generar_estereoisomeros, smiles_to_xyz, etc.)
    # sin modificar nada en las funciones, solo la estética visual.

    st.markdown("---")
    st.markdown(
        """
        <div style='text-align: center; color: #bbb;'>
            <small>🧬 <strong>Inchiral</strong> - Universidad Científica del Sur<br>
            Generador de Estereoisómeros | Desarrollado con Streamlit y RDKit</small>
        </div>
        """,
        unsafe_allow_html=True
    )

if __name__ == "__main__":
    main()
