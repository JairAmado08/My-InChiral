import streamlit as st
import itertools
import os
import zipfile
import tempfile
import io
import sys

# Configuración para evitar warnings de RDKit
import warnings
warnings.filterwarnings('ignore')

# Manejo de importación de RDKit
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    RDKIT_AVAILABLE = True
except ImportError:
    st.error("❌ RDKit no está instalado. Por favor instala RDKit para usar la funcionalidad de conversión a XYZ.")
    st.info("Instala con: pip install rdkit")
    RDKIT_AVAILABLE = False

# Configuración de CSS personalizado
def load_custom_css():
    st.markdown("""
    <style>
    /* Tema principal con gradiente */
    .main {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        padding: 0;
    }
    
    /* Sidebar personalizado */
    .css-1d391kg {
        background: linear-gradient(180deg, #2D3748 0%, #1A202C 100%);
        border-right: 3px solid #4FD1C7;
    }
    
    /* Títulos principales */
    .main-title {
        background: linear-gradient(45deg, #4FD1C7, #63B3ED);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        font-size: 3rem;
        font-weight: 800;
        text-align: center;
        margin-bottom: 1rem;
        text-shadow: 0 0 20px rgba(79, 209, 199, 0.3);
    }
    
    .subtitle {
        color: #FFFFFF;
        text-align: center;
        font-size: 1.2rem;
        margin-bottom: 2rem;
        font-weight: 400;
        text-shadow: 0 2px 4px rgba(0, 0, 0, 0.3);
        background: linear-gradient(45deg, #4FD1C7, #63B3ED);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        filter: drop-shadow(0 0 10px rgba(255, 255, 255, 0.5));
    }
    
    /* Tarjetas de información */
    .info-card {
        background: rgba(255, 255, 255, 0.95);
        padding: 1.5rem;
        border-radius: 15px;
        box-shadow: 0 10px 30px rgba(0, 0, 0, 0.1);
        border: 2px solid transparent;
        background-clip: padding-box;
        backdrop-filter: blur(10px);
        margin: 1rem 0;
        transition: all 0.3s ease;
    }
    
    .info-card:hover {
        transform: translateY(-5px);
        box-shadow: 0 20px 40px rgba(0, 0, 0, 0.15);
    }
    
    /* Botones personalizados */
    .stButton > button {
        background: linear-gradient(45deg, #4FD1C7, #63B3ED);
        color: white;
        border: none;
        border-radius: 25px;
        padding: 0.75rem 2rem;
        font-weight: 600;
        text-transform: uppercase;
        letter-spacing: 1px;
        transition: all 0.3s ease;
        box-shadow: 0 5px 15px rgba(79, 209, 199, 0.3);
    }
    
    .stButton > button:hover {
        transform: translateY(-2px);
        box-shadow: 0 10px 25px rgba(79, 209, 199, 0.4);
    }
    
    /* Input personalizado */
    .stTextInput > div > div > input {
        border-radius: 15px;
        border: 2px solid #E2E8F0;
        padding: 1rem;
        font-size: 1.1rem;
        transition: all 0.3s ease;
    }
    
    .stTextInput > div > div > input:focus {
        border-color: #4FD1C7;
        box-shadow: 0 0 0 3px rgba(79, 209, 199, 0.1);
    }
    
    /* Tabs personalizados */
    .stTabs [data-baseweb="tab-list"] {
        background: rgba(255, 255, 255, 0.1);
        border-radius: 15px;
        backdrop-filter: blur(10px);
    }
    
    .stTabs [data-baseweb="tab"] {
        color: #E2E8F0;
        font-weight: 600;
        padding: 1rem 2rem;
        border-radius: 10px;
        transition: all 0.3s ease;
    }
    
    .stTabs [aria-selected="true"] {
        background: linear-gradient(45deg, #4FD1C7, #63B3ED);
        color: white;
    }
    
    /* Código personalizado */
    .stCode {
        background: linear-gradient(135deg, #2D3748, #4A5568);
        border: 1px solid #4FD1C7;
        border-radius: 10px;
        color: #E2E8F0;
    }
    
    /* Alertas personalizadas */
    .stAlert {
        border-radius: 15px;
        backdrop-filter: blur(10px);
    }
    
    /* Progress bar */
    .stProgress > div > div > div {
        background: linear-gradient(45deg, #4FD1C7, #63B3ED);
        border-radius: 10px;
    }
    
    /* Footer */
    .footer {
        background: linear-gradient(45deg, #2D3748, #4A5568);
        color: #E2E8F0;
        padding: 2rem;
        border-radius: 15px;
        text-align: center;
        margin-top: 3rem;
        border: 1px solid #4FD1C7;
    }
    
    /* Sidebar content */
    .sidebar-content {
        color: #E2E8F0;
    }
    
    /* Métricas personalizadas */
    .metric-card {
        background: linear-gradient(135deg, rgba(79, 209, 199, 0.1), rgba(99, 179, 237, 0.1));
        padding: 1rem;
        border-radius: 10px;
        border: 1px solid rgba(79, 209, 199, 0.3);
        text-align: center;
        margin: 0.5rem 0;
    }
    
    /* Animaciones */
    @keyframes glow {
        0% { box-shadow: 0 0 5px rgba(79, 209, 199, 0.3); }
        50% { box-shadow: 0 0 20px rgba(79, 209, 199, 0.6); }
        100% { box-shadow: 0 0 5px rgba(79, 209, 199, 0.3); }
    }
    
    .glow-animation {
        animation: glow 2s infinite;
    }
    </style>
    """, unsafe_allow_html=True)

def detectar_quiralidad(smiles: str):
    if not RDKIT_AVAILABLE:
        return False, "RDKit no disponible", []
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, "SMILES inválido", []
        
        centros = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        
        if len(centros) == 0:
            return False, "Su molécula no es quiral", []
        else:
            return True, f"Su molécula es quiral. Se detectaron {len(centros)} posibles centros", centros
            
    except Exception as e:
        return False, f"Error al analizar la molécula: {str(e)}", []

def analizar_centros_existentes(smiles: str):
    centros_especificados = 0
    posiciones_at = []
    i = 0
    
    while i < len(smiles):
        if smiles[i] == "@":
            if i + 1 < len(smiles) and smiles[i+1] == "@":
                centros_especificados += 1
                posiciones_at.append(i)
                i += 2
            else:
                centros_especificados += 1
                posiciones_at.append(i)
                i += 1
        else:
            i += 1
    
    return centros_especificados, posiciones_at

def generar_estereoisomeros(smiles: str):
    posiciones = []
    i = 0
    while i < len(smiles):
        if smiles[i] == "@":
            if i + 1 < len(smiles) and smiles[i+1] == "@":
                posiciones.append((i, True))  # ya es @@
                i += 2
            else:
                posiciones.append((i, False))  # es @ simple
                i += 1
        else:
            i += 1
    
    n = len(posiciones)
    
    if n == 0:
        st.warning("⚠️ El SMILES no tiene centros quirales especificados con @ o @@. No se generarán isómeros.")
        return [], n
    elif n > 3:
        st.error("❌ El SMILES tiene más de 3 centros quirales. No se generarán isómeros.")
        return [], n
    
    combinaciones = list(itertools.product(["@", "@@"], repeat=n))
    resultados = []
    
    for comb in combinaciones:
        chars = list(smiles)
        offset = 0
        for (pos, era_doble), val in zip(posiciones, comb):
            real_pos = pos + offset
            if era_doble:
                chars[real_pos:real_pos+2] = list(val)
                offset += len(val) - 2
            else:
                chars[real_pos:real_pos+1] = list(val)
                offset += len(val) - 1
        resultados.append("".join(chars))
    
    return resultados, n

def smiles_to_xyz(smiles, mol_id):
    if not RDKIT_AVAILABLE:
        return None, "❌ RDKit no está disponible"
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None, f"❌ Error: SMILES inválido {smiles}"
        
        mol = Chem.AddHs(mol)
        
        params = AllChem.ETKDGv3()
        params.randomSeed = 42  
        
        embed_result = AllChem.EmbedMolecule(mol, params)
        if embed_result != 0:
            params.useRandomCoords = True
            embed_result = AllChem.EmbedMolecule(mol, params)
            if embed_result != 0:
                return None, f"⚠️ No se pudo generar conformación 3D para {smiles}"
        
        try:
            if AllChem.MMFFHasAllMoleculeParams(mol):
                AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
            else:
                AllChem.UFFOptimizeMolecule(mol, maxIters=500)
        except:
            pass
        
        conf = mol.GetConformer()
        xyz_content = f"{mol.GetNumAtoms()}\n{smiles}\n"
        
        for atom in mol.GetAtoms():
            pos = conf.GetAtomPosition(atom.GetIdx())
            xyz_content += f"{atom.GetSymbol()} {pos.x:.4f} {pos.y:.4f} {pos.z:.4f}\n"
        
        return xyz_content, f"✅ Molécula {mol_id} procesada correctamente"
        
    except Exception as e:
        return None, f"❌ Error procesando {smiles}: {str(e)}"

def crear_archivo_zip(archivos_xyz):
    zip_buffer = io.BytesIO()
    with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
        for filename, content in archivos_xyz.items():
            zip_file.writestr(filename, content)
    return zip_buffer.getvalue()

def main():
    st.set_page_config(
        page_title="Inchiral - Generador de Estereoisómeros",
        page_icon="🧬",
        layout="wide",
        initial_sidebar_state="expanded"
    )
    
    # Cargar CSS personalizado
    load_custom_css()
    
    # Título principal con estilo
    st.markdown('<h1 class="main-title">🧬 INCHIRAL</h1>', unsafe_allow_html=True)
    st.markdown('<p class="subtitle">Generador Avanzado de Estereoisómeros</p>', unsafe_allow_html=True)
    
    # Sidebar mejorado
    with st.sidebar:
        st.markdown('<div class="sidebar-content">', unsafe_allow_html=True)
        
        try:
            st.image("imagenes1/inchiralucsur.png", width=200)
        except:
            st.markdown("""
            <div style='text-align: center; padding: 2rem;'>
                <div style='font-size: 4rem; margin-bottom: 1rem;'>🧬</div>
                <h2 style='color: #4FD1C7; margin: 0;'>INCHIRAL</h2>
            </div>
            """, unsafe_allow_html=True)
        
        st.markdown("---")
        
        st.markdown("""
        <div class="info-card">
        <h3 style='color: #2D3748; margin-top: 0;'>📋 Instrucciones</h3>
        <ol style='color: #4A5568; line-height: 1.8;'>
            <li>Ingresa un código SMILES (con o sin quiralidad)</li>
            <li>Sistema detecta automáticamente centros quirales</li>
            <li>Genera todos los estereoisómeros posibles</li>
            <li>Máximo 3 centros quirales procesables</li>
            <li>Conversión opcional a formato XYZ 3D</li>
        </ol>
        </div>
        """, unsafe_allow_html=True)
        
        st.markdown("""
        <div class="info-card">
        <h3 style='color: #2D3748; margin-top: 0;'>💡 Ejemplos de SMILES</h3>
        <div style='color: #4A5568; line-height: 1.8;'>
            <strong>Sin quiralidad:</strong><br>
            <code>CCO</code><br><br>
            <strong>Molécula quiral:</strong><br>
            <code>CC(O)C(N)C</code><br><br>
            <strong>Con quiralidad:</strong><br>
            <code>C[C@H](O)[C@@H](N)C</code><br><br>
            <strong>Aminoácido:</strong><br>
            <code>N[C@@H](C)C(=O)O</code>
        </div>
        </div>
        """, unsafe_allow_html=True)
        
        st.markdown('</div>', unsafe_allow_html=True)
    
    # Sección de entrada de datos
    st.markdown("""
    <div class="info-card">
    <h2 style='color: #2D3748; margin-top: 0;'>📝 Entrada de Datos</h2>
    """, unsafe_allow_html=True)
    
    smiles_input = st.text_input(
        "👉 Ingresa el código SMILES:",
        placeholder="Ejemplo: C[C@H](O)[C@@H](N)C",
        help="Introduce tu molécula en formato SMILES"
    )
    st.markdown('</div>', unsafe_allow_html=True)
    
    if smiles_input:
        # Análisis de quiralidad
        st.markdown("""
        <div class="info-card">
        <h2 style='color: #2D3748; margin-top: 0;'>🔍 Análisis de Quiralidad</h2>
        """, unsafe_allow_html=True)
        
        es_quiral, mensaje_quiralidad, centros_detectados = detectar_quiralidad(smiles_input)
        centros_especificados, posiciones_at = analizar_centros_existentes(smiles_input)

        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("""
            <div class="metric-card">
            <h3 style='color: #2D3748; margin-top: 0;'>🔎 Análisis con RDKit</h3>
            """, unsafe_allow_html=True)
            
            if RDKIT_AVAILABLE:
                if es_quiral:
                    st.success(f"✅ {mensaje_quiralidad}")
                    if centros_detectados:
                        st.write("**Centros detectados:**")
                        for i, (idx, chirality) in enumerate(centros_detectados):
                            tipo_quiralidad = str(chirality) if chirality else "Sin asignar"
                            st.write(f"• Átomo {idx}: {tipo_quiralidad}")
                else:
                    if "inválido" in mensaje_quiralidad:
                        st.error(f"❌ {mensaje_quiralidad}")
                    else:
                        st.warning(f"⚠️ {mensaje_quiralidad}")
            else:
                st.warning("⚠️ RDKit no disponible para análisis")
            
            st.markdown('</div>', unsafe_allow_html=True)
        
        with col2:
            st.markdown("""
            <div class="metric-card">
            <h3 style='color: #2D3748; margin-top: 0;'>📋 Centros Especificados</h3>
            """, unsafe_allow_html=True)
            
            if centros_especificados > 0:
                st.success(f"✅ {centros_especificados} centros con @ o @@ especificados")
                for pos in posiciones_at:
                    st.write(f"• Posición {pos}")
            else:
                st.warning("⚠️ No hay centros especificados con @ o @@")
                
            st.markdown('</div>', unsafe_allow_html=True)
        
        st.markdown('</div>', unsafe_allow_html=True)
        
        if RDKIT_AVAILABLE and es_quiral and centros_especificados == 0:
            st.info("""
            💡 Tu molécula es quiral pero no tiene centros especificados con @ o @@.
            Ejemplo: `CC(O)C(N)C` → `C[C@H](O)[C@@H](N)C`
            """)
        
        # Generación de estereoisómeros
        isomeros, n_centros = [], 0
        if centros_especificados > 0:
            with st.spinner("🔄 Generando estereoisómeros..."):
                isomeros, n_centros = generar_estereoisomeros(smiles_input)
        
        # Tabs mejorados
        tab1, tab2, tab3 = st.tabs(["📋 Lista Completa", "💾 Descargar SMI", "🧪 Convertir a XYZ"])
        
        if isomeros:
            with tab1:
                st.markdown("""
                <div class="info-card">
                <h3 style='color: #2D3748; margin-top: 0;'>🧪 Estereoisómeros Generados</h3>
                """, unsafe_allow_html=True)
                
                col1, col2 = st.columns(2)
                for i, isomero in enumerate(isomeros):
                    if i % 2 == 0:
                        col1.code(f"{i+1}. {isomero}", language="text")
                    else:
                        col2.code(f"{i+1}. {isomero}", language="text")
                
                st.markdown('</div>', unsafe_allow_html=True)
            
            with tab2:
                st.markdown("""
                <div class="info-card">
                <h3 style='color: #2D3748; margin-top: 0;'>💾 Descarga de Archivos SMI</h3>
                """, unsafe_allow_html=True)
                
                smi_content = "\n".join(isomeros)
                st.download_button(
                    label="📥 Descargar archivo.smi",
                    data=smi_content,
                    file_name="estereoisomeros.smi",
                    mime="text/plain"
                )
                
                with st.expander("👀 Vista previa del archivo SMI"):
                    st.code(smi_content, language="text")
                    
                st.markdown('</div>', unsafe_allow_html=True)
            
            with tab3:
                st.markdown("""
                <div class="info-card">
                <h3 style='color: #2D3748; margin-top: 0;'>🧪 Conversión a Formato XYZ</h3>
                """, unsafe_allow_html=True)
                
                if st.button("🚀 Convertir todos a XYZ", type="primary"):
                    progress_bar = st.progress(0)
                    status_text = st.empty()
                    archivos_xyz = {}
                    mensajes_log = []
                    
                    for i, smiles in enumerate(isomeros):
                        try:
                            progress = (i + 1) / len(isomeros)
                            progress_bar.progress(progress)
                            status_text.text(f"Procesando molécula {i+1}/{len(isomeros)}: {smiles}")
                            
                            xyz_content, mensaje = smiles_to_xyz(smiles, i+1)
                            mensajes_log.append(mensaje)
                            
                            if xyz_content:
                                archivos_xyz[f"mol_{i+1}.xyz"] = xyz_content
                        except Exception as e:
                            mensajes_log.append(f"❌ Error procesando molécula {i+1}: {str(e)}")
                    
                    progress_bar.progress(1.0)
                    status_text.text("✅ Proceso completado!")
                    
                    with st.expander("📋 Log de procesamiento"):
                        for mensaje in mensajes_log:
                            if "❌" in mensaje or "⚠️" in mensaje:
                                st.error(mensaje)
                            else:
                                st.success(mensaje)
                    
                    if archivos_xyz:
                        zip_data = crear_archivo_zip(archivos_xyz)
                        st.download_button(
                            label="📦 Descargar archivos XYZ (ZIP)",
                            data=zip_data,
                            file_name="estereoisomeros_xyz.zip",
                            mime="application/zip"
                        )
                        with st.expander("👀 Vista previa del primer archivo XYZ"):
                            primer_archivo = list(archivos_xyz.values())[0]
                            st.code(primer_archivo, language="text")
                
                st.markdown('</div>', unsafe_allow_html=True)
        else:
            st.info("💡 Ingresa un SMILES con centros quirales especificados (@ o @@) para generar estereoisómeros")
    
    # Footer mejorado
    st.markdown("""
    <div class="footer">
        <div style='font-size: 2rem; margin-bottom: 1rem;'>🧬</div>
        <h3 style='color: #4FD1C7; margin: 0.5rem 0;'>INCHIRAL</h3>
        <p style='margin: 0.5rem 0; opacity: 0.8;'>Universidad Científica del Sur</p>
        <p style='margin: 0; font-size: 0.9rem; opacity: 0.7;'>
            Generador Avanzado de Estereoisómeros | Desarrollado con Streamlit y RDKit
        </p>
    </div>
    """, unsafe_allow_html=True)

if __name__ == "__main__":
    main()
