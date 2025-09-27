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
    
    # Logo principal grande con estilo
    st.markdown("""
    <style>
    .main-logo {
        display: flex;
        justify-content: center;
        width: 100%;
        margin-bottom: 2rem;
        animation: float 3s ease-in-out infinite;
    }
    .main-logo img {
        max-width: 400px;
        height: auto;
        border-radius: 20px;
        box-shadow: 0 10px 30px rgba(79, 209, 199, 0.4);
        transition: transform 0.3s ease;
        filter: drop-shadow(0 0 20px rgba(79, 209, 199, 0.3));
    }
    .main-logo img:hover {
        transform: scale(1.02);
        box-shadow: 0 15px 40px rgba(79, 209, 199, 0.6);
    }
    .main-logo-fallback {
        text-align: center;
        animation: float 3s ease-in-out infinite;
        margin-bottom: 2rem;
    }
    .main-logo-fallback .emoji {
        font-size: 8rem;
        margin-bottom: 1rem;
        display: block;
        filter: drop-shadow(0 0 20px rgba(79, 209, 199, 0.5));
    }
    .main-logo-fallback h1 {
        background: linear-gradient(45deg, #4FD1C7, #63B3ED);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        font-size: 4rem;
        font-weight: 800;
        margin: 0;
        text-shadow: 0 0 30px rgba(79, 209, 199, 0.5);
    }
    @keyframes float {
        0%, 100% { transform: translateY(0px); }
        50% { transform: translateY(-15px); }
    }
    </style>
    <div class="main-logo">
        <img src="https://raw.githubusercontent.com/JairAmado08/My-InChiral/main/imagenes1/inchiralucsur.png" alt="Inchiral Logo">
    </div>
    """, unsafe_allow_html=True)
    
    st.markdown('<p class="subtitle">Generador Avanzado de Estereoisómeros</p>', unsafe_allow_html=True)
    
    # Sidebar mejorado
    with st.sidebar:
        st.markdown('<div class="sidebar-content">', unsafe_allow_html=True)
        
        # Logo con animación
        st.markdown(
            """
            <style>
            .sidebar-logo {
                display: flex;
                justify-content: center;
                width: 100%;
                margin-bottom: 1rem;
                animation: float 3s ease-in-out infinite;
            }
            .sidebar-logo img {
                max-width: 150px;
                height: auto;
                border-radius: 10px;
                box-shadow: 0 5px 15px rgba(79, 209, 199, 0.3);
                transition: transform 0.3s ease;
            }
            .sidebar-logo img:hover {
                transform: scale(1.05);
                box-shadow: 0 8px 25px rgba(79, 209, 199, 0.5);
            }
            @keyframes float {
                0%, 100% { transform: translateY(0px); }
                50% { transform: translateY(-10px); }
            }
            .fallback-logo {
                text-align: center;
                padding: 2rem;
                animation: float 3s ease-in-out infinite;
            }
            .fallback-logo .emoji {
                font-size: 4rem;
                margin-bottom: 1rem;
                display: block;
            }
            .fallback-logo h2 {
                color: #4FD1C7;
                margin: 0;
                text-shadow: 0 0 10px rgba(79, 209, 199, 0.5);
            }
            </style>
            <div class="sidebar-logo">
                <img src="https://raw.githubusercontent.com/JairAmado08/My-InChiral/main/imagenes1/inchiralucsur.png" alt="Inchiral Logo">
            </div>
            """,
            unsafe_allow_html=True
        )
        
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
        
        # Tabs mejorados - SOLO AGREGUÉ LA CUARTA PESTAÑA
        tab1, tab2, tab3, tab4 = st.tabs(["📋 Lista Completa", "💾 Descargar SMI", "🧪 Convertir a XYZ", "🌐 Visualizar 3D"])
        
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
            
            # NUEVA PESTAÑA 3D - CON ETIQUETAS AGREGADAS
            with tab4:
                st.markdown("""
                <div class="info-card">
                <h3 style='color: #2D3748; margin-top: 0;'>🌐 Visualización Molecular 3D</h3>
                """, unsafe_allow_html=True)
                
                # Selector de molécula
                col1, col2 = st.columns([2, 1])
                
                with col1:
                    selected_idx = st.selectbox(
                        "Selecciona un estereoisómero para visualizar:",
                        range(len(isomeros)),
                        format_func=lambda x: f"Isómero {x+1}: {isomeros[x]}"
                    )
                
                with col2:
                    st.markdown(f"""
                    <div style='text-align: center; padding: 1rem; background: rgba(79, 209, 199, 0.1); border-radius: 10px; margin-top: 1.5rem;'>
                    <strong>SMILES Seleccionado:</strong><br>
                    <code>{isomeros[selected_idx]}</code>
                    </div>
                    """, unsafe_allow_html=True)
                
                if st.button("🚀 Generar Visualización 3D", type="primary", key="viz_3d"):
                    with st.spinner("Generando estructura 3D..."):
                        xyz_content, mensaje = smiles_to_xyz(isomeros[selected_idx], selected_idx + 1)
                        
                        if xyz_content:
                            st.success("✅ Estructura 3D generada correctamente")
                            
                            # Parsear coordenadas XYZ
                            lines = xyz_content.strip().split('\n')
                            num_atoms = int(lines[0])
                            atoms_data = []
                            
                            for i in range(2, 2 + num_atoms):
                                parts = lines[i].split()
                                element = parts[0]
                                x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                                atoms_data.append([element, x, y, z])
                            
                            # Generar HTML con Three.js para visualización 3D
                            atoms_js = []
                            for i, (element, x, y, z) in enumerate(atoms_data):
                                color_map = {
                                    'C': 0x404040, 'H': 0xFFFFFF, 'O': 0xFF0000, 
                                    'N': 0x0000FF, 'S': 0xFFFF00, 'P': 0xFFA500,
                                    'F': 0x00FF00, 'Cl': 0x00FF00, 'Br': 0xA52A2A
                                }
                                size_map = {
                                    'H': 0.5, 'C': 0.7, 'N': 0.65, 'O': 0.6, 'S': 1.0, 'P': 1.1,
                                    'F': 0.5, 'Cl': 0.9, 'Br': 1.2
                                }
                                
                                color = color_map.get(element, 0x808080)
                                size = size_map.get(element, 0.6)
                                
                                atoms_js.append({
                                    'element': element,
                                    'x': x, 'y': y, 'z': z,
                                    'color': color,
                                    'size': size
                                })
                            
                            html_3d = f"""
                            <div style="width: 100%; height: 600px; border: 2px solid #4FD1C7; border-radius: 15px; background: linear-gradient(135deg, #2D3748, #4A5568);">
                                <div id="molecule-3d" style="width: 100%; height: 100%;"></div>
                            </div>
                            
                            <script src="https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js"></script>
                            <script>
                            const scene = new THREE.Scene();
                            const camera = new THREE.PerspectiveCamera(75, 800/600, 0.1, 1000);
                            const renderer = new THREE.WebGLRenderer({{antialias: true, alpha: true}});
                            
                            const container = document.getElementById('molecule-3d');
                            renderer.setSize(container.clientWidth, container.clientHeight);
                            renderer.setClearColor(0x000000, 0);
                            container.appendChild(renderer.domElement);
                            
                            const ambientLight = new THREE.AmbientLight(0x404040, 0.6);
                            scene.add(ambientLight);
                            const directionalLight = new THREE.DirectionalLight(0xffffff, 0.8);
                            directionalLight.position.set(1, 1, 1);
                            scene.add(directionalLight);
                            
                            const atoms = {str(atoms_js).replace("'", '"')};
                            
                            atoms.forEach((atom, index) => {{
                                const geometry = new THREE.SphereGeometry(atom.size, 32, 32);
                                const material = new THREE.MeshPhongMaterial({{color: atom.color}});
                                const sphere = new THREE.Mesh(geometry, material);
                                
                                sphere.position.set(atom.x, atom.y, atom.z);
                                scene.add(sphere);
                                
                                const canvas = document.createElement('canvas');
                                const context = canvas.getContext('2d');
                                canvas.width = 64;
                                canvas.height = 64;
                                context.font = '32px Arial';
                                context.fillStyle = 'white';
                                context.textAlign = 'center';
                                context.fillText(atom.element, 32, 40);
                                
                                const texture = new THREE.CanvasTexture(canvas);
                                const spriteMaterial = new THREE.SpriteMaterial({{map: texture}});
                                const sprite = new THREE.Sprite(spriteMaterial);
                                sprite.position.set(atom.x, atom.y + atom.size + 0.5, atom.z);
                                sprite.scale.set(1, 1, 1);
                                scene.add(sprite);
                            }});

                            // Etiqueta de la molécula
                            const molCanvas = document.createElement('canvas');
                            const molContext = molCanvas.getContext('2d');
                            molCanvas.width = 256;
                            molCanvas.height = 64;
                            molContext.font = '18px Arial';
                            molContext.fillStyle = '#4FD1C7';
                            molContext.textAlign = 'center';
                            molContext.fillText('Isómero {selected_idx + 1}', 128, 25);
                            molContext.fillText('{isomeros[selected_idx]}', 128, 45);
                            
                            const molTexture = new THREE.CanvasTexture(molCanvas);
                            const molSpriteMaterial = new THREE.SpriteMaterial({{map: molTexture}});
                            const molSprite = new THREE.Sprite(molSpriteMaterial);
                            
                            const avgX = atoms.reduce((sum, atom) => sum + atom.x, 0) / atoms.length;
                            const avgZ = atoms.reduce((sum, atom) => sum + atom.z, 0) / atoms.length;
                            molSprite.position.set(avgX, -8, avgZ);
                            molSprite.scale.set(4, 1, 1);
                            scene.add(molSprite);

                            camera.position.z = 15;
                            
                            let mouseX = 0, mouseY = 0;
                            let targetRotationX = 0, targetRotationY = 0;
                            let rotationX = 0, rotationY = 0;
                            let isMouseDown = false;
                            
                            container.addEventListener('mousedown', (event) => {{
                                isMouseDown = true;
                                mouseX = event.clientX;
                                mouseY = event.clientY;
                            }});
                            
                            container.addEventListener('mousemove', (event) => {{
                                if (isMouseDown) {{
                                    targetRotationY += (event.clientX - mouseX) * 0.01;
                                    targetRotationX += (event.clientY - mouseY) * 0.01;
                                    mouseX = event.clientX;
                                    mouseY = event.clientY;
                                }}
                            }});
                            
                            container.addEventListener('mouseup', () => {{
                                isMouseDown = false;
                            }});
                            
                            container.addEventListener('wheel', (event) => {{
                                event.preventDefault();
                                camera.position.z += event.deltaY * 0.01;
                                camera.position.z = Math.max(5, Math.min(50, camera.position.z));
                            }});
                            
                            function animate() {{
                                requestAnimationFrame(animate);
                                
                                rotationX += (targetRotationX - rotationX) * 0.1;
                                rotationY += (targetRotationY - rotationY) * 0.1;
                                
                                scene.rotation.x = rotationX;
                                scene.rotation.y = rotationY;
                                
                                renderer.render(scene, camera);
                            }}
                            
                            animate();
                            
                            window.addEventListener('resize', () => {{
                                camera.aspect = container.clientWidth / container.clientHeight;
                                camera.updateProjectionMatrix();
                                renderer.setSize(container.clientWidth, container.clientHeight);
                            }});
                            </script>
                            """
                            
                            st.components.v1.html(html_3d, height=650)
                            
                            # Información adicional
                            col1, col2, col3 = st.columns(3)
                            
                            with col1:
                                st.metric("Total de átomos", num_atoms)
                            
                            with col2:
                                elements = [atom[0] for atom in atoms_data]
                                unique_elements = list(set(elements))
                                st.metric("Elementos únicos", len(unique_elements))
                            
                            with col3:
                                h_count = elements.count('H')
                                st.metric("Átomos de H", h_count)
                            
                            # Mostrar composición
                            st.markdown("### 📈 Composición Atómica")
                            elements = [atom[0] for atom in atoms_data]
                            composition = {}
                            for elem in elements:
                                composition[elem] = composition.get(elem, 0) + 1
                            
                            composition_text = " | ".join([f"{elem}: {count}" for elem, count in sorted(composition.items())])
                            st.info(f"**Fórmula molecular:** {composition_text}")
                            
                            # Opción de descarga del XYZ individual
                            st.download_button(
                                label=f"📥 Descargar XYZ - Isómero {selected_idx + 1}",
                                data=xyz_content,
                                file_name=f"isomero_{selected_idx + 1}.xyz",
                                mime="text/plain"
                            )
                            
                        else:
                            st.error(f"❌ {mensaje}")
                
                st.markdown('</div>', unsafe_allow_html=True)
                
        else:
            with tab4:
                st.info("💡 Genera estereoisómeros primero para acceder a la visualización 3D")
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
