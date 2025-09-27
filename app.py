import streamlit as st
import itertools
import zipfile
import io

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
                posiciones.append((i, True))
                i += 2
            else:
                posiciones.append((i, False))
                i += 1
        else:
            i += 1
    
    n = len(posiciones)
    
    if n == 0:
        st.warning("⚠️ El SMILES no tiene centros quirales especificados con @ o @@.")
        return [], n
    elif n > 3:
        st.error("❌ El SMILES tiene más de 3 centros quirales.")
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
        layout="wide"
    )
    
    # CSS personalizado simplificado
    st.markdown("""
    <style>
    .main {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        padding: 0;
    }
    .info-card {
        background: rgba(255, 255, 255, 0.95);
        padding: 1.5rem;
        border-radius: 15px;
        box-shadow: 0 10px 30px rgba(0, 0, 0, 0.1);
        margin: 1rem 0;
    }
    .stButton > button {
        background: linear-gradient(45deg, #4FD1C7, #63B3ED);
        color: white;
        border: none;
        border-radius: 25px;
        padding: 0.75rem 2rem;
        font-weight: 600;
    }
    .main-title {
        background: linear-gradient(45deg, #4FD1C7, #63B3ED);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        font-size: 3rem;
        font-weight: 800;
        text-align: center;
        margin-bottom: 2rem;
    }
    </style>
    """, unsafe_allow_html=True)
    
    # Título principal
    st.markdown('<h1 class="main-title">🧬 INCHIRAL</h1>', unsafe_allow_html=True)
    st.markdown('<p style="text-align: center; color: white; font-size: 1.2rem; margin-bottom: 2rem;">Generador Avanzado de Estereoisómeros</p>', unsafe_allow_html=True)
    
    # Sidebar
    with st.sidebar:
        st.markdown("## 📋 Instrucciones")
        st.write("1. Ingresa un código SMILES")
        st.write("2. Sistema detecta centros quirales")
        st.write("3. Genera estereoisómeros")
        st.write("4. Conversión opcional a XYZ")
        
        st.markdown("---")
        
        st.markdown("## 💡 Ejemplos")
        st.code("CCO")
        st.code("CC(O)C(N)C")
        st.code("C[C@H](O)[C@@H](N)C")
        st.code("N[C@@H](C)C(=O)O")
    
    # Entrada de datos
    st.markdown('<div class="info-card">', unsafe_allow_html=True)
    st.markdown("## 📝 Entrada de Datos")
    
    smiles_input = st.text_input(
        "Ingresa el código SMILES:",
        placeholder="Ejemplo: C[C@H](O)[C@@H](N)C",
        help="Introduce tu molécula en formato SMILES"
    )
    st.markdown('</div>', unsafe_allow_html=True)
    
    if smiles_input:
        # Análisis de quiralidad
        st.markdown('<div class="info-card">', unsafe_allow_html=True)
        st.markdown("## 🔍 Análisis de Quiralidad")
        
        es_quiral, mensaje_quiralidad, centros_detectados = detectar_quiralidad(smiles_input)
        centros_especificados, posiciones_at = analizar_centros_existentes(smiles_input)

        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("### 🔎 Análisis con RDKit")
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
                st.warning("⚠️ RDKit no disponible")
        
        with col2:
            st.markdown("### 📋 Centros Especificados")
            if centros_especificados > 0:
                st.success(f"✅ {centros_especificados} centros con @ o @@")
                for pos in posiciones_at:
                    st.write(f"• Posición {pos}")
            else:
                st.warning("⚠️ No hay centros especificados")
        
        st.markdown('</div>', unsafe_allow_html=True)
        
        # Generación de estereoisómeros
        isomeros = []
        if centros_especificados > 0:
            with st.spinner("🔄 Generando estereoisómeros..."):
                isomeros, n_centros = generar_estereoisomeros(smiles_input)
        
        # Tabs
        if isomeros:
            tab1, tab2, tab3, tab4 = st.tabs(["📋 Lista", "💾 Descargar SMI", "🧪 Convertir XYZ", "🌐 Visualizar 3D"])
            
            with tab1:
                st.markdown('<div class="info-card">', unsafe_allow_html=True)
                st.markdown("## 🧪 Estereoisómeros Generados")
                
                col1, col2 = st.columns(2)
                for i, isomero in enumerate(isomeros):
                    if i % 2 == 0:
                        col1.code(f"{i+1}. {isomero}")
                    else:
                        col2.code(f"{i+1}. {isomero}")
                
                st.markdown('</div>', unsafe_allow_html=True)
            
            with tab2:
                st.markdown('<div class="info-card">', unsafe_allow_html=True)
                st.markdown("## 💾 Descarga SMI")
                
                smi_content = "\n".join(isomeros)
                st.download_button(
                    label="📥 Descargar archivo.smi",
                    data=smi_content,
                    file_name="estereoisomeros.smi",
                    mime="text/plain"
                )
                
                with st.expander("Vista previa"):
                    st.code(smi_content)
                    
                st.markdown('</div>', unsafe_allow_html=True)
            
            with tab3:
                st.markdown('<div class="info-card">', unsafe_allow_html=True)
                st.markdown("## 🧪 Conversión a XYZ")
                
                if st.button("🚀 Convertir todos a XYZ", type="primary"):
                    progress_bar = st.progress(0)
                    archivos_xyz = {}
                    mensajes_log = []
                    
                    for i, smiles in enumerate(isomeros):
                        progress = (i + 1) / len(isomeros)
                        progress_bar.progress(progress)
                        
                        xyz_content, mensaje = smiles_to_xyz(smiles, i+1)
                        mensajes_log.append(mensaje)
                        
                        if xyz_content:
                            archivos_xyz[f"mol_{i+1}.xyz"] = xyz_content
                    
                    with st.expander("Log de procesamiento"):
                        for mensaje in mensajes_log:
                            if "❌" in mensaje or "⚠️" in mensaje:
                                st.error(mensaje)
                            else:
                                st.success(mensaje)
                    
                    if archivos_xyz:
                        zip_data = crear_archivo_zip(archivos_xyz)
                        st.download_button(
                            label="📦 Descargar XYZ (ZIP)",
                            data=zip_data,
                            file_name="estereoisomeros_xyz.zip",
                            mime="application/zip"
                        )
                
                st.markdown('</div>', unsafe_allow_html=True)
            
            with tab4:
                st.markdown('<div class="info-card">', unsafe_allow_html=True)
                st.markdown("## 🌐 Visualización 3D con Etiquetas")
                
                selected_idx = st.selectbox(
                    "Selecciona un estereoisómero:",
                    range(len(isomeros)),
                    format_func=lambda x: f"Isómero {x+1}: {isomeros[x]}"
                )
                
                st.info(f"**SMILES:** `{isomeros[selected_idx]}`")
                
                if st.button("🚀 Generar Visualización 3D", type="primary"):
                    with st.spinner("Generando estructura 3D..."):
                        xyz_content, mensaje = smiles_to_xyz(isomeros[selected_idx], selected_idx + 1)
                        
                        if xyz_content:
                            st.success("✅ Estructura 3D generada")
                            
                            # Parsear XYZ
                            lines = xyz_content.strip().split('\n')
                            num_atoms = int(lines[0])
                            atoms_data = []
                            
                            for i in range(2, 2 + num_atoms):
                                parts = lines[i].split()
                                element = parts[0]
                                x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                                atoms_data.append([element, x, y, z])
                            
                            # Crear datos para JavaScript
                            atoms_js = []
                            for element, x, y, z in atoms_data:
                                color_map = {
                                    'C': 0x404040, 'H': 0xFFFFFF, 'O': 0xFF0000, 
                                    'N': 0x0000FF, 'S': 0xFFFF00, 'P': 0xFFA500,
                                    'F': 0x00FF00, 'Cl': 0x00FF00, 'Br': 0xA52A2A
                                }
                                size_map = {
                                    'H': 0.5, 'C': 0.7, 'N': 0.65, 'O': 0.6, 'S': 1.0, 'P': 1.1
                                }
                                
                                atoms_js.append({
                                    'element': element,
                                    'x': x, 'y': y, 'z': z,
                                    'color': color_map.get(element, 0x808080),
                                    'size': size_map.get(element, 0.6)
                                })
                            
                            # HTML con Three.js
                            html_3d = f"""
                            <div style="width: 100%; height: 600px; border: 2px solid #4FD1C7; border-radius: 15px; background: linear-gradient(135deg, #2D3748, #4A5568);">
                                <div id="molecule-3d" style="width: 100%; height: 100%;"></div>
                            </div>
                            
                            <script src="https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js"></script>
                            <script>
                            const scene = new THREE.Scene();
                            const camera = new THREE.PerspectiveCamera(75, 1, 0.1, 1000);
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
                            const moleculeGroup = new THREE.Group();
                            
                            // Calcular centro y límites
                            let centerX = 0, centerY = 0, centerZ = 0, minY = Infinity;
                            atoms.forEach(atom => {{
                                centerX += atom.x;
                                centerY += atom.y;
                                centerZ += atom.z;
                                minY = Math.min(minY, atom.y);
                            }});
                            centerX /= atoms.length;
                            centerY /= atoms.length;
                            centerZ /= atoms.length;
                            
                            // Crear átomos
                            atoms.forEach(atom => {{
                                const geometry = new THREE.SphereGeometry(atom.size, 32, 32);
                                const material = new THREE.MeshPhongMaterial({{color: atom.color}});
                                const sphere = new THREE.Mesh(geometry, material);
                                sphere.position.set(atom.x, atom.y, atom.z);
                                moleculeGroup.add(sphere);
                                
                                // Etiqueta del átomo
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
                                moleculeGroup.add(sprite);
                            }});
                            
                            // Etiqueta de la molécula
                            const molCanvas = document.createElement('canvas');
                            const molContext = molCanvas.getContext('2d');
                            molCanvas.width = 400;
                            molCanvas.height = 80;
                            
                            molContext.fillStyle = 'rgba(45, 55, 72, 0.9)';
                            molContext.fillRect(0, 0, 400, 80);
                            molContext.strokeStyle = '#4FD1C7';
                            molContext.lineWidth = 2;
                            molContext.strokeRect(2, 2, 396, 76);
                            
                            molContext.font = 'bold 16px Arial';
                            molContext.fillStyle = '#4FD1C7';
                            molContext.textAlign = 'center';
                            molContext.fillText('Estereoisómero {selected_idx + 1}', 200, 25);
                            
                            molContext.font = '12px monospace';
                            molContext.fillStyle = '#FFFFFF';
                            molContext.fillText('{isomeros[selected_idx]}', 200, 45);
                            
                            molContext.font = '10px Arial';
                            molContext.fillStyle = '#63B3ED';
                            molContext.fillText('Átomos: {num_atoms}', 200, 65);
                            
                            const molTexture = new THREE.CanvasTexture(molCanvas);
                            const molSpriteMaterial = new THREE.SpriteMaterial({{map: molTexture}});
                            const molSprite = new THREE.Sprite(molSpriteMaterial);
                            molSprite.position.set(centerX, minY - 3, centerZ);
                            molSprite.scale.set(6, 1.2, 1);
                            moleculeGroup.add(molSprite);
                            
                            scene.add(moleculeGroup);
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
                                camera.position.z += event.deltaY * 0.02;
                                camera.position.z = Math.max(5, Math.min(50, camera.position.z));
                            }});
                            
                            function animate() {{
                                requestAnimationFrame(animate);
                                
                                rotationX += (targetRotationX - rotationX) * 0.1;
                                rotationY += (targetRotationY - rotationY) * 0.1;
                                
                                moleculeGroup.rotation.x = rotationX;
                                moleculeGroup.rotation.y = rotationY;
                                
                                if (!isMouseDown) {{
                                    targetRotationY += 0.005;
                                }}
                                
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
                                st.metric("Total átomos", num_atoms)
                            with col2:
                                elements = [atom[0] for atom in atoms_data]
                                st.metric("Elementos únicos", len(set(elements)))
                            with col3:
                                h_count = elements.count('H')
                                st.metric("Átomos de H", h_count)
                            
                            # Composición
                            composition = {}
                            for elem in elements:
                                composition[elem] = composition.get(elem, 0) + 1
                            
                            composition_text = " | ".join([f"{elem}: {count}" for elem, count in sorted(composition.items())])
                            st.info(f"**Composición:** {composition_text}")
                            
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
            st.info("💡 Ingresa un SMILES con centros quirales especificados (@ o @@)")
    
    # Footer
    st.markdown("---")
    st.markdown("""
    <div style='text-align: center; color: white; padding: 2rem;'>
        <h3>🧬 INCHIRAL</h3>
        <p>Universidad Científica del Sur</p>
        <p><em>Generador Avanzado de Estereoisómeros</em></p>
    </div>
    """, unsafe_allow_html=True)

if __name__ == "__main__":
    main()
