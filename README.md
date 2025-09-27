# 🌀 Inchiral - Generador de Estereoisómeros

**Inchiral** es una aplicación interactiva y educativa desarrollada en **Streamlit** con soporte de **RDKit**, diseñada para detectar quiralidad en moléculas, generar estereoisómeros y exportar sus estructuras en distintos formatos.

## 🎯 Funcionalidades principales
- Detectar si una molécula es **quiral o aquiral** a partir de un código SMILES.  
- Identificar y mostrar los **centros quirales detectados**.  
- Generar automáticamente todos los **estereoisómeros posibles** (hasta 3 centros quirales).  
- Descargar los resultados en formato **.smi**.  
- Convertir las estructuras a **archivos XYZ** para visualización 3D.  
- Interfaz moderna con **CSS personalizado** y ejemplos de uso.  

## 🧪 ¿Qué es la quiralidad?
La **quiralidad** es una propiedad de ciertas moléculas cuya estructura no puede superponerse con su imagen especular, como sucede con las manos humanas.  

Esto da lugar a **enantiómeros**, que aunque químicamente similares, pueden tener efectos biológicos muy diferentes. Un caso histórico es el de la **talidomida**.

## 🖥️ Tecnologías utilizadas
- [Python](https://www.python.org/) 🐍  
- [Streamlit](https://streamlit.io/) 🌐  
- [RDKit](https://www.rdkit.org/) ⚛️  
- [itertools](https://docs.python.org/3/library/itertools.html) 🔄  
- [ZIP/IO](https://docs.python.org/3/library/zipfile.html) 💾  

## 🚀 Ejemplo de uso
1. Ingresar un código **SMILES** en la aplicación.  
2. Verificar si la molécula es **quiral o aquiral**.  
3. Generar y explorar los **estereoisómeros** posibles.  
4. Descargar resultados en `.smi` o convertirlos a `.xyz` para uso en programas de modelado 3D.  

### 📑 Posibles licencias para Inchiral

| **Licencia** | **Justificación** |
|--------------|--------------------|
| **MIT**      | ✅ Permite máxima libertad de uso y modificación, fomentando colaboración abierta. |
| **GPLv3**    | ✅ Garantiza que las modificaciones sigan siendo libres, promoviendo reproducibilidad en la comunidad científica. |
