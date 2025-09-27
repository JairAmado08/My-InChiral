# 🌀 Inchiral - Detección de Quiralidad

**Inchiral** es una aplicación interactiva y educativa desarrollada en **Streamlit** cuyo propósito es detectar la quiralidad en moléculas y analizar sus posibles estructuras.

## 🎯 Propósito
- Comprender el concepto de **quiralidad** en química.  
- Detectar si una molécula es **quiral o aquiral**.  
- Visualizar moléculas en **3D**.  
- Servir como **apoyo didáctico** en el aprendizaje de química orgánica.  

## 🧪 ¿Qué es la quiralidad?
La **quiralidad** es una propiedad geométrica de ciertas moléculas en las que su estructura no es superponible con su imagen especular, como sucede con las manos humanas (izquierda y derecha).  

Esto da lugar a moléculas llamadas **enantiómeros**, que pueden tener propiedades químicas similares pero efectos biológicos muy diferentes.  
Un ejemplo histórico es el caso de la **talidomida**.

## 📂 Estructura del proyecto
Inchiral/  
│── app.py              # Archivo principal de la aplicación Streamlit  
│── requirements.txt    # Dependencias necesarias  
│── README.md           # Documentación del proyecto  
│── data/               # Archivos de ejemplo (moléculas, estructuras)  
│── assets/             # Imágenes y recursos gráficos  

## ⚙️ Instalación y uso

1. Clona este repositorio:
   git clone https://github.com/tu_usuario/inchiral.git  
   cd inchiral  

2. Crea un entorno virtual e instala las dependencias:
   pip install -r requirements.txt  

3. Ejecuta la aplicación:
   streamlit run app.py  

## 🖥️ Tecnologías utilizadas
- [Python](https://www.python.org/) 🐍  
- [Streamlit](https://streamlit.io/) 🌐  
- [RDKit](https://www.rdkit.org/) ⚛️  
- [Py3Dmol](https://3dmol.csb.pitt.edu/) 🔬  

## 🚀 Ejemplo de uso
1. Cargar una molécula desde un archivo `.mol` o `.sdf`.  
2. Visualizarla en 3D.  
3. Detectar automáticamente si es **quiral** o **aquiral**.  
4. Explorar posibles **enantiómeros**.  

## 🤝 Contribuciones
¡Las contribuciones son bienvenidas!  
Si deseas mejorar el proyecto, abre un **issue** o envía un **pull request**.  

## 📜 Licencia
Este proyecto está bajo la licencia **MIT**.  
Consulta el archivo [LICENSE](LICENSE) para más información.
