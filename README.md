# mtDNA-ONT Pipeline
Flujo de análisis bioinformático para lecturas largas de Nanopore aplicado al estudio mitocondrial en SARS-CoV-2

Desarrollado como parte del TFM:
"Análisis bioinformático y su implementación en la práctica clínica diaria de datos de NGS de lecturas largas mediante tecnología Nanopore"
- Pablo de la Calva Castiñeira, 2026

## Requisitos 
- Linux (Ubuntu >20.04)
- Conda (https://docs.conda.io/en/latest/miniconda.html)
- Java (versión =>11)
- Internet/Ethernet/WiFi (descarga y actualización de bases de datos)

## Instalación 
###
### 1. Clonar el repositorio
```bash
git clone https://github.com/pablodelacalva01/mtDNA-ONT-pipeline.git
cd mtDNA-ONT-pipeline

### 2. Crear entornos de conda
```bash
conda env create -f enviroment/ont_qc.yml
conda env create -f enviroment/ont_align.yml
conda env create -f enviroment/ont_vc.yml
conda env create -f enviroment/anot_annot.yml
```

### 3. Setup de bases de datos para las anotaciones (aplicar solo la primera vez)
````bash
conda activate ont_annot
bash scripts/04_annot_setup.sh
```
## Uso 

Edita las variables de entorno al inicio de cada script para adaptarlas a tu sistema
```bash
# Paso 1 - QC y filtrado
conda activate ont_qc
bash scripts/01_qc_script.sh

# Paso 2 - Alineamiento
conda activate ont_align
bash scripts/02_alignment_script.sh

# Paso 3 - Llamada de variantes
conda activate ont_vc
bash scripts/03_vc_script.sh

# Paso 4 - Anotación de variantes
conda activate ont_annot
bash scripts/04_annot_script.sh
```

## Resultados de ejemplo
En `results/example` se incluyen resultados de muestra para verificar que el pipeline funciona correctamente en el sistema.

## Herramientas utilizadas
- **QC**: Nanoplot, Porechop, NanoFilt y Filtlong
-**Alineamiento**: Minimap2, NGMLR, Samtools, GATK4 y Qualimap
-**Llamada de variantes**: Clair3, Mutserve2, Bcftools y Samtools
- **Anotación**: SnpEff, bcftools, HaploGrep3, ClinVar y gnomAD (mt)

## Licencia 
MIT License - ver [LICENSE](LICENSE)
