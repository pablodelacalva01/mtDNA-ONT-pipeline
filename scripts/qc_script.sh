#!/usr/bin/env bash
# # ---------------------------------------------------------------------------------------------------------------
# Autor:   Pablo de la Calva Castineira
# Fecha:   05/05/2026
# Titulo:  "Analisis bioinformatico y su implementacion en la practica clinica
#           diaria de datos de NGS de lecturas largas mediante tecnologia Nanopore"
# Version: 1.0
# Descripcion: Pipeline para analisis de variantes mitocondriales a partir
#              de lecturas largas Nanopore
#              Parte 1: Creacion de entorno de trabajo, directorios,
#              union de datos crudos, analisis de la calidad y recorte de las lecturas
# Entorno conda:
#   - ont_qc: porechop, nanoplot, nanofilt, filtlong
# Orden del procesamiento: 
#   - creacion de directorios → entrada de datos → nanoplot → porechop → nanofilt → filtlong → nanoplot 
# ---------------------------------------------------------------------------------------------------------------

set -euo pipefail

#0. CONFIGURACION DE ENTORNO: 
# Ubicacion del script y de la raiz del proyecto
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd -- "${SCRIPT_DIR}/.." && pwd)"

# Ubicacion de directorios
FASTQ_PASS="${PROJECT_DIR}/Data/Raw/fastq_pass"
OUTPUT="${PROJECT_DIR}/Data/Processed"
RESULTS="${PROJECT_DIR}/Results"
LOG_DIR="${OUTPUT}/logs"
TIMESTAMP=$(date +"%Y_%m_%d_%H:%M:%S")

# Creacion de directorios
mkdir -p "${OUTPUT}/00.Merged_Data" "${OUTPUT}/01.QC/nanoplot" \
    "${OUTPUT}/02.Filtered/porechop" "${OUTPUT}/02.Filtered/nanofilt" "${OUTPUT}/02.Filtered/filtlong" "${OUTPUT}/02.Filtered/nanoplot_filtered" \
    "${LOG_DIR}" \

#---------------------------------------------------------------------------------------------------------------------------------
#1. ENTRADA DE DATOS
echo "Iniciando pipeline de control de calidad y análisis de lecturas en Nanopore - [${TIMESTAMP}]"
echo "[Paso 1 - ${TIMESTAMP}] Identificación de muestras crudas (Merged por muestras)"
  
# Identificacion de archivos de entrada
BARCODES=($(ls -d ${FASTQ_PASS}/barcode* 2>/dev/null | xargs -n 1 basename))
if [[ ${#BARCODES[@]} -eq 0 ]]; then
    echo "ERROR: No se encontraron carpetas de barcodes en ${FASTQ_PASS}" >&2
    exit 1
fi
   
# Array: Almacen de archivos, 1 merged por muestra
MERGED_FILES=()

for BARCODE in "${BARCODES[@]}"; do
    BARCODE_DIR="${FASTQ_PASS}/${BARCODE}"
    MERGED_OUT="${OUTPUT}/00.Merged_Data/${BARCODE}_merged.fastq.gz"
    echo "Realizando merged: ${BARCODE}"

    # Comprobacion de archivos
    shopt -s nullglob
    FASTQ_FILES=( "${BARCODE_DIR}"/*.fastq.gz )
    shopt -u nullglob

    if [[ ${#FASTQ_FILES[@]} -eq 0 ]]; then
        echo "AVISO: ${BARCODE} no contiene archivos fastq.gz, se omite" >&2
        continue
    fi

    # Check de orden de los archivos de entrada
    mapfile -t FASTQ_FILES < <(printf '%s\n' "${FASTQ_FILES[@]}" | sort -V)
    cat "${FASTQ_FILES[@]}" > "${MERGED_OUT}"

    # Comprobaciones de archivos de salida
    gzip -t "${MERGED_OUT}"
    N_IN=$(zcat "${FASTQ_FILES[@]}" | awk 'END{print NR/4}')
    N_OUT=$(zcat "${MERGED_OUT}" | awk 'END{print NR/4}')
    if [[ "$N_IN" != "$N_OUT" ]]; then
        echo "ERROR: nº de lecturas distinto en ${BARCODE} (${N_IN} vs ${N_OUT})" >&2
        exit 1
    fi

    MERGED_FILES+=("${MERGED_OUT}")
    echo "Archivos correctamente juntados para la muestra ${BARCODE}, ubicado en: ${MERGED_OUT}"
done
echo "[Merge completado] Total de muestras: ${#MERGED_FILES[@]}"

#---------------------------------------------------------------------------------------------------------------------------------
#2. CONTROL DE CALIDAD
echo "[Paso 2 - ${TIMESTAMP}] Control de calidad inicial (NanoPlot)"

for MERGED in "${MERGED_FILES[@]}"; do
    SAMPLE=$(basename "${MERGED}" _merged.fastq.gz)
    # 2.1. NanoPlot: visualización de calidad y longitud de lecturas
    NanoPlot --fastq "${MERGED}" --outdir "${OUTPUT}/01.QC/nanoplot/${SAMPLE}" \
        --plots kde dot --N50 --threads 4 \
        2> >(grep -v "WARNING:root:" 2>/dev/null >&2)
     echo "Muestra ${SAMPLE} procesada QC con Nanoplot"
done
echo "[QC completado] Resultados en ${OUTPUT}/01.QC/"
#---------------------------------------------------------------------------------------------------------------------------------
#3. TRIMMING Y FILTRADO
echo "[Paso 3 - ${TIMESTAMP}] Trimming de adaptadores, filtrado de calidad y QC final (Porechop + NanoFilt + Filtlong + NanoPlot)"
# Parametros de filtrado: Ajustar umbrales
NANOFILT_HEADCROP=10   # Recorte de los primeros N bp
NANOFILT_Q=10          # Calidad minima media por read
NANOFILT_MINLEN=200    # Longitud minima (bp)
FILTLONG_MINLEN=200     # Longitud minima para Filtlong
FILTLONG_MINQ=80        # Calidad minima media para Filtlong (escala 0-100)

FILTERED_FILES=()
for MERGED in "${MERGED_FILES[@]}"; do
    SAMPLE=$(basename "${MERGED}" _merged.fastq.gz)
    echo "Filtrando muestra: ${SAMPLE}"
    PORECHOP_OUT="${OUTPUT}/02.Filtered/porechop/${SAMPLE}_porechop.fastq.gz"
    NANOFILT_OUT="${OUTPUT}/02.Filtered/nanofilt/${SAMPLE}_nanofilt.fastq.gz"
    FILTLONG_OUT="${OUTPUT}/02.Filtered/filtlong/${SAMPLE}_filtered.fastq.gz"

    # 3.1. Porechop: eliminacion de adaptadores ONT
    porechop --input "${MERGED}" --output "${PORECHOP_OUT}" \
        --threads 4 \
        2>&1 | tee "${LOG_DIR}/${SAMPLE}_porechop.log"

    # 3.2. NanoFilt:
    #   -q: Descarta reads con calidad media < Q10
    #   -l: Descarta reads < 500 bp
    #   --headcrop: Elimina los primeros 10 bp inestables de cada read
    gunzip -c "${PORECHOP_OUT}" \
        | NanoFilt -q "${NANOFILT_Q}" -l "${NANOFILT_MINLEN}" --headcrop "${NANOFILT_HEADCROP}" \
        | gzip > "${NANOFILT_OUT}" 
            
    # 3.3. Filtlong: optimización — selecciona las mejores lecturas por calidad
    #   --min_length: Descarta reads < 500 bp
    #   --min_mean_q: Descarta reads con calidad media < 80 (escala 0-100)
    filtlong --min_length "${FILTLONG_MINLEN}" --min_mean_q "${FILTLONG_MINQ}" "${NANOFILT_OUT}" \
        2> >(tee "${LOG_DIR}/${SAMPLE}_filtlong.log" >&2) | gzip > "${FILTLONG_OUT}"
    FILTERED_FILES+=("${FILTLONG_OUT}") 
    echo "Muestra ${SAMPLE} procesada con porechop y filtlong: ${FILTLONG_OUT}"

    # 3.4. Control de calidad post-filtrado (NanoPlot)
    SAMPLE=$(basename "${FILTLONG_OUT}" _filtered.fastq.gz)
    NanoPlot --fastq "${FILTLONG_OUT}" --outdir "${OUTPUT}/02.Filtered/nanoplot_filtered/${SAMPLE}" \
        --plots kde dot --N50 --threads 4 \
        2> >(grep -v "WARNING:root:" 2>/dev/null >&2)
done
echo "[Filtrado + Trimming + QC post-filtrado completado] Total muestras listas para alineamiento: ${#FILTERED_FILES[@]}"
#---------------------------------------------------------------------------------------------------------------------------------