#!/usr/bin/env bash
# --------------------------------------------------------------------------------------------
# Autor:   Pablo de la Calva Castineira
# Fecha:   05/05/2026
# Titulo:  "Analisis bioinformatico y su implementacion en la practica clinica
#           diaria de datos de NGS de lecturas largas mediante tecnologia Nanopore"
# Version: 1.0
# Descripcion: Pipeline para analisis de variantes mitocondriales a partir
#              de lecturas largas Nanopore
#              Parte 2: Alineamiento de secuencias pre-procesadas
# Entorno conda:
#   - ont_align: minimap2, ngmlr, gatk4, samtools, qualimap
# Orden del procesamiento:
#   indexado ref (samtools + minimap2) → BAM → sort → index 
#   minimap2/ngmlr → gatk4 → flagstat → qualimap
# --------------------------------------------------------------------------------------------

set -euo pipefail

#0. CONFIGURACION DE ENTORNO
# Activacion del entorno de Conda 
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate ont_align

# Variables de entorno
TFM_DIR="/home/pdelacalvac/TFM"
INPUT="${TFM_DIR}/Data/Processed/02.Filtered/filtlong"
OUTPUT="${TFM_DIR}/Data/Processed/03.Alignment"
RESULTS="${TFM_DIR}/Results"
REF="${TFM_DIR}/Data/Reference_Genome/mtDNA.fa"
LOG_DIR="${TFM_DIR}/Data/Processed/logs"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
THREADS=4
GATK4_MEM="3g" # Memoria para GATK4 (ajustar segun RAM disponible)
GATK4_TMP="${TFM_DIR}/Data/tmp"

# Creacion de directorios de salida
mkdir -p "${OUTPUT}/minimap2" "${OUTPUT}/ngmlr" "${OUTPUT}/QC" \
"${OUTPUT}/QC/minimap2" "${OUTPUT}/QC/ngmlr/" "${OUTPUT}/markdup/minimap2" "${OUTPUT}/markdup/ngmlr/" \
"${OUTPUT}/QC/markdup/minimap2" "${OUTPUT}/QC/markdup/ngmlr/" \
"${GATK4_TMP}"

LOGFILE="${LOG_DIR}/pipeline_align_${TIMESTAMP}.log"
exec > >(tee -a "${LOGFILE}") 2>&1

# ------------------------------------------------------------------------------
#1. IDENTIFICACION DE MUESTRAS
echo "Iniciando pipeline de alineamiento - [${TIMESTAMP}]"
FILTERED_FILES=("${INPUT}"/*_filtered.fastq.gz)
if [[ ${#FILTERED_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No se encontraron archivos filtrados en ${INPUT}" >&2
    exit 1
fi
echo "[Paso 1 - Identificación de muestras pre-procesadas] Muestras detectadas: ${#FILTERED_FILES[@]}"

# ------------------------------------------------------------------------------
#2. INDEXACION REFERENCIA
#2.1. INDEXACION DE LA REFERENCIA (SAMTOOLS)
echo "[Paso 2.1. - Indexación de la referencia] Indexando referencia (samtools): ${REF}"
if [[ ! -f "${REF}.fai" ]]; then
    samtools faidx "${REF}"
    echo "[Paso 2.1. - Indexación de la referencia] Índice samtools creado: ${REF}.fai"
else
    echo "[Paso 2.1. - Indexación de la referencia] Índice ya existe, se omite."
fi

#2.2. INDEXACIÓN DE LA REFERENCIA (MINIMAP2)
echo "[Paso 2.2. - Indexando referencia]: Minimap2"
MINIMAP2_INDEX="${REF%.fa}.mmi"
if [[ ! -f "${MINIMAP2_INDEX}" ]]; then
    minimap2 -x map-ont -d "${MINIMAP2_INDEX}" "${REF}"
    echo "[Paso 2.2. - Indexación de la referencia] Índice minimap2 creado: ${MINIMAP2_INDEX}"
else
    echo "[Paso 2.2. - Indexación de la referencia] Índice ya existe, se omite."
fi

#2.3. CREACION DEL DICCIONAR DE SECUENCIAS PARA GATK4 (MARKDUPLICATES)
REF_DICT="${REF%.fa}.dict"
echo "[Paso 2.3. - Diccionario de referencia para GATK4]"
if [[ ! -f "${REF_DICT}" ]]; then
    gatk --java-options "-Xmx${GATK4_MEM}" CreateSequenceDictionary \
        -R "${REF}" -O "${REF_DICT}"
    echo "[Paso 2.3. Diccionario creado]: ${REF_DICT}"
else
    echo "[Paso 2.3. Diccionario ya existe, se omite.]"
fi

#2.4. VALIDAR Y CONVERTIR SAM → BAM
sam_to_sorted_bam() {
    local SAM="$1"
    local BAM="$2"
    local SAMPLE="$3"
    local TOOL="$4"
    echo "[${TOOL} - ${SAMPLE}] Validando SAM..."
    if ! samtools quickcheck "${SAM}"; then
        echo "ERROR: SAM corrupto o truncado para ${SAMPLE} (${TOOL}): ${SAM}" >&2
        exit 1
    fi
    # Solución bug NGMLR 0.2.7 sólo para esta herramienta:
    if [[ "${TOOL}" == "ngmlr" ]]; then
        echo "[${TOOL} - ${SAMPLE}] Corrigiendo MAPQ inválidos (bug NGMLR 0.2.7)..."
        local SAM_FIXED="${SAM%.sam}_fixed.sam"
        awk 'BEGIN{OFS="\t"} /^@/{print; next} {if($5<0 || $5>255) $5=0; print}' \
            "${SAM}" > "${SAM_FIXED}"
        mv "${SAM_FIXED}" "${SAM}"
    fi

    local N_READS
    N_READS=$(samtools view -c "${SAM}")
    echo "[${TOOL} - ${SAMPLE}] Reads en SAM: ${N_READS}"
    if [[ "${N_READS}" -eq 0 ]]; then
        echo "ERROR: SAM sin reads para ${SAMPLE} (${TOOL})" >&2
        exit 1
    fi
    
    # Paso intermedio: SAM → BAM sin ordenar primero (evita bug de samtools sort) luego ordenar desde BAM
    local BAM_UNSORTED="${SAM%.sam}_unsorted.bam"
    echo "[${TOOL} - ${SAMPLE}] SAM → BAM sin ordenar..."
    samtools view -@ "${THREADS}" -bS "${SAM}" -o "${BAM_UNSORTED}"
    rm -f "${SAM}"

    echo "[${TOOL} - ${SAMPLE}] Ordenando BAM..."
    samtools sort -@ "${THREADS}" -o "${BAM}" "${BAM_UNSORTED}"
    rm -f "${BAM_UNSORTED}"

    samtools index "${BAM}"
    echo "[${TOOL} - ${SAMPLE}] BAM listo: ${BAM}"
}

# ------------------------------------------------------------------------------
# 3.1. ALINEAMIENTO MINIMAP2
echo "[Paso 3.1. - Alineamiento de las secuencias Minimap2]"
BAM_MINIMAP2=()

for FILTERED in "${FILTERED_FILES[@]}"; do
    SAMPLE=$(basename "${FILTERED}" _filtered.fastq.gz)
    MINIMAP2_DIR="${OUTPUT}/minimap2/${SAMPLE}"
    mkdir -p "${MINIMAP2_DIR}"

    SAM="${MINIMAP2_DIR}/${SAMPLE}_minimap2.sam"
    BAM="${MINIMAP2_DIR}/${SAMPLE}_minimap2.bam"
    #3.1. ALINEAMIENTO MINIMAP2
#   -ax map-ont: Parametro para lecturas largas ONT
#   -L: Evitar corrupcion BAM en UL reads por los CIGAR
#   --secondary=no: Descarta alineamientos secundarios 
#   -R: Read Group (para MarkDuplicates de GATK4)

    RG_TAG="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:mtDNA_lib\tPL:ONT\tPU:${SAMPLE}"
    minimap2 -ax map-ont -L -t "${THREADS}" --secondary=no -R "${RG_TAG}" "${MINIMAP2_INDEX}" "${FILTERED}" > "${SAM}" \
        2> "${LOG_DIR}/${SAMPLE}_minimap2.log"
    sam_to_sorted_bam "${SAM}" "${BAM}" "${SAMPLE}" "minimap2"
    BAM_MINIMAP2+=("${BAM}")

    echo "[Alineamiento Minimap2 - Completado]: ${SAMPLE} alineado"
done
echo "[Paso 3.1. - Alineamiento Minimap2 - Completado]"

# ------------------------------------------------------------------------------
# 3.2. ALINEAMIENTO NGMLR
echo "[Paso 3.2. - Alineamiento de las secuencias NGMLR]"
BAM_NGMLR=()

for FILTERED in "${FILTERED_FILES[@]}"; do
    SAMPLE=$(basename "${FILTERED}" _filtered.fastq.gz)
    NGMLR_DIR="${OUTPUT}/ngmlr/${SAMPLE}"
    mkdir -p "${NGMLR_DIR}"

    SAM="${NGMLR_DIR}/${SAMPLE}_ngmlr.sam"
    BAM="${NGMLR_DIR}/${SAMPLE}_ngmlr.bam"
    RG_TAG="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:mtDNA_lib\tPL:ONT\tPU:${SAMPLE}"

    echo "[NGMLR - ${SAMPLE}] Alineando..."
    # NGMLR usa la referencia directamente.
    # --rg-id, --rg-sm: flags de Read Group en NGMLR
    # -x ont: preset para Oxford Nanopore
    ngmlr -t "${THREADS}" -r "${REF}" -q "${FILTERED}" -x ont \
        --rg-id "${SAMPLE}" --rg-sm "${SAMPLE}" --rg-lb "mtDNA_lib" --rg-pl "ONT" \
        -o "${SAM}" \
        > "${LOG_DIR}/${SAMPLE}_ngmlr.log" 2>&1

    sam_to_sorted_bam "${SAM}" "${BAM}" "${SAMPLE}" "ngmlr"
    BAM_NGMLR+=("${BAM}")
    echo "[NGMLR - ${SAMPLE}] Completado"
done
echo "[Paso 3.2. - Alineamiento NGMLR - Completado]"

#---------------------------------------------------------------------------------------------------------------------------------
#4. DETECCION DE DUPLICADOS (GATK4 - MARKDUPLICATES)
echo "[Paso 4 - Detección de duplicados: GATK4 MarkDuplicates]"
BAM_MARKDUP_MINIMAP2=()
BAM_MARKDUP_NGMLR=()

run_markdup() {
    local BAM="$1"
    local SAMPLE="$2"
    local TOOL="$3"
    local MARKDUP_DIR="${OUTPUT}/markdup/${TOOL}/${SAMPLE}"
    mkdir -p "${MARKDUP_DIR}"

    local BAM_MD="${MARKDUP_DIR}/${SAMPLE}_${TOOL}_markdup.bam"
    local METRICS="${MARKDUP_DIR}/${SAMPLE}_${TOOL}_markdup_metrics.txt"
    echo "[MarkDuplicates - ${TOOL} - ${SAMPLE}] Procesando..." >&2
    
    gatk --java-options "-Xmx${GATK4_MEM} -Djava.io.tmpdir=${GATK4_TMP}" \
        MarkDuplicates \
        --INPUT "${BAM}" --OUTPUT "${BAM_MD}" \
        --METRICS_FILE "${METRICS}" \
        --ASSUME_SORT_ORDER coordinate \
        --READ_NAME_REGEX null \
        --TAGGING_POLICY All \
        --REMOVE_DUPLICATES false \
        --CREATE_INDEX true \
        --VALIDATION_STRINGENCY LENIENT \
        --TMP_DIR "${GATK4_TMP}" \
        > "${LOG_DIR}/${SAMPLE}_${TOOL}_markdup.log" 2>&1

    echo "[MarkDuplicates - ${TOOL} - ${SAMPLE}] Completado. Métricas: ${METRICS}" >&2
    echo "${BAM_MD}"
} 
for BAM in "${BAM_MINIMAP2[@]}"; do
    SAMPLE=$(basename "${BAM}" _minimap2.bam)
    BAM_MD=$(run_markdup "${BAM}" "${SAMPLE}" "minimap2")
    BAM_MARKDUP_MINIMAP2+=("${BAM_MD}")
done

for BAM in "${BAM_NGMLR[@]}"; do
    SAMPLE=$(basename "${BAM}" _ngmlr.bam)
    BAM_MD=$(run_markdup "${BAM}" "${SAMPLE}" "ngmlr")
    BAM_MARKDUP_NGMLR+=("${BAM_MD}")
done

echo "[Paso 4 - MarkDuplicates] Completado"
 
#---------------------------------------------------------------------------------------------------------------------------------
#5. INFORME DE ALINEAMIENTOS
#5.1. Samtools flagstat + Creacion informe global

echo "[Paso 5 - Control de calidad de los alineamientos]"
mkdir -p "${RESULTS}/Alignment"
REPORT_TXT="${RESULTS}/Alignment/alignment_report_${TIMESTAMP}.txt"
    #Cabecera del report
{
    echo "------------------------------------------------------------------------------"
    echo "INFORME GLOBAL DE CALIDAD"
    echo "Fecha: $(date)"
    echo " ------------------------------------------------------------------------------"
} > "${REPORT_TXT}"

run_qc_bam() {
    local BAM="$1"
    local SAMPLE="$2"
    local TOOL="$3"
    local QC_DIR="$4"

    mkdir -p "${QC_DIR}"
    {
        echo "------------------------------------------------------------"
        echo "SAMPLE: ${SAMPLE} | TOOL: ${TOOL}"
        echo "BAM: ${BAM}"
        echo "------------------------------------------------------------"
        # Conteo rapido de reads por categoria (mapped, paired, secondary...)
        echo "# FLAGSTAT"
        samtools flagstat "${BAM}"
        # Reads mapeados/no mapeados por referencia
        echo "# IDXSTATS"
        samtools idxstats "${BAM}"
        # Profundidad y cobertura por posicion
        echo "# COVERAGE"
        samtools coverage "${BAM}"
        # Estadisticas detalladas: longitud de lectura, GC, error rate, depth...
        echo "# STATS"
        samtools stats "${BAM}" | grep ^SN
    } >> "${REPORT_TXT}"

    #5.1. Qualimap bamqc
    echo "[Informe global de calidad - Qualimap]: ${SAMPLE}"
    # QC visual detallado
    # calidad de mapeo, sesgo GC 
    # --java-mem-size: ajustar segun RAM disponible
    qualimap bamqc -bam "${BAM}" -outdir "${QC_DIR}/qualimap_${SAMPLE}" \
        -outformat HTML --java-mem-size=2G -nt "${THREADS}" \
        > "${LOG_DIR}/${SAMPLE}_${TOOL}_qualimap.log" 2>&1
}

#5.1. QC para cada BAM de minimap2
echo "[Paso 5.1.] QC minimap2 (pre-markdup)"
for BAM in "${BAM_MINIMAP2[@]}"; do
    SAMPLE=$(basename "${BAM}" _minimap2.bam)
    run_qc_bam "${BAM}" "${SAMPLE}" "minimap2" "${OUTPUT}/QC/minimap2/${SAMPLE}"
done

echo "[Paso 5.2.] QC ngmlr (pre-markdup)"
for BAM in "${BAM_NGMLR[@]}"; do
    SAMPLE=$(basename "${BAM}" _ngmlr.bam)
    run_qc_bam "${BAM}" "${SAMPLE}" "ngmlr" "${OUTPUT}/QC/ngmlr/${SAMPLE}"
done

echo "[Paso 5.3.] QC minimap2 (post-markdup)"
for BAM in "${BAM_MARKDUP_MINIMAP2[@]}"; do
    SAMPLE=$(basename "${BAM}" _minimap2_markdup.bam)
    run_qc_bam "${BAM}" "${SAMPLE}" "minimap2_markdup" "${OUTPUT}/QC/markdup/minimap2/${SAMPLE}"
done

echo "[Paso 5.4.] QC ngmlr (post-markdup)"
for BAM in "${BAM_MARKDUP_NGMLR[@]}"; do
    SAMPLE=$(basename "${BAM}" _ngmlr_markdup.bam)
    run_qc_bam "${BAM}" "${SAMPLE}" "ngmlr_markdup" "${OUTPUT}/QC/markdup/ngmlr/${SAMPLE}"
done
echo "[Paso 5 completado] Reporte generado en: ${REPORT_TXT}"

echo "[Pipeline de alineamiento - Completado]: $(date)"

#--------------------------------------------------------------------------------------------------------------------------------