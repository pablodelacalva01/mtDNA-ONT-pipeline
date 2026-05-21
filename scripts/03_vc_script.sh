#!/usr/bin/env bash
# --------------------------------------------------------------------------------------------
# Autor:   Pablo de la Calva Castineira
# Fecha:   05/05/2026
# Titulo:  "Analisis bioinformatico y su implementacion en la practica clinica
#           diaria de datos de NGS de lecturas largas mediante tecnologia Nanopore"
# Version: 1.0
# Descripcion: Pipeline para analisis de variantes mitocondriales a partir
#              de lecturas largas Nanopore
#              Parte 3: Variant calling de variantes mitocondrial (ONT)
# Entorno conda:
#   - ont_vc: clair3, mutverse2, bcftools, samtools
# Callers:
#   - Clair3 → SNVs + indels (modelo r1041_e82_400bps_sup_v500)
#   - Mutserve2 → heteroplasmias mitocondriales
# Alineadores de entrada:
#   - minimap2 (BAMs post-MarkDuplicates)
#   - ngmlr    (BAMs post-MarkDuplicates)
# Matriz de comparacion (4 VCFs por muestra):
#   MC: minimap2 + Clair3   | MM: minimap2 + Mutserve2
#   NC: ngmlr    + Clair3   | NM: ngmlr    + Mutserve2
# Estrategia de comparacion:
#   - Merge de los 4 VCFs → union de todas las variantes detectadas
#   - Comparacion por pares (6 combinaciones) con indice de Jaccard
#   - Informe completo por muestra
# Orden del procesamiento:
#   BAM → Clair3 → VCF
#   BAM → Mutserve2 → VCF
#   VCFs → bcftools isec/merge → comparacion + informe
# --------------------------------------------------------------------------------------------

set -euo pipefail

#0. CONFIGURACION DE ENTORNO
# Activacion del entorno de Conda 
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate ont_vc


# Variables de entorno
TFM_DIR="/home/pdelacalvac/TFM"
INPUT_MARKDUP="${TFM_DIR}/Data/Processed/03.Alignment/markdup"
OUTPUT="${TFM_DIR}/Data/Processed/04.Variants"
RESULTS="${TFM_DIR}/Results"
REF="${TFM_DIR}/Data/Reference_Genome/mtDNA.fasta"
LOG_DIR="${TFM_DIR}/Data/Processed/logs"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
THREADS=4


# Requisitios para el VC (CLAIR3 + MUTSERVE)
CLAIR3_BIN="/home/pdelacalvac/miniconda3/envs/ont_vc/bin" #(modificar ruta absoluta segun usuario)
export PATH="${CLAIR3_BIN}:${PATH}"
#Forzar Java del entorno conda (openjdk=17) 
export JAVA_HOME=$CONDA_PREFIX
export PATH=$JAVA_HOME/bin:$PATH

# Modelo Clair3 — directorio que contiene pileup.pt y full_alignment.pt
CLAIR3_MODEL="${TFM_DIR}/Data/Processed/clair3_models/"

# Parametros Clair3
CLAIR3_PLATFORM="ont"
CLAIR3_MIN_SNP_AF="0.01"      # AF minimo para SNVs (variable)
CLAIR3_MIN_INDEL_AF="0.05"    # AF minimo para indels (variable)
 
# Parametros Mutserve2
MUTSERVE_MIN_HET="0.01"       # Frecuencia minima de heteroplasmia (variable)
MUTSERVE_MAPQ="20"            # MAPQ minimo de reads (variable)
MUTSERVE_BASEQ="15"           # Base quality minima (variable)
 
# Contig del mtDNA en la referencia (rCRS)
CONTIG="NC_012920.1"

# Alineadores a procesar
ALIGNERS=("minimap2" "ngmlr")

#CREACIÓN DE DIRECTORIOS DE SALIDA
mkdir -p "${OUTPUT}/clair3/minimap2" "${OUTPUT}/clair3/ngmlr" \
    "${OUTPUT}/mutserve2/minimap2" "${OUTPUT}/mutserve2/ngmlr" \
    "${OUTPUT}/comparison"
 
LOGFILE="${LOG_DIR}/pipeline_vc_${TIMESTAMP}.log"
exec > >(tee -a "${LOGFILE}") 2>&1

# FUNCIONES AUXILIARES 
#0.1. FX1: Indexar vcf.gz (forzado para evitar indices antiguos y/o corruptos)
index_vcf() {
    local VCF="$1"
    bcftools index -f --tbi "${VCF}"
}

#0.2. FX2: Filtracion de variantes
filter_normalize_vcf() {
    local VCF_IN="$1"
    local VCF_OUT="$2"
    local SAMPLE="$3"
    local TOOL="$4"
 
    echo "[${TOOL} - ${SAMPLE}] Filtrando (PASS) y normalizando VCF..."
    bcftools view --apply-filters PASS "${VCF_IN}" \
        | bcftools norm --fasta-ref "${REF}" --multiallelics -any \
            --output-type z --output "${VCF_OUT}"
    index_vcf "${VCF_OUT}"

    local N_VARS
    N_VARS=$(bcftools view -H "${VCF_OUT}" | wc -l)
    echo "[${TOOL} - ${SAMPLE}] Variantes PASS tras normalización: ${N_VARS}"
}

count_variants() {
    local VCF="$1"
    bcftools view -H "${VCF}" | wc -l
}

# --------------------------------------------------------------------------------------------
#1. IDENTIFICACION DE MUESTRAS
echo "Iniciando pipeline de variant calling - [${TIMESTAMP}]"
SAMPLES=()
for BAM in "${INPUT_MARKDUP}/minimap2/"/*/*_minimap2_markdup.bam; do
    SAMPLE=$(basename "${BAM}" _minimap2_markdup.bam)
    SAMPLES+=("${SAMPLE}")
done
 
if [[ ${#SAMPLES[@]} -eq 0 ]]; then
    echo "ERROR: No se encontraron BAMs post-MarkDuplicates en ${INPUT_MARKDUP}" >&2
    exit 1
fi
echo "[Paso 1 - Identificación de muestras BAM post-MarkDuplicates] Muestras detectadas: ${#SAMPLES[@]}"

# --------------------------------------------------------------------------------------------
#2. VERIFICAR ELEMENTOS PARA EL VC:
echo "[Paso 2 - Verificación de elementos para el variant calling]"

#2.1. Verificar modelo Clair3: Existe + archivos necesarios
echo "[Paso 2.1. - Verificación del modelo Clair3]"
if [[ ! -f "${CLAIR3_MODEL}/pileup.pt" || ! -f "${CLAIR3_MODEL}/full_alignment.pt" ]]; then
    echo "ERROR: Modelo Clair3 incompleto en ${CLAIR3_MODEL}" >&2
    exit 1
fi
echo "[Paso 2.1. - Modelo Clair3 verificado]: ${CLAIR3_MODEL}"

#2.2. Verificar indice de referencia para bcftools 
echo "[Paso 2.2. - Verificación de la referencia para bcftools]"
if [[ ! -f "${REF}.fai" ]]; then
    echo "[Paso 2.2. Indexando referencia para bcftools]"
    samtools faidx "${REF}"
fi
echo "[Paso 2 - Verificación completada]"

# --------------------------------------------------------------------------------------------
#3. VARIANT CALLING

#3.1. VC Clair3
#NOTA: Clair3 v2.0.1 tiene dos fases internas:
#   1. Pileup.pt
#   2. Full-alignment
# Ambos modelos deben estar en el mismo directorio (--model_path).
# El directorio de salida de Clair3 contiene:
#   - merge_output.vcf.gz      → VCF final con todas las variantes
#   - pileup.vcf.gz            → VCF fase pileup
#   - full_alignment.vcf.gz    → VCF fase full_alignment

echo "[Paso 3.1. Variant calling - Clair3]"
declare -A VCF_CLAIR3  # VCF_CLAIR3[aligner:sample]
 
for ALIGNER in "${ALIGNERS[@]}"; do
    echo "[Clair3 - Procesando muestras de ${ALIGNER}] "
    for SAMPLE in "${SAMPLES[@]}"; do
        BAM="${INPUT_MARKDUP}/${ALIGNER}/${SAMPLE}/${SAMPLE}_${ALIGNER}_markdup.bam"
        CLAIR3_OUT="${OUTPUT}/clair3/${ALIGNER}/${SAMPLE}"
        mkdir -p "${CLAIR3_OUT}"
 
        if [[ ! -f "${BAM}" ]]; then
            echo "ERROR: BAM no encontrado: ${BAM}" >&2
            exit 1
        fi
        echo "[Clair3 - ${ALIGNER} - ${SAMPLE}] Iniciando variant calling"
        # Parametros clair3:
        #   --ref:               Referencia FASTA
        #   --bam_fn:            BAM de entrada
        #   --output:            Directorio de salida
        #   --threads:           Hilos
        #   --platform:          "ont" para Oxford Nanopore
        #   --model_path:        Directorio con pileup.pt y full_alignment.pt
        #   --ctg_name:          Limitar al contig mtDNA
        #   --include_all_ctgs:  Desactivado (solo 1 contig especifico)
        #   --snp_min_af:        AF minimo para SNVs (0.01 = 1%, variable)
        #   --indel_min_af:      AF minimo para indels (variable)
        #   --enable_long_indel: Mejora deteccion de indels largos
        #   --no_phasing_for_fa: Desactiva phasing 
        #   --haploid_precise:   Modo "haploide"

        "${CLAIR3_BIN}/run_clair3.sh" \
            --ref_fn="${REF}" --bam_fn="${BAM}" --output="${CLAIR3_OUT}" \
            --threads="${THREADS}" \
            --platform="${CLAIR3_PLATFORM}" \
            --model_path="${CLAIR3_MODEL}" \
            --ctg_name="${CONTIG}" \
            --snp_min_af="${CLAIR3_MIN_SNP_AF}" --indel_min_af="${CLAIR3_MIN_INDEL_AF}" \
            --enable_long_indel \
            --no_phasing_for_fa \
            --haploid_precise \
            > "${LOG_DIR}/${SAMPLE}_${ALIGNER}_clair3.log" 2>&1
 
        # El VCF final de Clair3: merge_output.vcf.gz
        RAW_VCF="${CLAIR3_OUT}/merge_output.vcf.gz"
        if [[ ! -f "${RAW_VCF}" ]]; then
            echo "ERROR: Clair3 no generó merge_output.vcf.gz para ${SAMPLE} (${ALIGNER})" >&2
            exit 1
        fi

        # Filtrado de muestras validas + normalizar: VCF final
        FILTERED_VCF="${CLAIR3_OUT}/${SAMPLE}_${ALIGNER}_clair3_pass.vcf.gz"
        filter_normalize_vcf "${RAW_VCF}" "${FILTERED_VCF}" "${SAMPLE}" "clair3_${ALIGNER}"
        VCF_CLAIR3["${ALIGNER}:${SAMPLE}"]="${FILTERED_VCF}"
 
        echo "[Clair3 - ${ALIGNER} - ${SAMPLE}] Completado"
    done
done
echo "[Paso 3.1. Variant calling - Clair3]: Completado"

#--------------------------------------------------------------------------------------------
#3.2. VC Mutserve2
#NOTA: Mutserve2, especializado en heteroplasmia mitocondrial
# Trabaja sobre el BAM y produce:
#   - sample.txt: tabla de variantes con frecuencias de heteroplasmia
#   - sample.vcf: formato VCF estandar
echo "[Paso 3.2. Variant calling - Mutserve2]"
declare -A VCF_MUTSERVE  # VCF_MUTSERVE[aligner:sample]
 
for ALIGNER in "${ALIGNERS[@]}"; do
    echo "[Mutserve2 - Procesando muestras de ${ALIGNER}] "
    for SAMPLE in "${SAMPLES[@]}"; do
        BAM="${INPUT_MARKDUP}/${ALIGNER}/${SAMPLE}/${SAMPLE}_${ALIGNER}_markdup.bam"
        MUTSERVE_OUT="${OUTPUT}/mutserve2/${ALIGNER}/${SAMPLE}"
        mkdir -p "${MUTSERVE_OUT}"

    if [[ ! -f "${BAM}" ]]; then
            echo "ERROR: BAM no encontrado: ${BAM}" >&2
            exit 1
        fi
        echo "[Mutserve2 - ${ALIGNER} - ${SAMPLE}] Iniciando variant calling"
        # Parametros clave:
        #   call:               Subcomando de variant calling
        #   --output:           Prefijo de salida (genera .txt y .vcf)
        #   --reference:        Referencia mtDNA FASTA
        #   --threads:          Hilos
        #   --level:            Nivel de heteroplasmia minimo a reportar (0.01 = 1%, variable)
        #   --mapQ:             MAPQ minimo de reads (variable)
        #   --baseQ:            Base quality minima (variable)
        #Nota: El BAM de entrada o input no se declara explicitamente 
        mutserve2 call "${BAM}" \
            --output="${MUTSERVE_OUT}/${SAMPLE}_${ALIGNER}.vcf" --reference="${REF}" \
            --threads="${THREADS}" \
            --level="${MUTSERVE_MIN_HET}" --mapQ="${MUTSERVE_MAPQ}" --baseQ="${MUTSERVE_BASEQ}" \
            > "${LOG_DIR}/${SAMPLE}_${ALIGNER}_mutserve2.log" 2>&1
 
        # Mutserve2 genera el VCF
        RAW_VCF="${MUTSERVE_OUT}/${SAMPLE}_${ALIGNER}.vcf"
        if [[ ! -f "${RAW_VCF}" ]]; then
            echo "ERROR: Mutserve2 no generó VCF para ${SAMPLE} (${ALIGNER})" >&2
            exit 1
        fi
 
        # Comprimir, filtrar muestras validas + normalizar:
        bgzip -c "${RAW_VCF}" > "${RAW_VCF}.gz"
        index_vcf "${RAW_VCF}.gz"

        FILTERED_VCF="${MUTSERVE_OUT}/${SAMPLE}_${ALIGNER}_mutserve2_pass.vcf.gz"
        filter_normalize_vcf "${RAW_VCF}.gz" "${FILTERED_VCF}" "${SAMPLE}" "mutserve2_${ALIGNER}"
        VCF_MUTSERVE["${ALIGNER}:${SAMPLE}"]="${FILTERED_VCF}"
 
        echo "[Mutserve2 - ${ALIGNER} - ${SAMPLE}] Completado"
    done
done

echo "[Paso 3.2. Variant calling - Mutserve2]: Completado"
# --------------------------------------------------------------------------------------------
#4. COMPARACION VCFS
# Estrategia:
#   4.1. Merge de los 4 VCFs → union de TODAS las variantes detectadas por cualquier
#        combinacion caller/alineador. 
#   4.2. Comparacion por pares (6 combinaciones posibles) con bcftools isec.
#               ["MC_vs_MM"] → MC (minimap2+Clair3) vs MM (minimap2+Mutserve2)
#               ["MC_vs_NC"] → MC (minimap2+Clair3) vs NC (ngmlr+Clair3)
#               ["MC_vs_NM"] → MC (minimap2+Clair3) vs NM (ngmlr+Mutserve2)
#               ["MM_vs_NC"] → MM (minimap2+Mutserve2) vs NC(ngmlr+Clair3)
#               ["MM_vs_NM"] → MM (minimap2+Mutserve2) vs NM (ngmlr+Mutserve2)
#               ["NC_vs_NM"] → NC (ngmlr+Clair3) vs NM (ngmlr+Mutserve2)
#        Para cada par se calcula:
#          - N variantes en VCF1 y VCF2
#          - N variantes compartidas (presentes en ambos)
#          - N variantes privadas de cada uno
#          - Indice de Jaccard (concordancia entre callers)
#   4.3. Informe completo por muestra
# --------------------------------------------------------------------------------------------

echo "[Paso 4 Comparación de VCFs]"
mkdir -p "${RESULTS}/VariantCalling"
REPORT_VC="${RESULTS}/VariantCalling/variant_calling_report_${TIMESTAMP}.txt"
{
    echo "--------------------------------------------------------------------------------------------"
    echo "INFORME DE VARIANT CALLING - mtDNA ONT"
    echo "Fecha:            $(date)"
    echo "Referencia:       ${REF} (${CONTIG})"
    echo "Clair3 modelo:    ${CLAIR3_MODEL}"
    echo "--------------------------------------------------------------------------------------------"
} > "${REPORT_VC}"

for SAMPLE in "${SAMPLES[@]}"; do
    CMP_DIR="${OUTPUT}/comparison/${SAMPLE}"
    mkdir -p "${CMP_DIR}"
 
    VCF_MC="${VCF_CLAIR3[minimap2:${SAMPLE}]}"    # minimap2 + Clair3
    VCF_MM="${VCF_MUTSERVE[minimap2:${SAMPLE}]}"  # minimap2 + Mutserve2
    VCF_NC="${VCF_CLAIR3[ngmlr:${SAMPLE}]}"       # ngmlr    + Clair3
    VCF_NM="${VCF_MUTSERVE[ngmlr:${SAMPLE}]}"     # ngmlr    + Mutserve2
    
    # Verificar que todos los VCFs existen y están indexados
    for LABEL_VCF in "MC:${VCF_MC}" "MM:${VCF_MM}" "NC:${VCF_NC}" "NM:${VCF_NM}"; do
        LABEL="${LABEL_VCF%%:*}"
        VCF="${LABEL_VCF#*:}"
        if [[ -z "${VCF}" || ! -f "${VCF}" ]]; then
            echo "ERROR: VCF ${LABEL} no encontrado para ${SAMPLE}: '${VCF}'" >&2
            exit 1
        fi
        index_vcf "${VCF}"
    done

    # --------------------------------------------------------------------------------------------
    # 4.1. Merge de los 4 VCFs
    #   bcftools merge --merge none: no colapsa variantes, mantiene cada caller como columna separada. 
    #   Permite ver en qué callers aparece cada variante consultando las columnas de FORMAT/GT.
    # --------------------------------------------------------------------------------------------
    #4.1. Merge de los 4 VCFs
    echo "[Paso 4.1. - ${SAMPLE}] Generando merge de los 4 VCFs"
    MERGED_VCF="${CMP_DIR}/${SAMPLE}_all_callers_merged.vcf.gz"

    TMP_DIR="${CMP_DIR}/tmp_rename"
    mkdir -p "${TMP_DIR}"

    VCF_MC_R="${TMP_DIR}/MC.vcf.gz"
    VCF_MM_R="${TMP_DIR}/MM.vcf.gz"
    VCF_NC_R="${TMP_DIR}/NC.vcf.gz"
    VCF_NM_R="${TMP_DIR}/NM.vcf.gz"

    # --- Normalizar Clair3 VCFs: eliminar FORMAT/AF si existe, renombrar muestra a etiqueta ---
    for PAIR in "MC:${VCF_MC}:${VCF_MC_R}" "NC:${VCF_NC}:${VCF_NC_R}"; do
        LABEL="${PAIR%%:*}"
        REST="${PAIR#*:}"
        SRC="${REST%%:*}"
        DST="${REST#*:}"

        echo "${LABEL}" > "${TMP_DIR}/${LABEL}_name.txt"
        bcftools annotate -x FORMAT/AF --output-type u "${SRC}" \
            | bcftools reheader --samples "${TMP_DIR}/${LABEL}_name.txt" \
            | bcftools view --output-type z --output "${DST}"
        bcftools index -f --tbi "${DST}"
    done

    # --- Normalizar Mutserve2 VCFs: eliminar AF+BQ, recodificar GT, renombrar muestra ---
    normalize_mutserve() {
        local LABEL="$1"
        local SRC="$2"
        local DST="$3"

        echo "${LABEL}" > "${TMP_DIR}/${LABEL}_name.txt"

        # 1. Renombrar muestra
        # 2. Eliminar FORMAT/AF y FORMAT/BQ (incompatibles con Clair3)
        # 3. Recodificar GT: 0/1 → 0/1 se deja; 1/0 → 1 (homoplasmica/haploide)
        #    bcftools +fixploidy no existe en todas las versiones,
        #    usamos view + awk para reescribir el GT en el campo FORMAT
        bcftools reheader --samples "${TMP_DIR}/${LABEL}_name.txt" "${SRC}" \
            | bcftools annotate -x FORMAT/AF,FORMAT/BQ --output-type u \
            | bcftools view --output-type v \
            | awk '
                /^#/ { print; next }
                {
                    # Campo FORMAT es $9, datos de muestra es $10
                    n = split($9, fmt, ":")
                    split($10, val, ":")
                    # Buscar indice de GT
                    gt_idx = 0
                    for (i=1; i<=n; i++) if (fmt[i]=="GT") gt_idx=i
                    if (gt_idx > 0) {
                        g = val[gt_idx]
                        # 1/0 → 1  (homoplasmica detectada como haploide)
                        if (g == "1/0") val[gt_idx] = "1"
                        # 0/1 → 0/1 se mantiene (heteroplasmia real)
                        out = val[1]
                        for (i=2; i<=n; i++) out = out ":" val[i]
                        $10 = out
                    }
                    print
                }
            ' OFS='\t' \
            | bcftools view --output-type z --output "${DST}"
        bcftools index -f --tbi "${DST}"
    }

    normalize_mutserve "MM" "${VCF_MM}" "${VCF_MM_R}"
    normalize_mutserve "NM" "${VCF_NM}" "${VCF_NM_R}"

    # --- Merge de los 4 VCFs normalizados ---
    bcftools merge --output-type z --output "${MERGED_VCF}" \
        --merge none \
        --force-samples \
        "${VCF_MC_R}" "${VCF_MM_R}" "${VCF_NC_R}" "${VCF_NM_R}" \
        > "${LOG_DIR}/${SAMPLE}_merge.log" 2>&1

    index_vcf "${MERGED_VCF}"

    N_MERGED=$(count_variants "${MERGED_VCF}")
    N_SNV=$(bcftools view -H -v snps    "${MERGED_VCF}" 2>/dev/null | wc -l)
    N_INDEL=$(bcftools view -H -v indels "${MERGED_VCF}" 2>/dev/null | wc -l)
    echo "[Paso 4.1 - ${SAMPLE}] Merge completado: ${N_MERGED} variantes únicas."

    # -------------------------------------------------------------------------
    # 4.2. Comparacion por pares (6 combinaciones)
    #   0000.vcf.gz → privadas de V1 (solo en V1)
    #   0001.vcf.gz → privadas de V2 (solo en V2)
    #   0002.vcf.gz → compartidas vistas desde V1  (compartidas correctas)
    #   0003.vcf.gz → compartidas vistas desde V2
    # -------------------------------------------------------------------------
    echo "[Paso 4.2 - ${SAMPLE}] Comparaciones por pares"
 
    # Definir los 6 pares a comparar
    declare -A PARES=(
        ["MC_vs_MM"]="${VCF_MC}:${VCF_MM}"
        ["MC_vs_NC"]="${VCF_MC}:${VCF_NC}"
        ["MC_vs_NM"]="${VCF_MC}:${VCF_NM}"
        ["MM_vs_NC"]="${VCF_MM}:${VCF_NC}"
        ["MM_vs_NM"]="${VCF_MM}:${VCF_NM}"
        ["NC_vs_NM"]="${VCF_NC}:${VCF_NM}"
    )
 
    declare -A RESULTADOS_PARES  #para el informe final
 
    for PAR in "MC_vs_MM" "MC_vs_NC" "MC_vs_NM" "MM_vs_NC" "MM_vs_NM" "NC_vs_NM"; do
        IFS=':' read -r V1 V2 <<< "${PARES[$PAR]}"
        OUT_ISEC="${CMP_DIR}/isec_${PAR}"
        mkdir -p "${OUT_ISEC}"
 
        echo "[${SAMPLE}] Comparando ${PAR}"
        bcftools isec --output-type z \
            --prefix "${OUT_ISEC}" \
            "${V1}" "${V2}" \
            > "${LOG_DIR}/${SAMPLE}_isec_${PAR}.log" 2>&1
 
        # Indexar los archivos generados
        for F in "${OUT_ISEC}"/000*.vcf.gz; do
            [[ -f "$F" ]] && index_vcf "$F"
        done
 
        N1=$(count_variants "${V1}")
        N2=$(count_variants "${V2}")
        N_SHARED=$(count_variants "${OUT_ISEC}/0002.vcf.gz")
        N_PRIV1=$(count_variants "${OUT_ISEC}/0000.vcf.gz")
        N_PRIV2=$(count_variants "${OUT_ISEC}/0001.vcf.gz")
 
        # Indice de Jaccard = compartidas / (N1 + N2 - compartidas)
        DENOM=$((N1 + N2 - N_SHARED))
        JACCARD="0.0000"
        if [[ ${DENOM} -gt 0 ]]; then
            JACCARD=$(awk "BEGIN {printf \"%.4f\", ${N_SHARED}/${DENOM}}")
        fi
 
        RESULTADOS_PARES["${PAR}"]="${N1}:${N2}:${N_SHARED}:${N_PRIV1}:${N_PRIV2}:${JACCARD}"
        echo "[${SAMPLE} - ${PAR}] V1=${N1} V2=${N2} compartidas=${N_SHARED} Jaccard=${JACCARD}"
    done

    # -------------------------------------------------------------------------
    # 4.3. Informe por muestra
    # -------------------------------------------------------------------------
    N_MC=$(count_variants "${VCF_MC}")
    N_MM=$(count_variants "${VCF_MM}")
    N_NC=$(count_variants "${VCF_NC}")
    N_NM=$(count_variants "${VCF_NM}")
 
    {
        echo ""
        echo "-------------------------------------------------------------------------"
        echo "  MUESTRA: ${SAMPLE}"
        echo "-------------------------------------------------------------------------"
        echo ""
        echo "-------Variantes por caller/alineador -----------------------------------"
        printf "  %-30s %6s variantes\n" "minimap2 + Clair3  (MC):" "${N_MC}"
        printf "  %-30s %6s variantes\n" "minimap2 + Mutserve2 (MM):" "${N_MM}"
        printf "  %-30s %6s variantes\n" "ngmlr    + Clair3  (NC):" "${N_NC}"
        printf "  %-30s %6s variantes\n" "ngmlr    + Mutserve2 (NM):" "${N_NM}"
        echo ""
        echo "------- Merge de los 4 callers (union total) ----------------------------"
        printf "  %-30s %6s variantes\n" "Total variantes únicas:" "${N_MERGED}"
        printf "  %-30s %6s\n"           "  → SNVs:"               "${N_SNV}"
        printf "  %-30s %6s\n"           "  → Indels:"             "${N_INDEL}"
        printf "  %-30s %s\n"            "  → Archivo:"            "${MERGED_VCF}"
        echo ""
        echo "------- Comparación por pares -------------------------------------------"        
        printf "  %-12s %6s %6s %10s %8s %8s %8s\n" \
            "Par" "N_V1" "N_V2" "Compartidas" "Priv_V1" "Priv_V2" "Jaccard"
        printf "  %s\n" "$(printf '%.0s-' {1..72})"
 
        for PAR in "MC_vs_MM" "MC_vs_NC" "MC_vs_NM" "MM_vs_NC" "MM_vs_NM" "NC_vs_NM"; do
            IFS=':' read -r N1 N2 NS NP1 NP2 JI <<< "${RESULTADOS_PARES[$PAR]}"
            printf "  %-12s %6s %6s %10s %8s %8s %8s\n" \
                "${PAR}" "${N1}" "${N2}" "${NS}" "${NP1}" "${NP2}" "${JI}"
        done
 
        echo ""
        echo "------- Estadísticas detalladas (bcftools stats) ------------------------"    
        echo "  ## minimap2 + Clair3 (MC)"
        bcftools stats "${VCF_MC}" 2>/dev/null | grep ^SN | sed 's/^/  /'
        echo "  ## minimap2 + Mutserve2 (MM)"
        bcftools stats "${VCF_MM}" 2>/dev/null | grep ^SN | sed 's/^/  /'
        echo "  ## ngmlr + Clair3 (NC)"
        bcftools stats "${VCF_NC}" 2>/dev/null | grep ^SN | sed 's/^/  /'
        echo "  ## ngmlr + Mutserve2 (NM)"
        bcftools stats "${VCF_NM}" 2>/dev/null | grep ^SN | sed 's/^/  /'
        echo ""
        echo "------- Archivos generados ----------------------------------------------"    
        printf "  %-20s %s\n" "Merge total:"     "${MERGED_VCF}"
        for PAR in "MC_vs_MM" "MC_vs_NC" "MC_vs_NM" "MM_vs_NC" "MM_vs_NM" "NC_vs_NM"; do
            printf "  %-20s %s\n" "isec ${PAR}:" "${CMP_DIR}/isec_${PAR}/"
        done
 
    } >> "${REPORT_VC}"
 
    unset PARES RESULTADOS_PARES
    echo "[Paso 4 - ${SAMPLE}] Completado."
done
 
{
    echo ""
    echo "-------------------------------------------------------------------------"
    echo "  FIN DEL INFORME"
    echo "  Completado: $(date)"
    echo "-------------------------------------------------------------------------"
} >> "${REPORT_VC}"

echo "[Paso 4 - Informe generado]: ${REPORT_VC}"
echo "[Pipeline de variant calling - Completado]: $(date)"

# --------------------------------------------------------------------------------------------