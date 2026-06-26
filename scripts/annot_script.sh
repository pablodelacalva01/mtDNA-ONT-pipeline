#!/usr/bin/env bash
# --------------------------------------------------------------------------------------------
# Autor:   Pablo de la Calva Castineira
# Fecha:   06/05/2026
# Titulo:  "Analisis bioinformatico y su implementacion en la practica clinica
#           diaria de datos de NGS de lecturas largas mediante tecnologia Nanopore"
# Version: 1.2
# Descripcion: Pipeline de anotacion de variantes mitocondriales de alta confianza (ONT)
#              Parte 4: Anotacion de variantes mtDNA
#
# Input:   VCF consenso por muestra (Clair3 minimap2 + Clair3 ngmlr basado en Jaccard)
#          Nota: Generado en Parte 3: isec_MC_vs_NC/0002.vcf.gz
#
# Pipeline:
#   1. Identificacion de muestras y extraccion del VCF consenso MC_vs_NC
#   2. Filtrado AF >= 0.01 (variable)
#   3. Anotacion funcional con SnpEff (gen, region, impacto, HGVS)
#   4. Anotacion con ClinVar (CLNSIG, CLNDN)
#   5. Anotacion con gnomAD mt (AF_hom, AF_het, AC_hom, AC_het)
#   6. Clasificacion de haplogrupo con HaploGrep3
#   7. Exportacion a TSV anotado por muestra
#   8. Generacion del informe HTML clinico (via script Python)
#   9. Informe de texto resumen por muestra
#
# Entorno conda: ont_annot
# Prerequisitos: parte4_setup.sh ejecutado (HaploGrep3 + DBs instalados)
# --------------------------------------------------------------------------------------------

set -euo pipefail

# --------------------------------------------------------------------------------------------
# 0. CONFIGURACION
# Activacion del entorno de Conda 
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate ont_annot

# Variables de entorno
TFM_DIR="/home/pdelacalvac/TFM"
VARIANTS_DIR="${TFM_DIR}/Data/Processed/04.Variants"
OUTPUT="${TFM_DIR}/Data/Processed/05.Annotation"
RESULTS="${TFM_DIR}/Results"
DB_DIR="${TFM_DIR}/Data/Annotation_DBs"
LOG_DIR="${TFM_DIR}/Data/Processed/logs"
DB_VERSION_LOG="${DB_DIR}/db_versions.log"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
THREADS=4

# Referencia
CONTIG="NC_012920.1"

# Bases de datos
CLINVAR_DB="${DB_DIR}/clinvar_mt_rCRS.vcf.gz"
GNOMAD_DB="${DB_DIR}/gnomad_mt_rCRS.vcf.gz"

# HaploGrep3
HAPLOGREP_JAR="${CONDA_PREFIX}/bin/haplogrep3.jar"

# SnpEff
SNPEFF_DB="GRCh38.99"                                                        
SNPEFF_DIR="${CONDA_PREFIX}/share/snpeff-5.2-3"                              
SNPEFF_CFG="${SNPEFF_DIR}/snpEff.config"
SNPEFF_JAR="${SNPEFF_DIR}/snpEff.jar"
SNPEFF_DATA_DIR="${SNPEFF_DIR}/data"    

# Filtro AF (variable)
AF_MIN="0.01"

# Campos gnomAD a extraer (VCF de gnomAD mt)
# AF_hom: frecuencia homoplasmica | AF_het: frecuencia heteroplasmica
# AC_hom: conteo homoplasias      | AC_het: conteo heteroplasmias
GNOMAD_FIELDS="AF_hom,AF_het,AC_hom,AC_het"

#0.1. FUNCIONES AUXILIARES 

log()  { echo "[$(date '+%H:%M:%S')] $*"; }
err()  { echo "[ERROR] $*" >&2; exit 1; }

require() {
    command -v "$1" &>/dev/null \
        || err "Comando no encontrado: '$1'. Esta el entorno ont_annot activo?"
}

index_vcf() {
    bcftools index -f --tbi "$1"
}

count_variants() {
    bcftools view -H "$1" | wc -l
}

# --------------------------------------------------------------------------------------------
#1. VERIFICACIONES
echo "Iniciando pipeline de anotacion del variant calling - [${TIMESTAMP}]"
echo "[Paso 1 - Verificación de elementos para la anotación]"
# Cabecera del log
log "Iniciando pipeline de anotacion del variant calling - [${TIMESTAMP}]"
log "[Paso 1: Verificación de elementos para la anotación]"

require bcftools
require snpEff
require java
require bgzip

# Verificar bases de datos
for DB_FILE in "${CLINVAR_DB}" "${GNOMAD_DB}"; do
    [[ -f "${DB_FILE}" ]] \
        || err "Base de datos no encontrada: ${DB_FILE}. Ejecuta primero annot_setup.sh"
    [[ -f "${DB_FILE}.tbi" ]] \
        || err "Indice .tbi no encontrado para: ${DB_FILE}. Re-ejecuta annot_setup.sh"
done

# Verificar HaploGrep3
[[ -f "${HAPLOGREP_JAR}" ]] \
    || err "HaploGrep3 no encontrado en ${HAPLOGREP_JAR}. Ejecuta parte4_setup.sh"

# Verificar SnpEff DB
if [[ ! -d ${SNPEFF_DATA_DIR}/${SNPEFF_DB} ]]; then
    err "Base de datos SnpEff '${SNPEFF_DB}' no encontrada. Ejecuta: snpEff download ${SNPEFF_DB}"
fi
log "Base de datos SnpEff verificada: ${SNPEFF_DB}"

# --------------------------------------------------------------------------------------------
# 2. IDENTIFICACION DE MUESTRAS
echo "[Paso 2 - Verificación de elementos para la anotación]"
log "[Paso 2 - Verificación de elementos para la anotación]"

SAMPLES=()
for ISEC_DIR in "${VARIANTS_DIR}/comparison/"*/; do
    SAMPLE=$(basename "${ISEC_DIR}")
    CONSENSUS_VCF="${ISEC_DIR}/isec_MC_vs_NC/0002.vcf.gz"
    if [[ -f "${CONSENSUS_VCF}" ]]; then
        SAMPLES+=("${SAMPLE}")
    else
        log "AVISO: No se encontro VCF consenso MC_vs_NC para ${SAMPLE}, omitiendo."
    fi
done

if [[ ${#SAMPLES[@]} -eq 0 ]]; then
    err "No se encontraron muestras con VCF consenso MC_vs_NC en ${VARIANTS_DIR}/comparison/"
fi

log "Muestras detectadas: ${#SAMPLES[@]}"

# Crear directorios de salida
mkdir -p "${RESULTS}/Variants_Annotated" "${RESULTS}/Haplogroups" "${RESULTS}/TSV" "${RESULTS}/HTML" "${RESULTS}/AnnotationReports"

# --------------------------------------------------------------------------------------------
# 3. ANOTACION POR MUESTRA
echo "[Paso 3 - Anotación para cada muestra + generación de informe]"
log "[Paso 3 - Anotación para cada muestra + generación de informe]"

REPORT_ANNOT="${RESULTS}/AnnotationReports/annotation_report_${TIMESTAMP}.txt"
{
    echo "--------------------------------------------------------------------------------------------"
    echo "INFORME DE ANOTACION - mtDNA ONT"
    echo "Fecha:          $(date)"
    echo "Referencia:     ${CONTIG} (rCRS)"
    echo "--------------------------------------------------------------------------------------------"
    echo "SnpEff DB:      ${SNPEFF_DB}"
    echo "ClinVar:        ${CLINVAR_DB}"
    echo "gnomAD:         ${GNOMAD_DB}"
    echo "AF minimo:      ${AF_MIN}"
    echo "Versiones de bases de datos utilizadas:"
    # Extraer solo el bloque de versiones mas reciente del log (evita duplicados entre ejecuciones)
    tac "${DB_VERSION_LOG}" \
        | awk '/^$/{if(found) exit} /\t/{found=1; print}' \
        | tac
    echo ""
} > "${REPORT_ANNOT}"

for SAMPLE in "${SAMPLES[@]}"; do
    log "--------------------------------------------------------------------------------------------"
    log "Procesando muestra: ${SAMPLE}"
    log "--------------------------------------------------------------------------------------------"

    SAMPLE_OUT="${RESULTS}/Variants_Annotated/${SAMPLE}"
    mkdir -p "${SAMPLE_OUT}"

    CONSENSUS_VCF="${VARIANTS_DIR}/comparison/${SAMPLE}/isec_MC_vs_NC/0002.vcf.gz"

    # -------------------------------------------------------------------------
    # 3.1. Filtrado AF >= 0.01
    # El VCF consenso MC_vs_NC viene de Clair3 (formato GT:GQ:DP:AD)
    # AF se calcula con +fill-tags desde el campo AD (ref,alt)
    # bcftools view con expresion sobre [muestra:subcampo]
    # -------------------------------------------------------------------------
    echo "Paso 3.1 - Filtrado AF >= ${AF_MIN}: [${SAMPLE}]"
    log "Paso 3.1 - Filtrado AF >= ${AF_MIN}: [${SAMPLE}]"
    VCF_FILTERED="${SAMPLE_OUT}/${SAMPLE}_consensus_af_filtered.vcf.gz"

    bcftools view "${CONSENSUS_VCF}" \
        | bcftools +fill-tags -- -t AF \
        | bcftools view -i "FORMAT/AF[0:0] >= ${AF_MIN}" \
        -Oz -o "${VCF_FILTERED}"
    index_vcf "${VCF_FILTERED}"

    N_FILT=$(count_variants "${VCF_FILTERED}")
    N_ORIG=$(count_variants "${CONSENSUS_VCF}")
    log "[${SAMPLE}] Variantes antes del filtro AF: ${N_ORIG} → tras filtro: ${N_FILT}"

    # -------------------------------------------------------------------------
    # 3.2. Anotacion funcional con SnpEff
    # -noStats: sin HTML de estadisticas
    # -noLog:   sin telemetria a servidores de SnpEff
    # -------------------------------------------------------------------------
    echo "Paso 3.2 - Anotación funcional SnpEff: [${SAMPLE}]"
    log "Paso 3.2 - Anotación funcional SnpEff: [${SAMPLE}]"

    VCF_SNPEFF="${SAMPLE_OUT}/${SAMPLE}_snpeff.vcf.gz"
    VCF_FOR_SNPEFF="${SAMPLE_OUT}/${SAMPLE}_MT.vcf.gz"
    SNPEFF_CHRM="${SAMPLE_OUT}/rCRS_to_MT.txt"
    SNPEFF_RCRS="${SAMPLE_OUT}/MT_to_rCRS.txt"
    SNPEFF_LOG="${LOG_DIR}/${SAMPLE}_snpeff_${TIMESTAMP}.log"

    echo "NC_012920.1 MT" > "${SNPEFF_CHRM}"
    echo "MT NC_012920.1" > "${SNPEFF_RCRS}"

    # Asegura el uso de codones mitocondriales
    CODON_LINE="GRCh38.99.MT.codonTable : Vertebrate_Mitochondrial"
    grep -qF "${CODON_LINE}" "${SNPEFF_CFG}" \
        || echo "${CODON_LINE}" >> "${SNPEFF_CFG}"

    # Renombrar NC_012920.1 → MT para SnpEff
    bcftools annotate --rename-chrs "${SNPEFF_CHRM}" \
        "${VCF_FILTERED}" -Oz -o "${VCF_FOR_SNPEFF}"
    index_vcf "${VCF_FOR_SNPEFF}"

    # Anotar con SnpEff y restaurar NC_012920.1
    java -Xmx8g -jar "${SNPEFF_JAR}" ann \
        -noStats \
        -noLog \
        -config "${SNPEFF_CFG}" \
        "${SNPEFF_DB}" \
        "${VCF_FOR_SNPEFF}" \
        2> "${SNPEFF_LOG}" \
        | bcftools annotate --rename-chrs "${SNPEFF_RCRS}" \
        | bgzip -c > "${VCF_SNPEFF}"
    index_vcf "${VCF_SNPEFF}"

    rm -f "${VCF_FOR_SNPEFF}" "${VCF_FOR_SNPEFF}.tbi" \
          "${SNPEFF_CHRM}" "${SNPEFF_RCRS}"
    log "[${SAMPLE}] SnpEff completado"

    # -------------------------------------------------------------------------
    # 3.3. Anotacion ClinVar (CLNSIG + CLNDN)
    # -------------------------------------------------------------------------
    echo "Paso 3.3 - Anotación por ClinVar: [${SAMPLE}]"
    log "Paso 3.3 - Anotación por ClinVar: [${SAMPLE}]"
    VCF_CLINVAR="${SAMPLE_OUT}/${SAMPLE}_clinvar.vcf.gz"

    bcftools annotate \
        --annotations "${CLINVAR_DB}" \
        --columns "INFO/CLNSIG,INFO/CLNDN" \
        --output-type z \
        --output "${VCF_CLINVAR}" \
        "${VCF_SNPEFF}"
    index_vcf "${VCF_CLINVAR}"

    N_CLINVAR=$(bcftools view -H "${VCF_CLINVAR}" \
        | awk -F'\t' '$8 ~ /CLNSIG/' | wc -l)
    log "[${SAMPLE}] Variantes con anotacion ClinVar: ${N_CLINVAR}"

    # -------------------------------------------------------------------------
    # 3.4. Anotacion gnomAD mitocondrial (AF_hom, AF_het, AC_hom, AC_het)
    # -------------------------------------------------------------------------
    echo "Paso 3.4 - Anotación con gnomAD mitocondrial: [${SAMPLE}]"
    log "Paso 3.4 - Anotación con gnomAD mitocondrial: [${SAMPLE}]"

    VCF_GNOMAD="${SAMPLE_OUT}/${SAMPLE}_gnomad.vcf.gz"

    bcftools annotate \
        --annotations "${GNOMAD_DB}" \
        --columns "INFO/AF_hom,INFO/AF_het,INFO/AC_hom,INFO/AC_het" \
        --output-type z \
        --output "${VCF_GNOMAD}" \
        "${VCF_CLINVAR}"
    index_vcf "${VCF_GNOMAD}"

    N_GNOMAD=$(bcftools view -H "${VCF_GNOMAD}" \
        | awk -F'\t' '$8 ~ /AF_hom/' | wc -l)
    log "[${SAMPLE}] Variantes con frecuencia gnomAD: ${N_GNOMAD}"

    # -------------------------------------------------------------------------
    # 3.5. Clasificacion de haplogrupo con HaploGrep3
    # HaploGrep3 acepta VCF como input y devuelve clasificacion en TSV
    # -------------------------------------------------------------------------
    echo "Paso 3.5 - Clasificación de haplogrupos (HaploGrep3): [${SAMPLE}]"
    log "Paso 3.5 - Clasificación de haplogrupos (HaploGrep3): [${SAMPLE}]"

    HAPLO_OUT="${RESULTS}/Haplogroups/${SAMPLE}_haplogroup.txt"
    HAPLO_LOG="${LOG_DIR}/${SAMPLE}_haplogrep_${TIMESTAMP}.log"

    # Filtrar indels complejos antes de clasificar con Haplogrep3
    VCF_FOR_HAPLO="${SAMPLE_OUT}/${SAMPLE}_haplo_input.vcf.gz"
    bcftools view -v snps "${VCF_GNOMAD}" -Oz -o "${VCF_FOR_HAPLO}"
    index_vcf "${VCF_FOR_HAPLO}"

    java -jar "${HAPLOGREP_JAR}" classify \
        --in "${VCF_FOR_HAPLO}" \
        --out "${HAPLO_OUT}" \
        --tree phylotree-rcrs@17.2 \
        --extend-report \
        > "${HAPLO_LOG}" 2>&1

    rm -f "${VCF_FOR_HAPLO}" "${VCF_FOR_HAPLO}.tbi"

    HAPLOGROUP=$(awk -F'\t' 'NR==2 {gsub(/"/, "", $2); print $2}' "${HAPLO_OUT}" 2>/dev/null || echo "N/A")
    HAPLO_QUALITY=$(awk -F'\t' 'NR==2 {gsub(/"/, "", $4); print $4}' "${HAPLO_OUT}" 2>/dev/null || echo "N/A")
    log "[${SAMPLE}] Haplogrupo: ${HAPLOGROUP} (calidad: ${HAPLO_QUALITY})"

    # -------------------------------------------------------------------------
    # 3.6. Exportacion a TSV anotado
    # Extrae campos clave del VCF final anotado a formato tabular para el informe
    # Campos: CHROM, POS, REF, ALT, AF, DP, ANN (SnpEff), CLNSIG, CLNDN,
    #         AF_hom, AF_het, AC_hom, AC_het
    # -------------------------------------------------------------------------
    echo "Paso 3.6 - Exportar a TSV: [${SAMPLE}]"
    log "Paso 3.6 - Exportar a TSV: [${SAMPLE}]"
    TSV_OUT="${RESULTS}/TSV/${SAMPLE}_annotated.tsv"

    # Cabecera TSV
     {
        echo "# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"
        echo "# TSV ANOTADO - mtDNA ONT | Muestra: ${SAMPLE}"
        echo "# Generado: $(date)"
        echo "# Referencia: ${CONTIG} (rCRS)"
        echo "#"
        echo "# VERSIONES DE BASES DE DATOS:"
        printf "# %-25s %-15s %-55s %s\n" "Base_de_datos" "Version" "Fuente" "Fecha_descarga"
            # Solo el bloque mas reciente del log (evita duplicados de ejecuciones anteriores)
        tac "${DB_VERSION_LOG}" \
            | awk '/^$/{if(found) exit} /\t/{found=1; print}' \
            | tac \
            | awk -F'\t' 'NF>=4 {printf "# %-25s %-15s %-55s %s\n", $1, $2, $3, $4}'
        echo "#"
        echo "# LEYENDA DE COLUMNAS:"
        echo "# CHROM      : Cromosoma / contig de referencia (${CONTIG})"
        echo "# POS        : Posicion en la referencia (rCRS)"
        echo "# REF        : Alelo de referencia"
        echo "# ALT        : Alelo alternativo"
        echo "# AF_ALT     : Frecuencia del alelo alternativo (heteroplasma)= AD_ALT / DP_TOTAL"
        echo "# DP_TOTAL   : Cobertura total en la posicion"
        echo "# GENE       : Gen mitocondrial afectado (anotacion SnpEff)"
        echo "# EFFECT     : Efecto funcional predicho (SnpEff)"
        echo "# IMPACT     : Impacto funcional: HIGH / MODERATE / LOW / MODIFIER"
        echo "# HGVS_C     : Nomenclatura HGVS a nivel de ADN (c.)"
        echo "# HGVS_P     : Nomenclatura HGVS a nivel de proteina (p.); '.' si no aplica"
        echo "# CLNSIG     : Clasificacion clinica ClinVar (Pathogenic, Benign, VUS, etc.)"
        echo "# CLNDN      : Nombre de la enfermedad asociada en ClinVar"
        echo "# AF_hom     : Frecuencia de homoplasmia en gnomAD v3.1"
        echo "# AF_het     : Frecuencia de heteroplasmia en gnomAD v3.1"
        echo "# AC_hom     : Numero de individuos homoplasmaticos en gnomAD v3.1"
        echo "# AC_het     : Numero de individuos heteroplasmaticos en gnomAD v3.1"
        echo "# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"
        echo -e "CHROM\tPOS\tREF\tALT\tAF_ALT\tDP_TOTAL\tGENE\tEFFECT\tIMPACT\tHGVS_C\tHGVS_P\tCLNSIG\tCLNDN\tAF_hom\tAF_het\tAC_hom\tAC_het"
    } > "${TSV_OUT}"

    # Extraer campos con bcftools query
    # ANN de SnpEff, formato: Allele|Effect|Impact|GeneName|...|HGVS.c|HGVS.p|...
    # Se extraen campos con awk para las posiciones 4,2,3,10,11:
    bcftools query \
        -f '%CHROM\t%POS\t%REF\t%ALT\t[%AD]\t%INFO/ANN\t%INFO/CLNSIG\t%INFO/CLNDN\t%INFO/AF_hom\t%INFO/AF_het\t%INFO/AC_hom\t%INFO/AC_het\n' \
        "${VCF_GNOMAD}" \
    | awk 'BEGIN {OFS="\t"} {
        # Calcular AF_ALT y DP_TOTAL desde el campo AD (ref,alt)
        split($5, ad, ",")
        ad_ref = (ad[1] != "." && ad[1] != "") ? ad[1]+0 : 0
        ad_alt = (ad[2] != "." && ad[2] != "") ? ad[2]+0 : 0
        dp_total = ad_ref + ad_alt
        af_alt = (dp_total > 0) ? (ad_alt / dp_total) : "."
        if (af_alt != ".") af_alt = sprintf("%.4f", af_alt)
 
        # Parsear ANN de SnpEff (primer transcrito)
        n = split($6, transcripts, ",")
        split(transcripts[1], ann, "|")
        gene   = (ann[4] != "") ? ann[4] : "."
        effect = (ann[2] != "") ? ann[2] : "."
        impact = (ann[3] != "") ? ann[3] : "."
        hgvs_c = (ann[10] != "") ? ann[10] : "."
        hgvs_p = (ann[11] != "") ? ann[11] : "."
 
        print $1,$2,$3,$4,af_alt,dp_total,gene,effect,impact,hgvs_c,hgvs_p,$7,$8,$9,$10,$11,$12
    }' >> "${TSV_OUT}"

    N_TSV=$(tail -n +2 "${TSV_OUT}" | wc -l)
    log "[${SAMPLE}] TSV exportado: ${N_TSV} variantes → ${TSV_OUT}"

    # -------------------------------------------------------------------------
    # 3.7. Estadisticas de anotacion para el informe

    N_HIGH=$(awk -F'\t' 'NR>1 && $9=="HIGH"'    "${TSV_OUT}" | wc -l)
    N_MOD=$(awk  -F'\t' 'NR>1 && $9=="MODERATE"' "${TSV_OUT}" | wc -l)
    N_LOW=$(awk  -F'\t' 'NR>1 && $9=="LOW"'     "${TSV_OUT}" | wc -l)
    N_MOD_IMPACT=$(awk  -F'\t' 'NR>1 && $9=="MODIFIER"' "${TSV_OUT}" | wc -l)
    N_PATH=$(awk -F'\t' 'NR>1 && ($12 ~ /Pathogenic/ || $12 ~ /Likely_pathogenic/)' \
        "${TSV_OUT}" | wc -l)
    N_BENIGN=$(awk -F'\t' 'NR>1 && ($12 ~ /Benign/ || $12 ~ /Likely_benign/)' \
        "${TSV_OUT}" | wc -l)
    N_VUS=$(awk -F'\t' 'NR>1 && $12 ~ /Uncertain/' "${TSV_OUT}" | wc -l)
    N_COMMON=$(awk -F'\t' 'NR>1 && $14!="." && $14+0 >= 0.01' "${TSV_OUT}" | wc -l)

    {
        echo "-------------------------------------------------------------------------"
        echo "  MUESTRA: ${SAMPLE}"
        echo "-------------------------------------------------------------------------"
        echo ""
        printf "  %-35s %s\n" "Haplogrupo:"              "${HAPLOGROUP}"
        echo ""
        echo "  --- Clasificacion HaploGrep3 (matriz de polimorfismos) ---"
        # Parsear el archivo de haplogrupo y formatear como tabla
        awk -F'\t' 'NR==1 {
            # Imprimir cabecera de tabla
            printf "  %-15s %-8s %-8s %-12s %s\n", "SampleID", "Haplogr.", "Calidad", "Rango", "Polimorfismos encontrados"
            printf "  %-15s %-8s %-8s %-12s %s\n", "---------------", "--------", "--------", "------------", "------------------------"
        }
        NR==2 {
            # Limpiar comillas
            gsub(/"/, "")
            sampleid = $1
            haplogroup = $2
            quality   = $4
            range_val = $5
            found     = $7
            remaining = $8
            not_found = $6
 
            printf "  %-15s %-8s %-8s %-12s\n", sampleid, haplogroup, quality, range_val
            printf "\n"
            if (found != "") {
                printf "  Polimorfismos haplogrupo encontrados:\n"
                n = split(found, poly, " ")
                line = "    "
                for (i=1; i<=n; i++) {
                    line = line poly[i] "  "
                    if (length(line) > 75) { printf "%s\n", line; line = "    " }
                }
                if (line != "    ") printf "%s\n", line
            }
            if (not_found != "") {
                printf "\n  Polimorfismos NO encontrados (ausentes en muestra):\n"
                printf "    %s\n", not_found
            }
            if (remaining != "") {
                printf "\n  Polimorfismos privados / restantes:\n"
                n = split(remaining, rem, " ")
                line = "    "
                for (i=1; i<=n; i++) {
                    line = line rem[i] "  "
                    if (length(line) > 75) { printf "%s\n", line; line = "    " }
                }
                if (line != "    ") printf "%s\n", line
            }
        }' "${HAPLO_OUT}"
        echo ""
        printf "  %-35s %s\n" "Calidad clasificacion:"   "${HAPLO_QUALITY}"
        echo ""
        echo "  --- Variantes anotadas ---"
        printf "  %-35s %6s\n" "Total variantes (AF >= ${AF_MIN}):" "${N_TSV}"
        printf "  %-35s %6s\n" "Con anotacion ClinVar:"             "${N_CLINVAR}"
        printf "  %-35s %6s\n" "Con frecuencia gnomAD:"             "${N_GNOMAD}"
        echo ""
        echo "  --- Impacto funcional (SnpEff) ---"
        printf "  %-35s %6s\n" "HIGH:"     "${N_HIGH}"
        printf "  %-35s %6s\n" "MODERATE:" "${N_MOD}"
        printf "  %-35s %6s\n" "LOW:"      "${N_LOW}"
        printf "  %-35s %6s\n" "MODIFIER:" "${N_MOD_IMPACT}"
        echo ""
        echo "  --- Clasificacion ClinVar ---"
        printf "  %-35s %6s\n" "Patogenicas / Prob. patogenicas:" "${N_PATH}"
        printf "  %-35s %6s\n" "Benignas / Prob. benignas:"       "${N_BENIGN}"
        printf "  %-35s %6s\n" "VUS (significado incierto):"      "${N_VUS}"
        echo ""
        echo "  --- Frecuencia poblacional (gnomAD) ---"
        printf "  %-35s %6s\n" "Variantes comunes (AF_hom >= 1%):" "${N_COMMON}"
        echo ""
        echo "  --- Archivos generados ---"
        printf "  %-35s %s\n" "VCF anotado final:"   "${VCF_GNOMAD}"
        printf "  %-35s %s\n" "TSV clinico:"          "${TSV_OUT}"
        printf "  %-35s %s\n" "Haplogrupo:"           "${HAPLO_OUT}"
        echo ""
    } >> "${REPORT_ANNOT}"

    log "[${SAMPLE}] Anotacion completada"
done

# --------------------------------------------------------------------------------------------
# 4. CIERRE DEL INFORME
# --------------------------------------------------------------------------------------------

{
    echo "-------------------------------------------------------------------------"
    echo "  FIN DEL INFORME DE ANOTACION"
    echo "  Completado: $(date)"
    echo "  Informes HTML: ${OUTPUT}/html/"
    echo "-------------------------------------------------------------------------"
} >> "${REPORT_ANNOT}"

echo ""
log "Pipeline de anotacion completado: $(date)"
log "Informe resumen: ${REPORT_ANNOT}"
log "Informes HTML:   ${OUTPUT}/html/"

# --------------------------------------------------------------------------------------------