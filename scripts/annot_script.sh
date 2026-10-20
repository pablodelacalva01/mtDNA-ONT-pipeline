#!/usr/bin/env bash
# --------------------------------------------------------------------------------------------
# Autor:   Pablo de la Calva Castineira
# Fecha:   06/05/2026
# Titulo:  "Analisis bioinformatico y su implementacion en la practica clinica
#           diaria de datos de NGS de lecturas largas mediante tecnologia Nanopore"
# Version: 1.3
# Descripcion: Pipeline de anotacion de variantes mitocondriales de alta confianza (ONT)
#              Parte 4: Anotacion de variantes mtDNA
#
# Input:   VCF union de las 4 combinaciones por muestra (MC, MM, NC, NM)
#          Nota: Generado en Parte 3: comparison/<muestra>/<muestra>_all_callers_merged.vcf.gz
#
# Pipeline:
#   1. Identificacion de muestras (VCF union + BAMs de minimap2 y ngmlr)
#   2. Soporte por combinacion (N_SUPPORT, CALLED_BY) y AF/DP por alineador desde los BAM
#   3. Anotacion funcional con SnpEff (gen, region, impacto, HGVS)
#   4. Anotacion con ClinVar (CLNSIG, CLNDN)
#   5. Anotacion con gnomAD mt (AF_hom, AF_het, AC_hom, AC_het)
#   6. Clasificacion de haplogrupo con HaploGrep3
#   7. Exportacion a TSV anotado por muestra
#   8. Informe de texto resumen por muestra
#
# Cambios v1.3 (respecto a v1.2):
#   - Input: union de las 4 combinaciones (sin eleccion de par por Jaccard, sin filtro de AF)
#   - Nuevas etiquetas INFO / columnas TSV: soporte por combinacion y AF/DP por alineador
#   - HaploGrep3 recibe un VCF de una sola muestra (SNVs con soporte >= HAPLO_MIN_SUPPORT)
#
# Entorno conda: ont_annot
# Prerequisitos: parte4_setup.sh ejecutado (HaploGrep3 + DBs instalados)
# --------------------------------------------------------------------------------------------

set -euo pipefail

# --------------------------------------------------------------------------------------------
# 0. CONFIGURACION
# Variables de entorno
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
VARIANTS_DIR="${PROJECT_DIR}/Data/Processed/04.Variants"
OUTPUT="${PROJECT_DIR}/Data/Processed/05.Annotation"
RESULTS="${PROJECT_DIR}/Results"
DB_DIR="${PROJECT_DIR}/Data/Annotation_DBs"
LOG_DIR="${PROJECT_DIR}/Data/Processed/logs"
DB_VERSION_LOG="${DB_DIR}/db_versions.log"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
THREADS=4

# Referencia
CONTIG="NC_012920.1"

# Bases de datos
CLINVAR_DB="${DB_DIR}/clinvar_mt_rCRS.vcf.gz"
GNOMAD_DB="${DB_DIR}/gnomad_mt_rCRS.vcf.gz"

# HaploGrep3
HAPLOGREP_JAR="${PROJECT_DIR}/Tools/haplogrep3/haplogrep3.jar"

# SnpEff
SNPEFF_DB="GRCh38.99"                                                        
SNPEFF_DIR="${CONDA_PREFIX}/share/snpeff-5.2-3"                              
SNPEFF_CFG="${SNPEFF_DIR}/snpEff.config"
SNPEFF_JAR="${SNPEFF_DIR}/snpEff.jar"
SNPEFF_DATA_DIR="${SNPEFF_DIR}/data"    

# Alineamientos (BAM post-MarkDuplicates) y referencia: AF/DP por alineador
ALIGN_DIR="${PROJECT_DIR}/Data/Processed/03.Alignment/markdup"
REF_FASTA="${PROJECT_DIR}/Data/Reference_Genome/mtDNA.fasta"

# Filtros del pileup (igualar a los usados en Mutserve2: MAPQ / BASEQ)
PILEUP_MAPQ="20"
PILEUP_BASEQ="12"

# Soporte minimo (n.o de combinaciones) para que un SNV entre en la clasificacion de haplogrupo.
# Solo afecta a HaploGrep3; la anotacion y el TSV incluyen TODAS las variantes de la union.
HAPLO_MIN_SUPPORT="2"

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

# Pileup por posicion (bcftools mpileup) sobre un BAM. Salida TSV: POS REF ALT_list DP AD_list
#   -d / -L altos: el limite por defecto (250) submuestrearia la cobertura del mtDNA
#   -B: sin BAQ | -I: sin llamada de indels | -q / -Q: filtros de MAPQ y BASEQ
run_pileup() {
    local BAM="$1" REGIONS="$2" OUT="$3" LOG_FILE="$4"
    if [[ ! -s "${REGIONS}" ]]; then
        : > "${OUT}"
        return 0
    fi
    bcftools mpileup \
        -f "${REF_FASTA}" -R "${REGIONS}" \
        -q "${PILEUP_MAPQ}" -Q "${PILEUP_BASEQ}" \
        -B -I -d 1000000 -L 1000000 \
        -a FORMAT/AD,FORMAT/DP \
        -Ou "${BAM}" 2>> "${LOG_FILE}" \
        | bcftools query -f '%POS\t%REF\t%ALT\t[%DP]\t[%AD]\n' - > "${OUT}" \
        || err "bcftools mpileup fallo sobre ${BAM} (ver ${LOG_FILE})"
}

# Cabecera INFO de las etiquetas que se anaden al VCF union
write_tag_header() {
cat > "$1" <<'EOF'
##INFO=<ID=N_SUPPORT,Number=1,Type=Integer,Description="Numero de combinaciones alineador+caller (MC,MM,NC,NM) que detectan la variante">
##INFO=<ID=CALLED_BY,Number=1,Type=String,Description="Combinaciones que detectan la variante, separadas por +">
##INFO=<ID=CROSS_CALLER,Number=1,Type=String,Description="SI si la variante es detectada por Clair3 y por Mutserve2">
##INFO=<ID=VAR_TYPE,Number=1,Type=String,Description="SNV, INDEL u OTRO">
##INFO=<ID=AF_MINIMAP2,Number=1,Type=Float,Description="Fraccion alelica del ALT en el BAM de minimap2">
##INFO=<ID=DP_MINIMAP2,Number=1,Type=Integer,Description="Profundidad en el BAM de minimap2">
##INFO=<ID=REFC_MINIMAP2,Number=1,Type=Integer,Description="Lecturas REF en el BAM de minimap2">
##INFO=<ID=ALTC_MINIMAP2,Number=1,Type=Integer,Description="Lecturas ALT en el BAM de minimap2">
##INFO=<ID=AF_NGMLR,Number=1,Type=Float,Description="Fraccion alelica del ALT en el BAM de ngmlr">
##INFO=<ID=DP_NGMLR,Number=1,Type=Integer,Description="Profundidad en el BAM de ngmlr">
##INFO=<ID=REFC_NGMLR,Number=1,Type=Integer,Description="Lecturas REF en el BAM de ngmlr">
##INFO=<ID=ALTC_NGMLR,Number=1,Type=Integer,Description="Lecturas ALT en el BAM de ngmlr">
##INFO=<ID=AF_SOURCE,Number=1,Type=String,Description="Origen de AF/DP: pileup (SNVs), Clair3 (indels) o . si no disponible">
EOF
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
require tabix
require awk

# Verificar bases de datos
for DB_FILE in "${CLINVAR_DB}" "${GNOMAD_DB}"; do
    [[ -f "${DB_FILE}" ]] \
        || err "Base de datos no encontrada: ${DB_FILE}. Ejecuta primero annot_setup.sh"
    [[ -f "${DB_FILE}.tbi" ]] \
        || err "Indice .tbi no encontrado para: ${DB_FILE}. Re-ejecuta annot_setup.sh"
done

# Verificar HaploGrep3
[[ -f "${HAPLOGREP_JAR}" ]] \
    || err "HaploGrep3 no encontrado en ${HAPLOGREP_JAR}. Ejecuta primero annot_setup.sh"

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
for CMP_SAMPLE_DIR in "${VARIANTS_DIR}/comparison/"*/; do
    SAMPLE=$(basename "${CMP_SAMPLE_DIR}")
    MERGED_VCF="${CMP_SAMPLE_DIR%/}/${SAMPLE}_all_callers_merged.vcf.gz"
    MISSING=""
    [[ -f "${MERGED_VCF}" ]] || MISSING+=" VCF_union"
    for ALN in minimap2 ngmlr; do
        BAM="${ALIGN_DIR}/${ALN}/${SAMPLE}/${SAMPLE}_${ALN}_markdup.bam"
        [[ -f "${BAM}" ]] || MISSING+=" BAM_${ALN}"
        [[ -f "${BAM}.bai" || -f "${BAM%.bam}.bai" || -f "${BAM}.csi" ]] || MISSING+=" indice_${ALN}"
    done
    if [[ -z "${MISSING}" ]]; then
        SAMPLES+=("${SAMPLE}")
    else
        log "AVISO: ${SAMPLE} omitida, falta:${MISSING}"
    fi
done

if [[ ${#SAMPLES[@]} -eq 0 ]]; then
    err "No se encontraron muestras con VCF union (all_callers_merged) y BAMs en ${VARIANTS_DIR}/comparison/"
fi

mkdir -p "${LOG_DIR}"

log "Muestras detectadas: ${#SAMPLES[@]}"

# Crear directorios de salida
mkdir -p "${RESULTS}/Variants_Annotated" "${RESULTS}/Haplogroups" "${RESULTS}/TSV" "${RESULTS}/AnnotationReports"

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
    echo "Input:          union de las 4 combinaciones (MC, MM, NC, NM), sin filtro de AF"
    echo "AF / DP:        pileup por alineador (MAPQ >= ${PILEUP_MAPQ}, BASEQ >= ${PILEUP_BASEQ})"
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

    MERGED_VCF="${VARIANTS_DIR}/comparison/${SAMPLE}/${SAMPLE}_all_callers_merged.vcf.gz"

    # >>> BLOQUE_SOPORTE_AF_INICIO
    # -------------------------------------------------------------------------
    # 3.1. VCF union de las 4 combinaciones (MC, MM, NC, NM) + etiquetas INFO:
    #   - N_SUPPORT / CALLED_BY / CROSS_CALLER: desde FORMAT/GT del merge
    #   - AF/DP/AD por ALINEADOR (minimap2, ngmlr):
    #       SNVs   -> pileup del BAM con los mismos filtros para todas las variantes
    #       indels -> FORMAT/AD de Clair3 (MC para minimap2, NC para ngmlr)
    # NO se filtra ninguna variante: se anota la union completa.
    # -------------------------------------------------------------------------
    echo "Paso 3.1 - Union de las 4 combinaciones, soporte y AF por alineador: [${SAMPLE}]"
    log "Paso 3.1 - Union de las 4 combinaciones, soporte y AF por alineador: [${SAMPLE}]"

    VCF_UNION="${SAMPLE_OUT}/${SAMPLE}_union_tagged.vcf.gz"
    TAGS_TSV="${SAMPLE_OUT}/${SAMPLE}_support_tags.tsv"
    TAGS_HDR="${SAMPLE_OUT}/${SAMPLE}_tags_header.tmp"
    SITES_TMP="${SAMPLE_OUT}/${SAMPLE}_sites.tmp"
    REGIONS_TMP="${SAMPLE_OUT}/${SAMPLE}_regions.tmp"
    PILEUP_MM2="${SAMPLE_OUT}/${SAMPLE}_pileup_minimap2.tmp"
    PILEUP_NGMLR="${SAMPLE_OUT}/${SAMPLE}_pileup_ngmlr.tmp"
    PILEUP_LOG="${LOG_DIR}/${SAMPLE}_mpileup_${TIMESTAMP}.log"
    : > "${PILEUP_LOG}"

    # Las columnas del merge se usan por posicion: deben ser MC MM NC NM
    MERGED_SAMPLES=$(bcftools query -l "${MERGED_VCF}" | paste -sd' ' -)
    [[ "${MERGED_SAMPLES}" == "MC MM NC NM" ]] \
        || err "[${SAMPLE}] Orden de muestras inesperado en ${MERGED_VCF}: '${MERGED_SAMPLES}' (esperado: MC MM NC NM)"

    # Sitios: CHROM POS REF ALT + [GT DP AD] por combinacion (12 columnas a partir de la 5)
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT\t%DP\t%AD]\n' "${MERGED_VCF}" > "${SITES_TMP}"

    # Posiciones de los SNVs (para el pileup)
    awk -F'\t' '$3 ~ /^[ACGTacgt]$/ && $4 ~ /^[ACGTacgt]$/ {print $1 "\t" $2}' "${SITES_TMP}" \
        | sort -u -k2,2n > "${REGIONS_TMP}"

    run_pileup "${ALIGN_DIR}/minimap2/${SAMPLE}/${SAMPLE}_minimap2_markdup.bam" "${REGIONS_TMP}" "${PILEUP_MM2}" "${PILEUP_LOG}"
    run_pileup "${ALIGN_DIR}/ngmlr/${SAMPLE}/${SAMPLE}_ngmlr_markdup.bam"       "${REGIONS_TMP}" "${PILEUP_NGMLR}" "${PILEUP_LOG}"

    # Tabla de etiquetas: CHROM POS REF ALT N_SUPPORT CALLED_BY CROSS_CALLER VAR_TYPE
    #                     AF DP REFC ALTC (minimap2) | AF DP REFC ALTC (ngmlr) | AF_SOURCE
    awk -v f_mm2="${PILEUP_MM2}" -v f_ng="${PILEUP_NGMLR}" '
        BEGIN {
            FS = OFS = "\t"
            lab[1] = "MC"; lab[2] = "MM"; lab[3] = "NC"; lab[4] = "NM"
        }
        # Pileups: POS REF ALT_list DP AD_list  (M = minimap2, N = ngmlr)
        FILENAME == f_mm2 || FILENAME == f_ng {
            st = (FILENAME == f_mm2) ? "M" : "N"
            pos = $1
            has[st, pos] = 1
            dp[st, pos]  = $4
            nA = split($3, alts, ",")
            nD = split($5, ads, ",")
            rc[st, pos]  = ads[1]
            cnt[st, pos, $2] = ads[1]
            for (k = 1; k <= nA; k++)
                if (alts[k] !~ /^</) cnt[st, pos, alts[k]] = ads[k + 1]
            next
        }
        function pile(st, pos, alt) {
            if (!((st, pos) in has)) { DP[st] = "."; RC[st] = "."; AC[st] = "."; AF[st] = "."; return }
            DP[st] = dp[st, pos]
            RC[st] = rc[st, pos]
            AC[st] = ((st, pos, alt) in cnt) ? cnt[st, pos, alt] : 0
            AF[st] = (DP[st] + 0 > 0) ? sprintf("%.4f", AC[st] / DP[st]) : "."
        }
        function clair(st, i,    ad, x, d) {
            DP[st] = "."; RC[st] = "."; AC[st] = "."; AF[st] = "."
            ad = $(7 + 3 * (i - 1)); d = $(6 + 3 * (i - 1))
            if (ad == "." || ad == "" || d == "." || d == "" || d + 0 <= 0) return
            split(ad, x, ",")
            if (x[2] == "") return
            DP[st] = d; RC[st] = x[1]; AC[st] = x[2]; AF[st] = sprintf("%.4f", x[2] / d)
        }
        # Sitios del merge
        {
            chrom = $1; pos = $2; r = $3; a = $4
            n = 0; cb = ""; cl = 0; mu = 0
            for (i = 1; i <= 4; i++) {
                if ($(5 + 3 * (i - 1)) ~ /1/) {
                    n++
                    cb = cb (cb == "" ? "" : "+") lab[i]
                    if (i == 1 || i == 3) cl = 1; else mu = 1
                }
            }
            if (cb == "") cb = "."
            cross = (cl && mu) ? "SI" : "NO"

            if (index(a, ",") > 0 || a ~ /^</)          type = "OTRO"
            else if (length(r) == 1 && length(a) == 1)  type = "SNV"
            else                                        type = "INDEL"

            if (type == "SNV") {
                pile("M", pos, a); pile("N", pos, a); src = "pileup"
            } else if (type == "INDEL") {
                clair("M", 1); clair("N", 3); src = "Clair3"
            } else {
                DP["M"] = RC["M"] = AC["M"] = AF["M"] = "."
                DP["N"] = RC["N"] = AC["N"] = AF["N"] = "."
                src = "."
            }
            if (AF["M"] == "." && AF["N"] == ".") src = "."

            print chrom, pos, r, a, n, cb, cross, type, \
                  AF["M"], DP["M"], RC["M"], AC["M"], \
                  AF["N"], DP["N"], RC["N"], AC["N"], src
        }' "${PILEUP_MM2}" "${PILEUP_NGMLR}" "${SITES_TMP}" > "${TAGS_TSV}"

    [[ -s "${TAGS_TSV}" ]] || err "[${SAMPLE}] El merge ${MERGED_VCF} no contiene variantes."

    bgzip -c "${TAGS_TSV}" > "${TAGS_TSV}.gz"
    tabix -f -s1 -b2 -e2 "${TAGS_TSV}.gz"
    write_tag_header "${TAGS_HDR}"

    bcftools annotate \
        --annotations "${TAGS_TSV}.gz" \
        --header-lines "${TAGS_HDR}" \
        --columns "CHROM,POS,REF,ALT,INFO/N_SUPPORT,INFO/CALLED_BY,INFO/CROSS_CALLER,INFO/VAR_TYPE,INFO/AF_MINIMAP2,INFO/DP_MINIMAP2,INFO/REFC_MINIMAP2,INFO/ALTC_MINIMAP2,INFO/AF_NGMLR,INFO/DP_NGMLR,INFO/REFC_NGMLR,INFO/ALTC_NGMLR,INFO/AF_SOURCE" \
        --output-type z \
        --output "${VCF_UNION}" \
        "${MERGED_VCF}"
    index_vcf "${VCF_UNION}"

    N_MERGED=$(count_variants "${MERGED_VCF}")
    N_UNION=$(count_variants "${VCF_UNION}")
    [[ "${N_UNION}" -eq "${N_MERGED}" ]] \
        || err "[${SAMPLE}] La union etiquetada (${N_UNION}) no coincide con el merge (${N_MERGED})."
    log "[${SAMPLE}] Variantes en la union (sin filtrar): ${N_UNION}"

    rm -f "${SITES_TMP}" "${REGIONS_TMP}" "${PILEUP_MM2}" "${PILEUP_NGMLR}" "${TAGS_HDR}"
    # <<< BLOQUE_SOPORTE_AF_FIN

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
        "${VCF_UNION}" -Oz -o "${VCF_FOR_SNPEFF}"
    index_vcf "${VCF_FOR_SNPEFF}"

    # Anotar con SnpEff y restaurar NC_012920.1
    java -Xmx2g -jar "${SNPEFF_JAR}" ann \
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
    # 3.3. Anotacion ClinVar (ALLELEID + CLNSIG + CLNDN)
    # -------------------------------------------------------------------------
    echo "Paso 3.3 - Anotación por ClinVar: [${SAMPLE}]"
    log "Paso 3.3 - Anotación por ClinVar: [${SAMPLE}]"
    VCF_CLINVAR="${SAMPLE_OUT}/${SAMPLE}_clinvar.vcf.gz"

    bcftools annotate \
        --annotations "${CLINVAR_DB}" \
        --columns "INFO/ALLELEID,INFO/CLNSIG,INFO/CLNDN" \
        --threads "${THREADS}" \
        --output-type z \
        --output "${VCF_CLINVAR}" \
        "${VCF_SNPEFF}"
    index_vcf "${VCF_CLINVAR}"

    N_CLINVAR=$(bcftools query -f '%INFO/ALLELEID\n' "${VCF_CLINVAR}" \
        | awk '$1 != "." && $1 != "" {count++} END {print count+0}')
    log "[${SAMPLE}] Variantes con anotacion ClinVar: ${N_CLINVAR}"
    if [[ "${N_CLINVAR}" -eq 0 ]]; then
        log "AVISO: [${SAMPLE}] ninguna variante casa con ClinVar. Es esperable si solo hay polimorfismos de haplogrupo; si no, revisa la BD."
    fi

    # -------------------------------------------------------------------------
    # 3.4. Anotacion gnomAD mitocondrial (AF_hom, AF_het, AC_hom, AC_het)
    # -------------------------------------------------------------------------
    echo "Paso 3.4 - Anotación con gnomAD mitocondrial: [${SAMPLE}]"
    log "Paso 3.4 - Anotación con gnomAD mitocondrial: [${SAMPLE}]"

    VCF_GNOMAD="${SAMPLE_OUT}/${SAMPLE}_gnomad.vcf.gz"

    bcftools annotate \
        --annotations "${GNOMAD_DB}" \
        --columns "INFO/AF_hom,INFO/AF_het,INFO/AC_hom,INFO/AC_het" \
        --threads "${THREADS}" \
        --output-type z \
        --output "${VCF_GNOMAD}" \
        "${VCF_CLINVAR}"
    index_vcf "${VCF_GNOMAD}"

    N_GNOMAD=$(bcftools query -f '%INFO/AF_hom\t%INFO/AF_het\t%INFO/AC_hom\t%INFO/AC_het\n' "${VCF_GNOMAD}" \
        | awk '{for (i=1; i<=4; i++) if ($i != "." && $i != "") {count++; break}} END {print count+0}')
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
    # La union tiene 4 columnas de muestra (MC MM NC NM): HaploGrep3 las trataria como 4 muestras.
    # Se genera un VCF de UNA muestra (GT=1) con los SNVs de la union con soporte >= HAPLO_MIN_SUPPORT.
    bcftools view -G -v snps -i "INFO/N_SUPPORT>=${HAPLO_MIN_SUPPORT}" "${VCF_GNOMAD}" \
        | awk -v s="${SAMPLE}" 'BEGIN { OFS = "\t" }
            /^##/ { print; next }
            /^#CHROM/ { print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"; print $0, "FORMAT", s; next }
            { print $0, "GT", "1" }' \
        | bgzip -c > "${VCF_FOR_HAPLO}"
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
    # Campos: CHROM, POS, REF, ALT, DP, AD, AF, ANN, ClinVar, gnomAD.
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
        echo "# CHROM        : Contig de referencia (${CONTIG})"
        echo "# POS          : Posicion en la referencia (rCRS)"
        echo "# REF / ALT    : Alelo de referencia / alternativo (una fila por alelo ALT)"
        echo "# AF_ALT       : Fraccion alelica del ALT (heteroplasmia): media de AF_MINIMAP2 y AF_NGMLR"
        echo "# COVERAGE     : Profundidad total en la posicion: media de DP_MINIMAP2 y DP_NGMLR"
        echo "# AD           : Lecturas ref,alt: media de ambos alineadores, redondeada"
        echo "# GENE         : Gen del efecto SnpEff seleccionado (mayor impacto; upstream/downstream penalizados)"
        echo "# EFFECT       : Efecto funcional predicho (SnpEff)"
        echo "# IMPACT       : HIGH / MODERATE / LOW / MODIFIER"
        echo "# HGVS_C       : HGVS a nivel de ADN codificante (c./n.) segun SnpEff"
        echo "# HGVS_P       : HGVS a nivel de proteina (p.); '.' si no aplica"
        echo "# CLINVAR_ID   : ClinVar Allele ID (INFO/ALLELEID)"
        echo "# CLNSIG       : Clasificacion clinica ClinVar"
        echo "# CLNDN        : Enfermedad asociada en ClinVar"
        echo "# AF_hom       : Frecuencia de homoplasmia en gnomAD v3.1"
        echo "# AF_het       : Frecuencia de heteroplasmia en gnomAD v3.1"
        echo "# AC_hom       : Individuos homoplasmaticos en gnomAD v3.1"
        echo "# AC_het       : Individuos heteroplasmaticos en gnomAD v3.1"
        echo "# VAR_TYPE     : SNV / INDEL / OTRO (multialelico o simbolico)"
        echo "# N_SUPPORT    : Numero de combinaciones (0-4) que detectan la variante"
        echo "# CALLED_BY    : Combinaciones que la detectan: MC=minimap2+Clair3, MM=minimap2+Mutserve2, NC=ngmlr+Clair3, NM=ngmlr+Mutserve2"
        echo "# CROSS_CALLER : SI si la detectan Clair3 y Mutserve2 (con cualquier alineador), NO en caso contrario"
        echo "# AF_MINIMAP2 / DP_MINIMAP2 / AD_MINIMAP2 : fraccion alelica, profundidad y lecturas ref,alt en el BAM de minimap2"
        echo "# AF_NGMLR    / DP_NGMLR    / AD_NGMLR    : idem en el BAM de ngmlr"
        echo "# AF_SOURCE    : pileup = bcftools mpileup (MAPQ>=${PILEUP_MAPQ}, BASEQ>=${PILEUP_BASEQ}, SNVs) | Clair3 = FORMAT/AD de MC y NC (indels) | . = no disponible"
        printf "CHROM\tPOS\tREF\tALT\tAF_ALT\tCOVERAGE\tAD\tGENE\tEFFECT\tIMPACT\tHGVS_C\tHGVS_P\tCLINVAR_ID\tCLNSIG\tCLNDN\tAF_hom\tAF_het\tAC_hom\tAC_het\tVAR_TYPE\tN_SUPPORT\tCALLED_BY\tCROSS_CALLER\tAF_MINIMAP2\tDP_MINIMAP2\tAD_MINIMAP2\tAF_NGMLR\tDP_NGMLR\tAD_NGMLR\tAF_SOURCE\n"
    } > "${TSV_OUT}"

    # >>> BLOQUE_TSV_INICIO
    # Extraer campos con bcftools query
    # ANN de SnpEff, formato: Allele|Effect|Impact|GeneName|...|HGVS.c|HGVS.p|...
    # Las columnas 1-19 mantienen el orden y significado de v1.2; las nuevas van al final.
    # AF_ALT / COVERAGE / AD = media de los dos alineadores (minimap2, ngmlr).
    bcftools query \
        -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/ANN\t%INFO/ALLELEID\t%INFO/CLNSIG\t%INFO/CLNDN\t%INFO/AF_hom\t%INFO/AF_het\t%INFO/AC_hom\t%INFO/AC_het\t%INFO/VAR_TYPE\t%INFO/N_SUPPORT\t%INFO/CALLED_BY\t%INFO/CROSS_CALLER\t%INFO/AF_MINIMAP2\t%INFO/DP_MINIMAP2\t%INFO/REFC_MINIMAP2\t%INFO/ALTC_MINIMAP2\t%INFO/AF_NGMLR\t%INFO/DP_NGMLR\t%INFO/REFC_NGMLR\t%INFO/ALTC_NGMLR\t%INFO/AF_SOURCE\n' \
        "${VCF_GNOMAD}" \
    | awk '
        function mean(a, b, fmt,    s, k) {
            s = 0; k = 0
            if (a != "" && a != ".") { s += a; k++ }
            if (b != "" && b != ".") { s += b; k++ }
            return (k > 0) ? sprintf(fmt, s / k) : "."
        }
        function pair(r, a) {
            return (r != "" && r != "." && a != "" && a != ".") ? r "," a : "."
        }
        BEGIN { OFS = "\t" }
        {
            # Parsear ANN de SnpEff (primer transcrito)
            n = split($5, transcripts, ",")
            split(transcripts[1], ann, "|")
            gene   = (ann[4] != "") ? ann[4] : "."
            effect = (ann[2] != "") ? ann[2] : "."
            impact = (ann[3] != "") ? ann[3] : "."
            hgvs_c = (ann[10] != "") ? ann[10] : "."
            hgvs_p = (ann[11] != "") ? ann[11] : "."

            clinvar_id = ($6 != "" && $6 != ".") ? $6 : "."
            clnsig = ($7 != "" && $7 != ".") ? $7 : "."
            clndn = ($8 != "" && $8 != ".") ? $8 : "."
            af_hom = ($9 != "" && $9 != ".") ? $9 : "."
            af_het = ($10 != "" && $10 != ".") ? $10 : "."
            ac_hom = ($11 != "" && $11 != ".") ? $11 : "."
            ac_het = ($12 != "" && $12 != ".") ? $12 : "."

            # Etiquetas de soporte y AF por alineador (3.1)
            vtype = $13; nsup = $14; cby = $15; cross = $16
            af_m = $17; dp_m = $18; r_m = $19; a_m = $20
            af_n = $21; dp_n = $22; r_n = $23; a_n = $24
            src = $25

            af = mean(af_m, af_n, "%.4f")
            dp = mean(dp_m, dp_n, "%.0f")
            ad = pair(mean(r_m, r_n, "%.0f"), mean(a_m, a_n, "%.0f"))

            print $1,$2,$3,$4,af,dp,ad,gene,effect,impact,hgvs_c,hgvs_p,clinvar_id,clnsig,clndn,af_hom,af_het,ac_hom,ac_het,vtype,nsup,cby,cross,af_m,dp_m,pair(r_m,a_m),af_n,dp_n,pair(r_n,a_n),src
        }' >> "${TSV_OUT}"
    # <<< BLOQUE_TSV_FIN

    N_TSV=$(awk -F'\t' '$1 != "CHROM" && $1 !~ /^#/ {count++} END {print count+0}' "${TSV_OUT}")
    log "[${SAMPLE}] TSV exportado: ${N_TSV} variantes → ${TSV_OUT}"

    # -------------------------------------------------------------------------
    # 3.7. Estadisticas de anotacion para el informe

    N_HIGH=$(awk -F'\t' '$1 != "CHROM" && $1 !~ /^#/ && $10=="HIGH"' "${TSV_OUT}" | wc -l)
    N_MOD=$(awk -F'\t' '$1 != "CHROM" && $1 !~ /^#/ && $10=="MODERATE"' "${TSV_OUT}" | wc -l)
    N_LOW=$(awk -F'\t' '$1 != "CHROM" && $1 !~ /^#/ && $10=="LOW"' "${TSV_OUT}" | wc -l)
    N_MOD_IMPACT=$(awk -F'\t' '$1 != "CHROM" && $1 !~ /^#/ && $10=="MODIFIER"' "${TSV_OUT}" | wc -l)
    N_PATH=$(awk -F'\t' '$1 != "CHROM" && $1 !~ /^#/ && ($14 ~ /Pathogenic/ || $14 ~ /Likely_pathogenic/)' \
        "${TSV_OUT}" | wc -l)
    N_BENIGN=$(awk -F'\t' '$1 != "CHROM" && $1 !~ /^#/ && ($14 ~ /Benign/ || $14 ~ /Likely_benign/)' \
        "${TSV_OUT}" | wc -l)
    N_VUS=$(awk -F'\t' '$1 != "CHROM" && $1 !~ /^#/ && $14 ~ /Uncertain/' "${TSV_OUT}" | wc -l)
    N_COMMON=$(awk -F'\t' '$1 != "CHROM" && $1 !~ /^#/ && $16!="." && $16+0 >= 0.01' "${TSV_OUT}" | wc -l)
    N_SUP4=$(awk -F'\t' '$1 != "CHROM" && $1 !~ /^#/ && $21==4' "${TSV_OUT}" | wc -l)
    N_SUP3=$(awk -F'\t' '$1 != "CHROM" && $1 !~ /^#/ && $21==3' "${TSV_OUT}" | wc -l)
    N_SUP2=$(awk -F'\t' '$1 != "CHROM" && $1 !~ /^#/ && $21==2' "${TSV_OUT}" | wc -l)
    N_SUP1=$(awk -F'\t' '$1 != "CHROM" && $1 !~ /^#/ && $21==1' "${TSV_OUT}" | wc -l)
    N_CROSS=$(awk -F'\t' '$1 != "CHROM" && $1 !~ /^#/ && $23=="SI"' "${TSV_OUT}" | wc -l)
    N_INDEL_T=$(awk -F'\t' '$1 != "CHROM" && $1 !~ /^#/ && $20=="INDEL"' "${TSV_OUT}" | wc -l)

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
        printf "  %-35s %6s\n" "Total variantes anotadas (union 4):" "${N_TSV}"
        printf "  %-35s %6s\n" "Con anotacion ClinVar:"             "${N_CLINVAR}"
        printf "  %-35s %6s\n" "Con frecuencia gnomAD:"             "${N_GNOMAD}"
        echo ""
        echo "  --- Soporte entre combinaciones (MC, MM, NC, NM) ---"
        printf "  %-35s %6s\n" "Detectadas por las 4:"              "${N_SUP4}"
        printf "  %-35s %6s\n" "Detectadas por 3:"                  "${N_SUP3}"
        printf "  %-35s %6s\n" "Detectadas por 2:"                  "${N_SUP2}"
        printf "  %-35s %6s\n" "Detectadas por 1 (revisar en IGV):" "${N_SUP1}"
        printf "  %-35s %6s\n" "Clair3 + Mutserve2 (cross-caller):" "${N_CROSS}"
        printf "  %-35s %6s\n" "Indels:"                            "${N_INDEL_T}"
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
    echo "-------------------------------------------------------------------------"
} >> "${REPORT_ANNOT}"

echo ""
log "Pipeline de anotacion completado: $(date)"
log "Informe resumen: ${REPORT_ANNOT}"

# --------------------------------------------------------------------------------------------