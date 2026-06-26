#!/usr/bin/env bash
# --------------------------------------------------------------------------------------------
# Autor:   Pablo de la Calva Castineira
# Fecha:   05/05/2026
# Titulo:  "Analisis bioinformatico y su implementacion en la practica clinica
#           diaria de datos de NGS de lecturas largas mediante tecnologia Nanopore"
# Version: 1.0
# Descripcion: Script de instalacion y actualizacion del entorno de anotacion (ont_annot)
#              Parte 4: Setup de HaploGrep3 y bases de datos de anotacion mtDNA
#
# Metodo de uso:
#   Instalacion completa:
#     bash annot_setup.sh
#   Actualizacion de bases de datos:
#     bash annot_setup.sh --update
#       Nota: Descarga las versiones mas recientes de ClinVar, gnomAD y HmtVar,
#       aplica el renombrado de contigs chrM → NC_012920.1 y actualiza el log de versiones.
#       Las bases de datos anteriores se sobreescriben directamente.
#
# Requisitos:
#   - Entorno conda ont_annot creado y activo:
#       conda env create -f ont_annot.yml
#       conda activate ont_annot
#       snpEff download GRCh38.99
#   - bcftools disponible en el entorno activo (verificacion automatica)
# --------------------------------------------------------------------------------------------

set -euo pipefail

#0. CONFIGURACION DE ENTORNO
TFM_DIR="/home/pdelacalvac/TFM"
DB_DIR="${TFM_DIR}/Data/Annotation_DBs"
TOOLS_DIR="${TFM_DIR}/Tools"
HAPLOGREP_DIR="${TOOLS_DIR}/haplogrep3"
SCRIPTS_DIR="${TFM_DIR}/Scripts"
LOG_DIR="${TFM_DIR}/Data/Processed/logs"
DB_VERSION_LOG="${DB_DIR}/db_versions.log"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")

# Contig mtDNA en la referencia rCRS
CONTIG="NC_012920.1"

# Enlace simbolico de HaploGrep3 en el entorno conda activo
CONDA_BIN="${CONDA_PREFIX}/bin"

# Modo de ejecucion
UPDATE_ONLY=false
if [[ "${1:-}" == "--update" ]]; then
    UPDATE_ONLY=true
fi

#0.1. FUNCIONES AUXILIARES 
log() { echo "[$(date '+%H:%M:%S')] $*"; }
err() { echo "[ERROR] $*" >&2; exit 1; }

#0.1.1. Verificacion de comandos disponibles:
require() {
    command -v "$1" &>/dev/null || err "Comando no encontrado: '$1'. ¿Está el entorno ont_annot activo?"
}

#01.1.2. Funcion para descargas mostrando progreso
download() {
    local URL="$1"
    local DEST="$2"
    log "Descargando: $(basename "${DEST}")"
    wget -q -O "${DEST}" "${URL}" \
        || err "Fallo en la descarga de: ${URL}"
}

#0.1.3. Renombra contig chrM → NC_012920.1 en un VCF.gz e indexa .tbi
    #Entrada:  archivo .vcf.gz o .vcf.bgz con chrM
    #Salida:   archivo _rCRS.vcf.gz con NC_012920.1 + indice .tbi
rename_and_index() {
    local INPUT="$1"
    local OUTPUT="$2"
    local RENAME_TMP="${DB_DIR}/chrM_to_rCRS.txt"

    echo "chrM ${CONTIG}" > "${RENAME_TMP}"
    bcftools annotate --rename-chrs "${RENAME_TMP}" "${INPUT}" \
        -Oz -o "${OUTPUT}"
    bcftools index --tbi "${OUTPUT}"
    rm -f "${RENAME_TMP}"
}

#0.1.4. Registro de fecha y version debase de datos en db_versions.log
log_db_version() {
    local DB_NAME="$1"
    local VERSION="$2"
    local SOURCE="$3"
    printf "%-20s %-30s %-50s %s\n" \
        "${DB_NAME}" "${VERSION}" "${SOURCE}" "$(date '+%Y-%m-%d %H:%M:%S')" \
        >> "${DB_VERSION_LOG}"
}

# --------------------------------------------------------------------------------------------
#1. VERIFICACIONES
# Cabecera del log
log "--------------------------------------------------------------------------------------------"
if [[ "${UPDATE_ONLY}" == true ]]; then
    log "Modo: ACTUALIZACIÓN de bases de datos (--update)"
else
    log "Modo: INSTALACIÓN COMPLETA"
fi
log "--------------------------------------------------------------------------------------------"

# Requisitos 
require bcftools
require wget
require java
require bgzip

# Crear directorios necesarios
mkdir -p "${DB_DIR}" "${TOOLS_DIR}" 

# Inicializar log de versiones si no existe
if [[ ! -f "${DB_VERSION_LOG}" ]]; then
    printf "%-20s %-30s %-50s %s\n" \
        "BASE_DE_DATOS" "VERSION" "FUENTE" "FECHA_DESCARGA" \
        > "${DB_VERSION_LOG}"
    printf "%s\n" "$(printf '%.0s-' {1..110})" >> "${DB_VERSION_LOG}"
fi

# --------------------------------------------------------------------------------------------
#2. HAPLOGREP3 (instalacion completa)

if [[ "${UPDATE_ONLY}" == false ]]; then
    log "-------------------------- Descargando Haplogrep3 -------------------------- "
    mkdir -p "${HAPLOGREP_DIR}"

    # Uso de variables para comodidad en las actualizaciones:
    HAPLOGREP_VERSION_TAG="v3.2.2"
    HAPLOGREP_ZIP="haplogrep3-3.2.2-linux.zip"
    HAPLOGREP_URL="https://github.com/genepi/haplogrep3/releases/download/${HAPLOGREP_VERSION_TAG}/${HAPLOGREP_ZIP}"
    HAPLOGREP_JAR="${HAPLOGREP_DIR}/haplogrep3"

    # Descargar, descomprimir 
    HAPLOGREP_ZIP_PATH="${HAPLOGREP_DIR}/${HAPLOGREP_ZIP}"
    download "${HAPLOGREP_URL}" "${HAPLOGREP_ZIP_PATH}"
    unzip -q "${HAPLOGREP_ZIP_PATH}" -d "${HAPLOGREP_DIR}"
    rm -f "${HAPLOGREP_ZIP_PATH}"

    # Enlace simbolico en el entorno conda para uso en el pipeline
    SYMLINK="${CONDA_BIN}/haplogrep3"
    [[ -L "${SYMLINK}" ]] && rm -f "${SYMLINK}"
    ln -s "${HAPLOGREP_JAR}.jar" "${SYMLINK}.jar"
    
    # Verificar que el jar es valido + version
    HP_VERSION=$(java -jar "${HAPLOGREP_JAR}.jar" 2>&1 | head -2)
    log_db_version "HaploGrep3" "${HP_VERSION}" "${HAPLOGREP_URL}"
    log "HaploGrep3 instalado correctamente: ${HP_VERSION}"

fi

# --------------------------------------------------------------------------------------------
#3. BASE DE DATOS ClinVar (mtDNA, GRCh38)

log " -------------------------- Descargando ClinVar -------------------------- "

# Uso de variables para comodidad en las actualizaciones:
CLINVAR_RAW="${DB_DIR}/clinvar_raw.vcf.gz"
CLINVAR_TMP="${DB_DIR}/clinvar_mt_raw.vcf.gz"
CLINVAR_OUT="${DB_DIR}/clinvar_mt_rCRS.vcf.gz"
CLINVAR_URL="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
CLINVAR_TBI_URL="${CLINVAR_URL}.tbi"

download "${CLINVAR_URL}"     "${CLINVAR_RAW}"
download "${CLINVAR_TBI_URL}" "${CLINVAR_RAW}.tbi"

# Filtrar solo variantes mitocondriales (chrM) antes de renombrar
log "Filtrando variantes mitocondriales de ClinVar..."
bcftools view -r chrM "${CLINVAR_RAW}" -Oz -o "${CLINVAR_TMP}"
bcftools index --tbi "${CLINVAR_TMP}"
 
# Renombrar chrM → NC_012920.1 e indexar
rename_and_index "${CLINVAR_TMP}" "${CLINVAR_OUT}"
 
# Registrar version (##fileDate, cabecera del VCF)
CLINVAR_DATE=$(bcftools view -h "${CLINVAR_OUT}" \
    | grep "^##fileDate" | head -1 | cut -d= -f2 || echo "desconocida")
log_db_version "ClinVar_mt" "${CLINVAR_DATE}" "${CLINVAR_URL}"

# Limpiar intermedios
rm -f "${CLINVAR_RAW}" "${CLINVAR_RAW}.tbi" "${CLINVAR_TMP}" "${CLINVAR_TMP}.tbi"
log "ClinVar procesado: ${CLINVAR_OUT} (fecha DB: ${CLINVAR_DATE})"

# --------------------------------------------------------------------------------------------
# 4. BASE DE DATOS gnomAD mitocondrial v3.1
log " -------------------------- Descargando gnomAD mitocondrial -------------------------- "

GNOMAD_RAW="${DB_DIR}/gnomad_mt_raw.vcf.bgz"
GNOMAD_OUT="${DB_DIR}/gnomad_mt_rCRS.vcf.gz"
GNOMAD_URL="https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1/vcf/genomes/gnomad.genomes.v3.1.sites.chrM.vcf.bgz"
GNOMAD_TBI_URL="${GNOMAD_URL}.tbi"

download "${GNOMAD_URL}"     "${GNOMAD_RAW}"
download "${GNOMAD_TBI_URL}" "${GNOMAD_RAW}.tbi"

log "Renombrando contig gnomAD chrM → ${CONTIG}"
rename_and_index "${GNOMAD_RAW}" "${GNOMAD_OUT}"

log_db_version "gnomAD_mt" "v3.1" "${GNOMAD_URL}"

rm -f "${GNOMAD_RAW}" "${GNOMAD_RAW}.tbi"
log "gnomAD mitocondrial descargado y procesado → ${GNOMAD_OUT}"

# --------------------------------------------------------------------------------------------
# 5. RESUMEN FINAL

{
    echo ""
    echo "--------------------------------------------------------------------------------------------"
    if [[ "${UPDATE_ONLY}" == true ]]; then
        echo "  ACTUALIZACIÓN COMPLETADA — $(date)"
    else
        echo "  INSTALACIÓN COMPLETADA — $(date)"
    fi
    echo "--------------------------------------------------------------------------------------------"
    echo ""
    echo "  Bases de datos instaladas en: ${DB_DIR}"
    printf "  %-30s %s\n" "ClinVar mtDNA:"   "${DB_DIR}/clinvar_mt_rCRS.vcf.gz"
    printf "  %-30s %s\n" "gnomAD mt v3.1:"  "${DB_DIR}/gnomad_mt_rCRS.vcf.gz"
    if [[ "${UPDATE_ONLY}" == false ]]; then
        printf "  %-30s %s\n" "HaploGrep3:"  "${HAPLOGREP_DIR}/haplogrep3.jar"
    fi
    echo ""
    echo " Log de versiones: ${DB_VERSION_LOG} "
    echo ""
    echo " Siguiente paso: "
    echo "    bash ${SCRIPTS_DIR}/annot_script4.sh"
    echo "--------------------------------------------------------------------------------------------"
} | tee -a "${LOG_DIR}/setup_${TIMESTAMP}.log"

# --------------------------------------------------------------------------------------------