# mtDNA-ONT-pipeline

Flujo de análisis bioinformático de lecturas largas obtenidas mediante tecnología Oxford Nanopore Technologies (ONT), orientado al estudio y análisis de variantes del ADN mitocondrial.

Desarrollado en inicio como parte de mi Trabajo Fin de Máster:

**«Análisis bioinformático y su implementación en la práctica clínica diaria de datos de NGS de lecturas largas mediante tecnología Nanopore»**

  Avanzado y finalizado en mi periodo predoctoral en el Hospital Universitario Central de Asturias.

**Autor:** Pablo de la Calva Castiñeira  
**Año:** 2026
---

## Descripción del proyecto

`mtDNA-ONT-pipeline` integra distintas etapas de análisis bioinformático para procesar lecturas de secuenciación Nanopore, desde el control de calidad y el alineamiento hasta la identificación y anotación de variantes.

El proyecto incluye scripts para automatizar las distintas etapas del análisis, entornos de Conda para gestionar las dependencias, resultados de ejemplo y un visor HTML para facilitar la consulta de los resultados.

## Características principales

- Control de calidad y filtrado de lecturas.
- Alineamiento de lecturas frente al genoma de referencia.
- Llamada de variantes.
- Anotación de variantes mitocondriales.
- Clasificación de haplogrupos.
- Generación y organización de resultados intermedios y finales.
- Visualización de resultados mediante una interfaz HTML.

## Estructura del repositorio que debe seguir el usuario (Buenas Prácticas)

La siguiente estructura representa la organización de los principales directorios del proyecto:

```text
mtDNA-ONT-pipeline/
├── Code/
│   ├── script.sh
│   └── env.yml
├── Data/
│   ├── Raw/
│   │   └── fastq_pass/
│   ├── Reference_Genome/
│   │   └── NC_012920.1.fasta.gz
│   ├── Annotation_DBs/
│   └── Processed/
├── Tools/
├── Results/
│   ├── Variants_Annotated/
│   ├── AnnotationReports/
│   ├── Haplogroups/
│   └── TSV/
├── visor_mtdna/
├── enviroment/
├── scripts/
├── results/
├── README.md
└── LICENSE
```

**Nota:** `Code/`, `Data/`, `Tools/` y `Results/` representan la estructura de trabajo mostrada en el esquema. Los directorios `scripts/`, `enviroment/` y `results/` corresponden a las rutas utilizadas en las instrucciones del repositorio. Adapta este árbol a la estructura definitiva que publiques en GitHub.

### Descripción de los directorios

| Directorio | Descripción |
|---|---|
| `Code/` | Scripts y archivos de configuración asociados al procesamiento. |
| `Data/Raw/` | Datos de secuenciación originales, organizados por muestra o barcode. |
| `Data/Reference_Genome/` | Genoma de referencia utilizado en el análisis. |
| `Data/Annotation_DBs/` | Bases de datos necesarias para la anotación. |
| `Data/Processed/` | Datos generados durante el procesamiento. |
| `Tools/` | Herramientas y recursos auxiliares utilizados por el proyecto. |
| `Results/` | Resultados organizados por tipo de análisis. |
| `Results/Variants_Annotated/` | Variantes anotadas. |
| `Results/AnnotationReports/` | Informes de anotación. |
| `Results/Haplogroups/` | Resultados de asignación de haplogrupos. |
| `Results/TSV/` | Tablas de variantes en formato TSV. |
| `visor_mtdna/` | Archivos del visor HTML para consultar los resultados. |
| `enviroment/` | Archivos YAML con los entornos de Conda. |
| `scripts/` | Scripts para ejecutar las distintas etapas del pipeline. |
| `results/` | Resultados generados por los scripts del pipeline. |

La disponibilidad de cada directorio y su contenido dependerá de la versión del proyecto y de las etapas ejecutadas.

## Requisitos

Para ejecutar el pipeline bioinformático:

- Linux, preferiblemente Ubuntu 20.04 o posterior.
- Conda o Miniconda.
- Java 11 o posterior, cuando lo requieran las herramientas correspondientes.
- Conexión a Internet para descargar o actualizar bases de datos y recursos externos.
- Espacio de almacenamiento suficiente para los datos de secuenciación y los resultados intermedios.

Para consultar el visor HTML en Windows, sigue las instrucciones de la sección **Visor HTML en Windows**.

## Instalación del pipeline bioinformático

### 1. Clonar el repositorio

```bash
git clone https://github.com/pablodelacalva01/mtDNA-ONT-pipeline.git
cd mtDNA-ONT-pipeline
```

### 2. Crear los entornos de Conda

Los entornos se crean a partir de los archivos YAML incluidos en `enviroment/`.

```bash
conda env create -f enviroment/ont_qc.yml
conda env create -f enviroment/ont_align.yml
conda env create -f enviroment/ont_vc.yml
conda env create -f enviroment/ont_annot.yml
```

Si ya existe alguno de los entornos, comprueba si necesita actualizarse antes de volver a crearlo.

### 3. Configurar las bases de datos de anotación

La preparación de las bases de datos debe realizarse la primera vez, siempre que el script correspondiente no se haya ejecutado previamente.

```bash
conda activate ont_annot
bash scripts/04_annot_setup.sh
```

Es necesario comprobar que el nombre y la ruta del script coinciden con los archivos presentes en la versión instalada del repositorio.

## Uso del pipeline

Antes de ejecutar los scripts, revisa las rutas de entrada y salida, los parámetros y las variables de configuración definidos al principio de cada archivo.

### Paso 1. Control de calidad y filtrado

```bash
conda activate ont_qc
bash scripts/01_qc_script.sh
```

Esta etapa realiza el control de calidad y el filtrado de las lecturas.

### Paso 2. Alineamiento

```bash
conda activate ont_align
bash scripts/02_alignment_script.sh
```

Esta etapa alinea las lecturas frente al genoma de referencia y genera los archivos necesarios para las etapas posteriores.

### Paso 3. Llamada de variantes

```bash
conda activate ont_vc
bash scripts/03_vc_script.sh
```

Esta etapa identifica las variantes a partir de los datos alineados.

### Paso 4. Anotación de variantes

```bash
conda activate ont_annot
bash scripts/04_annot_script.sh
```

Esta etapa añade información funcional y clínica a las variantes, según las herramientas y bases de datos configuradas.

Los comandos anteriores deben ajustarse si los nombres reales de los scripts difieren de los indicados.

## Visor HTML en Windows

El proyecto incluye un visor HTML para facilitar la consulta de los resultados del análisis.

### Cómo iniciar el visor

1. Localiza la carpeta que contiene el visor HTML y el archivo `start.cmd`.
2. Abre la consola de Windows (`cmd`).
3. Accede a la carpeta del visor utilizando el comando `cd`. Por ejemplo:

   ```cmd
   cd /d "C:\ruta\completa\hasta\visor_mtdna"
   ```

   Sustituye la ruta del ejemplo por la ubicación real del directorio en tu equipo.

4. Ejecuta el archivo de inicio:

   ```cmd
   start.cmd
   ```

5. Tras la ejecución deberá esperar unos segundos a que se abra un navegador automáticamente en un puerto abierto (ej.: http://127.0.0.1:8765/).

### Consideraciones

- Mantén los archivos del visor y los recursos que necesita en las ubicaciones esperadas por la aplicación.
- Si el visor utiliza archivos de resultados externos, comprueba que las rutas sean correctas.
- Si el visor no se inicia, ejecuta `start.cmd` desde `cmd` para poder leer los mensajes de error.
- No asumas que basta con abrir el archivo HTML directamente en el navegador: si el visor requiere un servidor local o un proceso de inicio, debe ejecutarse mediante el mecanismo proporcionado.
- Para visualizar los resultados finales `TSV` en el visor debe cargar la carpeta completa en la esquina superior izquierda.

## Resultados

Los resultados generados por el pipeline se organizan según la etapa de análisis y el tipo de archivo.

Entre los resultados que pueden generarse se incluyen:

- Archivos de alineamiento.
- Llamadas de variantes.
- Tablas de variantes anotadas en formato TSV.
- Informes de anotación.
- Archivos de asignación de haplogrupos.
- Resultados preparados para su visualización.

Los archivos concretos dependerán de los datos de entrada, de la configuración y de las herramientas utilizadas.

## Herramientas utilizadas

Las herramientas empleadas en el pipeline incluyen:

**Control de calidad**
- NanoPlot (v.1.42.0)
- Porechop (v.0.2.4)
- NanoFilt (v.2.8.0)
- Filtlong (v.0.2.1)

**Alineamiento y procesamiento**
- Minimap2 (v.2.28)
- NGMLR (v.0.2.7)
- Samtools (v.1.19)
- GATK4 (v.4.4.0.0)
- Qualimap (v.2.3)

**Llamada de variantes**
- Clair3 (v.2.0.1)
- Mutserve2 (v.2.0.1)
- BCFtools (v.1.23)
- Samtools (v.1.23)

**Anotación**
- SnpEff (v.5.2)
- BCFtools (v.1.23)
- HaploGrep3 (v3.2.2)
- ClinVar
- gnomAD mitochondrial (v.3.1)

La selección y disponibilidad de las herramientas dependen de los entornos y de las versiones configuradas en el proyecto.

## Reproducibilidad y limitaciones

Para facilitar la reproducibilidad de los análisis, se recomienda:

- Comprobar las versiones de las herramientas utilizadas.
- Conservar los parámetros de ejecución y la configuración de los entornos.
- Documentar el origen y la versión del genoma de referencia.
- Indicar las versiones de las bases de datos de anotación.
- Diferenciar los resultados de ejemplo de los resultados obtenidos en análisis reales.

Los resultados bioinformáticos deben interpretarse en función de la calidad de los datos, los parámetros empleados y las limitaciones de cada herramienta.

## Autoría

Proyecto desarrollado por Pablo de la Calva Castiñeira en el marco de su Trabajo Fin de Máster (2026).

## Licencia

Este proyecto se distribuye bajo la licencia MIT. Consulta el archivo [LICENSE](LICENSE) para conocer sus condiciones de uso.
