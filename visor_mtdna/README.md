# Mitoview

Visor local para explorar y revisar los TSV de variantes mitocondriales.

## Iniciar en Windows

1. Instala Node.js 22 y Python 3.12.
2. Abre `start.cmd`. La primera ejecución instala las dependencias y compila la web.
3. Mitoview abre `http://127.0.0.1:8765`. Cierra la ventana de consola para detener el servicio.

Por defecto lee la carpeta `TSV results` situada un nivel por encima de esta aplicación. Puedes cambiarla desde el botón de ajustes dentro de la web. Los TSV se vuelven a comprobar cada cinco segundos.

Las notas, vistas guardadas, versiones importadas y anotaciones de ClinVar se guardan en `.local-data/mitoview.sqlite3`. Esa carpeta está excluida del control de versiones.

## Desarrollo

Instala los paquetes de Python de `requirements.txt` y los de Node con `npm install`. Ejecuta `python -m uvicorn server.main:app --host 127.0.0.1 --port 8765` y `npm run dev`; la web de desarrollo escucha en `http://127.0.0.1:5173` y envía `/api` al servicio local.

ClinVar se consulta desde el servicio local usando NCBI E-utilities y Variation Services. Si no hay internet, los TSV y las anotaciones guardadas siguen disponibles.
