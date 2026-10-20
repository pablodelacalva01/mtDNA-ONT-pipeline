@echo off
setlocal
cd /d "%~dp0"

if not exist .venv\Scripts\python.exe (
  python -m venv --system-site-packages .venv
  if errorlevel 1 exit /b 1
  call .venv\Scripts\activate.bat
  python -m pip install -r requirements.txt
  if errorlevel 1 exit /b 1
) else (
  call .venv\Scripts\activate.bat
)

if not exist node_modules call npm install
if errorlevel 1 exit /b 1

python -m server.run

