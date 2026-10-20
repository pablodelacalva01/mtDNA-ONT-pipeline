from __future__ import annotations

import shutil
import subprocess
import threading
import time
import urllib.request
import webbrowser
from pathlib import Path

import uvicorn

ROOT = Path(__file__).resolve().parents[1]
URL = "http://127.0.0.1:8765"


def open_when_ready() -> None:
    deadline = time.monotonic() + 60
    while time.monotonic() < deadline:
        try:
            with urllib.request.urlopen(f"{URL}/api/health", timeout=1) as response:
                if response.status == 200:
                    webbrowser.open(URL, new=2)
                    return
        except Exception:
            time.sleep(0.4)


def main() -> None:
    npm = shutil.which("npm") or shutil.which("npm.cmd")
    if not npm:
        raise SystemExit("No se encontró Node.js y npm. Instala Node.js 22 para iniciar Mitoview.")
    subprocess.run([npm, "run", "build"], cwd=ROOT, check=True)
    threading.Thread(target=open_when_ready, daemon=True).start()
    uvicorn.run("server.main:app", host="127.0.0.1", port=8765, log_level="warning")


if __name__ == "__main__":
    main()
