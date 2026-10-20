from __future__ import annotations

import asyncio
import csv
import hashlib
import json
import math
import os
import re
import sqlite3
import time
import uuid
from contextlib import asynccontextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Literal

import httpx
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, Response
from pydantic import BaseModel, Field
from starlette.staticfiles import StaticFiles


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / ".local-data"
DB_PATH = DATA_DIR / "mitoview.sqlite3"
DEFAULT_TSV_DIR = ROOT.parent / "TSV results"
TSV_DIR = Path(os.environ.get("MITOVIEW_TSV_DIR", str(DEFAULT_TSV_DIR))).expanduser().resolve()
POLL_SECONDS = 5
CLINVAR_TTL_SECONDS = 7 * 24 * 60 * 60
NCBI_EMAIL = os.environ.get("NCBI_EMAIL", "")

FIELDS = {
    "sample": "s.id",
    "sampleLabel": "s.label",
    "chrom": "o.chrom",
    "pos": "o.pos",
    "ref": "o.ref",
    "alt": "o.alt",
    "change": "o.alt",
    "afAlt": "o.af_alt",
    "coverage": "o.coverage",
    "gene": "o.gene",
    "effect": "o.effect",
    "impact": "o.impact",
    "hgvsC": "o.hgvs_c",
    "hgvsP": "o.hgvs_p",
    "afHom": "o.af_hom",
    "afHet": "o.af_het",
    "clinvarStatus": "c.status",
    "reviewStatus": "COALESCE(r.status, 'pendiente')",
    "note": "COALESCE(r.note, '')",
    "flagged": "COALESCE(r.flagged, 0)",
}
NUMERIC_FIELDS = {"pos", "afAlt", "coverage", "afHom", "afHet"}
TEXT_FIELDS = set(FIELDS) - NUMERIC_FIELDS
OPERATORS = {"eq", "neq", "gt", "gte", "lt", "lte", "between", "contains", "in", "isMissing", "isPresent"}
HEADER_PREFIX = "CHROM\tPOS\t"


def now_iso() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def connect() -> sqlite3.Connection:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    db = sqlite3.connect(DB_PATH, timeout=10)
    db.row_factory = sqlite3.Row
    db.execute("PRAGMA journal_mode=WAL")
    db.execute("PRAGMA foreign_keys=ON")
    return db


def initialize() -> None:
    with connect() as db:
        db.executescript("""
        CREATE TABLE IF NOT EXISTS samples (
            id TEXT PRIMARY KEY, path TEXT NOT NULL UNIQUE, label TEXT NOT NULL,
            active INTEGER NOT NULL DEFAULT 1, current_version_id TEXT, updated_at TEXT NOT NULL
        );
        CREATE TABLE IF NOT EXISTS sample_versions (
            id TEXT PRIMARY KEY, sample_id TEXT NOT NULL REFERENCES samples(id),
            sha256 TEXT NOT NULL, imported_at TEXT NOT NULL, row_count INTEGER NOT NULL,
            active INTEGER NOT NULL DEFAULT 1
        );
        CREATE TABLE IF NOT EXISTS observations (
            id TEXT PRIMARY KEY, version_id TEXT NOT NULL REFERENCES sample_versions(id),
            chrom TEXT, pos INTEGER, ref TEXT, alt TEXT, af_alt REAL, coverage INTEGER,
            ad_ref INTEGER, ad_alt INTEGER, gene TEXT, effect TEXT, impact TEXT,
            hgvs_c TEXT, hgvs_p TEXT, clinvar_id TEXT, clnsig TEXT, clndn TEXT,
            af_hom REAL, af_het REAL, ac_hom INTEGER, ac_het INTEGER, raw_json TEXT NOT NULL
        );
        CREATE INDEX IF NOT EXISTS observations_version ON observations(version_id);
        CREATE INDEX IF NOT EXISTS observations_variant ON observations(chrom,pos,ref,alt);
        CREATE TABLE IF NOT EXISTS reviews (
            observation_id TEXT PRIMARY KEY REFERENCES observations(id), status TEXT NOT NULL,
            note TEXT NOT NULL DEFAULT '', flagged INTEGER NOT NULL DEFAULT 0, updated_at TEXT NOT NULL
        );
        CREATE TABLE IF NOT EXISTS saved_views (
            id TEXT PRIMARY KEY, name TEXT NOT NULL, definition_json TEXT NOT NULL, updated_at TEXT NOT NULL
        );
        CREATE TABLE IF NOT EXISTS clinvar_cache (
            variant_key TEXT PRIMARY KEY, status TEXT NOT NULL, canonical_spdi TEXT,
            records_json TEXT NOT NULL DEFAULT '[]', checked_at TEXT, error TEXT
        );
        CREATE TABLE IF NOT EXISTS import_issues (
            path TEXT NOT NULL, line INTEGER NOT NULL DEFAULT 0, message TEXT NOT NULL, updated_at TEXT NOT NULL,
            PRIMARY KEY(path,line,message)
        );
        CREATE TABLE IF NOT EXISTS settings (key TEXT PRIMARY KEY, value TEXT NOT NULL);
        """)
        db.execute("INSERT OR IGNORE INTO settings(key,value) VALUES('tsv_dir',?)", (str(TSV_DIR),))


def parse_optional_float(value: str | None) -> float | None:
    if value is None or value.strip() in ("", "."):
        return None
    return float(value.strip().replace(",", "."))


def parse_optional_int(value: str | None) -> int | None:
    if value is None or value.strip() in ("", "."):
        return None
    return int(value.strip())


def split_ad(value: str | None) -> tuple[int | None, int | None]:
    if value is None or value.strip() in ("", "."):
        return None, None
    parts = value.strip().split(",")
    if len(parts) != 2:
        raise ValueError("AD debe tener dos valores separados por coma")
    return parse_optional_int(parts[0]), parse_optional_int(parts[1])


def read_tsv(path: Path) -> tuple[list[dict[str, str]], list[tuple[int, str]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        lines = handle.readlines()
    header_idx = next((i for i, line in enumerate(lines) if line.rstrip("\r\n").startswith(HEADER_PREFIX)), None)
    if header_idx is None:
        raise ValueError("No se encontró la cabecera CHROM / POS")
    reader = csv.DictReader(lines[header_idx:], delimiter="\t")
    required = {"CHROM", "POS", "REF", "ALT", "AF_ALT", "COVERAGE", "AD", "GENE", "EFFECT", "IMPACT"}
    missing = required - set(reader.fieldnames or [])
    if missing:
        raise ValueError("Faltan columnas: " + ", ".join(sorted(missing)))
    rows: list[dict[str, str]] = []
    issues: list[tuple[int, str]] = []
    for line_number, row in enumerate(reader, header_idx + 2):
        if not row or not row.get("CHROM"):
            continue
        if None in row or any(value is None for value in row.values()):
            issues.append((line_number, "Fila con más columnas que la cabecera; se omitió"))
            continue
        try:
            pos = parse_optional_int(row.get("POS"))
            af = parse_optional_float(row.get("AF_ALT"))
            depth = parse_optional_int(row.get("COVERAGE"))
            ad_ref, ad_alt = split_ad(row.get("AD"))
            if pos is None or pos < 1:
                raise ValueError("POS debe ser un entero positivo")
            if af is not None and (not math.isfinite(af) or not 0 <= af <= 1):
                raise ValueError("AF_ALT debe estar entre 0 y 1")
            if depth is not None and depth < 0:
                raise ValueError("COVERAGE no puede ser negativo")
            if ad_ref is not None and ad_ref < 0 or ad_alt is not None and ad_alt < 0:
                raise ValueError("AD no puede ser negativo")
            if not row.get("REF") or not row.get("ALT"):
                raise ValueError("REF y ALT no pueden estar vacíos")
            if row.get("CHROM") != "NC_012920.1":
                raise ValueError("La referencia debe ser NC_012920.1 (rCRS)")
            rows.append({k: (v or "").strip() for k, v in row.items()})
        except (ValueError, TypeError) as error:
            issues.append((line_number, f"{error}; se omitió la fila"))
    return rows, issues


def num(value: str | None, parser: Any) -> Any:
    return parser(value)


def import_file(path: Path) -> dict[str, Any]:
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    stat = path.stat()
    rows, issues = read_tsv(path)
    full_path = str(path.resolve())
    label = path.stem.removesuffix("_annotated")
    match = re.search(r"Muestra:\s*([^\s]+)", path.read_text(encoding="utf-8-sig", errors="replace")[:1500])
    if match:
        label = match.group(1)
    with connect() as db:
        db.execute("DELETE FROM import_issues WHERE path=?", (full_path,))
        for line_number, message in issues:
            db.execute("INSERT OR REPLACE INTO import_issues VALUES(?,?,?,?)", (full_path, line_number, message, now_iso()))
        sample = db.execute("SELECT * FROM samples WHERE path=?", (full_path,)).fetchone()
        if sample is None:
            sample_id = str(uuid.uuid4())
            db.execute("INSERT INTO samples(id,path,label,active,updated_at) VALUES(?,?,?,1,?)", (sample_id, full_path, label, now_iso()))
            current_hash = None
        else:
            sample_id = sample["id"]
            current_hash = None
            if sample["current_version_id"]:
                found = db.execute("SELECT sha256 FROM sample_versions WHERE id=?", (sample["current_version_id"],)).fetchone()
                current_hash = found["sha256"] if found else None
            db.execute("UPDATE samples SET label=?,active=1,updated_at=? WHERE id=?", (label, now_iso(), sample_id))
        if current_hash == digest:
            return {"path": full_path, "status": "sin_cambios", "rows": len(rows), "issues": len(issues)}
        version_id = str(uuid.uuid4())
        db.execute("UPDATE sample_versions SET active=0 WHERE sample_id=?", (sample_id,))
        db.execute("INSERT INTO sample_versions(id,sample_id,sha256,imported_at,row_count,active) VALUES(?,?,?,?,?,1)", (version_id, sample_id, digest, now_iso(), len(rows)))
        for row in rows:
            ad_ref, ad_alt = split_ad(row.get("AD"))
            db.execute("""INSERT INTO observations VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""", (
                str(uuid.uuid4()), version_id, row.get("CHROM") or None, num(row.get("POS"), parse_optional_int),
                row.get("REF") or None, row.get("ALT") or None, num(row.get("AF_ALT"), parse_optional_float),
                num(row.get("COVERAGE"), parse_optional_int), ad_ref, ad_alt,
                row.get("GENE") or None, row.get("EFFECT") or None, row.get("IMPACT") or None,
                row.get("HGVS_C") or None, row.get("HGVS_P") or None,
                row.get("CLINVAR_ID") or None, row.get("CLNSIG") or None, row.get("CLNDN") or None,
                num(row.get("AF_hom"), parse_optional_float), num(row.get("AF_het"), parse_optional_float),
                num(row.get("AC_hom"), parse_optional_int), num(row.get("AC_het"), parse_optional_int),
                json.dumps(row, ensure_ascii=False),
            ))
        db.execute("UPDATE samples SET current_version_id=?,active=1,updated_at=? WHERE id=?", (version_id, now_iso(), sample_id))
    return {"path": full_path, "status": "importado", "rows": len(rows), "issues": len(issues)}


def current_tsv_dir() -> Path:
    with connect() as db:
        setting = db.execute("SELECT value FROM settings WHERE key='tsv_dir'").fetchone()
    return Path(setting["value"] if setting else str(TSV_DIR)).expanduser().resolve()


sync_lock = asyncio.Lock()


async def sync_folder() -> dict[str, Any]:
    async with sync_lock:
        folder = current_tsv_dir()
        if not folder.is_dir():
            return {"folder": str(folder), "files": [], "error": "La carpeta configurada no existe o no está disponible."}
        paths = sorted(p for p in folder.glob("*.tsv") if p.is_file())
        normalized: list[dict[str, Any]] = []
        for path in paths:
            try:
                result = await asyncio.to_thread(import_file, path)
            except Exception as error:
                with connect() as db:
                    db.execute("INSERT OR REPLACE INTO import_issues VALUES(?,?,?,?)", (str(path.resolve()), 0, str(error), now_iso()))
                normalized.append({"path": str(path), "status": "error", "error": str(error)})
            else:
                normalized.append(result)
        present = {str(path.resolve()) for path in paths}
        with connect() as db:
            for row in db.execute("SELECT id,path FROM samples").fetchall():
                if row["path"] not in present:
                    db.execute("UPDATE samples SET active=0,updated_at=? WHERE id=?", (now_iso(), row["id"]))
        return {"folder": str(folder), "files": normalized, "errors": sum(x.get("status") == "error" for x in normalized)}


def compile_filter(node: dict[str, Any] | None, params: list[Any]) -> str:
    if not node:
        return "1=1"
    if node.get("type") == "group":
        op = node.get("op")
        children = node.get("children", [])
        if op not in ("and", "or", "not") or not isinstance(children, list) or len(children) > 20:
            raise HTTPException(422, "Grupo de filtros no válido")
        if op == "not":
            if len(children) != 1:
                raise HTTPException(422, "El grupo NOT debe tener una condición")
            return f"NOT ({compile_filter(children[0], params)})"
        glue = " AND " if op == "and" else " OR "
        return "(" + glue.join(compile_filter(child, params) for child in children) + ")" if children else "1=1"
    field = node.get("field")
    operator = node.get("operator")
    if field not in FIELDS or operator not in OPERATORS:
        raise HTTPException(422, "Campo u operador de filtro no permitido")
    sql_field = FIELDS[field]
    value = node.get("value")
    if operator == "isMissing":
        return f"({sql_field} IS NULL OR {sql_field} = '')"
    if operator == "isPresent":
        return f"({sql_field} IS NOT NULL AND {sql_field} != '')"
    if operator == "in":
        if not isinstance(value, list) or len(value) > 100:
            raise HTTPException(422, "El filtro IN requiere una lista de hasta 100 valores")
        if not value:
            return "0=1"
        params.extend(value)
        return f"{sql_field} IN ({','.join('?' for _ in value)})"
    if operator == "between":
        if not isinstance(value, list) or len(value) != 2:
            raise HTTPException(422, "El intervalo requiere dos extremos")
        params.extend(value)
        return f"{sql_field} BETWEEN ? AND ?"
    if operator == "contains":
        params.append(f"%{str(value)}%")
        return f"LOWER(COALESCE({sql_field},'')) LIKE LOWER(?) ESCAPE '\\'"
    sql_op = {"eq": "=", "neq": "!=", "gt": ">", "gte": ">=", "lt": "<", "lte": "<="}.get(operator)
    if sql_op is None:
        raise HTTPException(422, "Operador no válido")
    if field in NUMERIC_FIELDS:
        try:
            value = float(value)
        except (TypeError, ValueError):
            raise HTTPException(422, "Este filtro necesita un valor numérico")
    params.append(value)
    return f"{sql_field} {sql_op} ?"


BASE_SELECT = """SELECT o.*, s.id AS sample_id, s.label AS sample_label, s.path AS sample_path,
       v.id AS version_id, v.sha256 AS version_sha256, v.imported_at,
       COALESCE(r.status,'pendiente') AS review_status, COALESCE(r.note,'') AS note,
       COALESCE(r.flagged,0) AS flagged, r.updated_at AS review_updated_at,
       c.status AS clinvar_status, c.canonical_spdi, c.records_json AS clinvar_records,
       c.checked_at AS clinvar_checked_at, c.error AS clinvar_error
       FROM observations o JOIN sample_versions v ON v.id=o.version_id
       JOIN samples s ON s.id=v.sample_id LEFT JOIN reviews r ON r.observation_id=o.id
       LEFT JOIN clinvar_cache c ON c.variant_key=(o.chrom || ':' || o.pos || ':' || o.ref || ':' || o.alt)"""
BASE_COUNT = """FROM observations o JOIN sample_versions v ON v.id=o.version_id
       JOIN samples s ON s.id=v.sample_id LEFT JOIN reviews r ON r.observation_id=o.id
       LEFT JOIN clinvar_cache c ON c.variant_key=(o.chrom || ':' || o.pos || ':' || o.ref || ':' || o.alt)"""


def serialize(row: sqlite3.Row) -> dict[str, Any]:
    item = dict(row)
    raw = json.loads(item.pop("raw_json", "{}"))
    records_json = item.pop("clinvar_records", "[]")
    item["raw"] = raw
    item["clinvar_records"] = normalize_clinvar_records(json.loads(records_json or "[]"))
    item["flagged"] = bool(item.get("flagged"))
    return item


def normalize_clinvar_records(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Keep records cached by earlier builds and new ESummary records in one UI shape."""
    normalized = []
    for record in records:
        raw = record.get("raw") or {}
        clinical = raw.get("germline_classification") or raw.get("clinical_significance") or {}
        if isinstance(clinical, str):
            clinical = {"description": clinical}
        if not isinstance(clinical, dict):
            clinical = {}
        traits = clinical.get("trait_set") or raw.get("trait_set") or raw.get("traits") or []
        if isinstance(traits, str):
            conditions = [traits]
        else:
            conditions = [x.get("trait_name") or x.get("name") or str(x) for x in traits if isinstance(x, (dict, str))]
        normalized.append({
            **record,
            "classification": clinical.get("description") or record.get("classification") or "",
            "review_status": clinical.get("review_status") or record.get("review_status") or "",
            "conditions": conditions or record.get("conditions") or [],
        })
    return normalized


def query_rows(filter_tree: dict[str, Any] | None, order_by: str = "pos", direction: str = "asc", limit: int | None = 5000) -> tuple[list[dict[str, Any]], int]:
    params: list[Any] = []
    where = compile_filter(filter_tree, params)
    order_field = FIELDS.get(order_by, "o.pos")
    direction_sql = "DESC" if direction.lower() == "desc" else "ASC"
    with connect() as db:
        count = db.execute(f"SELECT COUNT(*) {BASE_COUNT} WHERE v.active=1 AND s.active=1 AND {where}", params).fetchone()[0]
        query = f"{BASE_SELECT} WHERE v.active=1 AND s.active=1 AND {where} ORDER BY {order_field} {direction_sql}, o.pos ASC"
        rows = db.execute(query, params).fetchall() if limit is None else db.execute(query + " LIMIT ?", [*params, limit]).fetchall()
    return [serialize(row) for row in rows], count


class SearchRequest(BaseModel):
    filter: dict[str, Any] | None = None
    page: int = Field(default=1, ge=1)
    pageSize: int = Field(default=50, ge=1, le=500)
    sortBy: str = "pos"
    sortDirection: Literal["asc", "desc"] = "asc"


class ReviewRequest(BaseModel):
    status: Literal["pendiente", "en_revision", "revisada"]
    note: str = Field(default="", max_length=10000)
    flagged: bool = False


class ViewRequest(BaseModel):
    id: str | None = None
    name: str = Field(min_length=1, max_length=80)
    definition: dict[str, Any]


class FolderRequest(BaseModel):
    path: str = Field(min_length=1, max_length=1000)


class FilterRequest(BaseModel):
    filter: dict[str, Any] | None = None


class CompareRequest(BaseModel):
    sampleIds: list[str] = Field(min_length=1, max_length=12)
    filter: dict[str, Any] | None = None


request_lock = asyncio.Lock()
last_ncbi_request = 0.0


async def ncbi_json(client: httpx.AsyncClient, url: str, params: dict[str, str] | None = None) -> dict[str, Any]:
    global last_ncbi_request
    if NCBI_EMAIL and params is not None:
        params = {**params, "email": NCBI_EMAIL}
    async with request_lock:
        delay = 0.52 - (time.monotonic() - last_ncbi_request)
        if delay > 0:
            await asyncio.sleep(delay)
        response = await client.get(url, params=params, timeout=20)
        last_ncbi_request = time.monotonic()
        response.raise_for_status()
        return response.json()


async def fetch_clinvar_record(client: httpx.AsyncClient, canonical: str) -> tuple[str, list[dict[str, Any]]]:
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    search = await ncbi_json(client, f"{base}/esearch.fcgi", {"db": "clinvar", "term": f'"{canonical}"[cspdi]', "retmode": "json", "retmax": "20"})
    ids = search.get("esearchresult", {}).get("idlist", [])
    if not ids:
        return "sin_coincidencia", []
    summary = await ncbi_json(client, f"{base}/esummary.fcgi", {"db": "clinvar", "id": ",".join(ids), "retmode": "json"})
    result = summary.get("result", {})
    records: list[dict[str, Any]] = []
    for uid in result.get("uids", ids):
        rec = result.get(uid, {})
        clinical = rec.get("germline_classification") or rec.get("clinical_significance") or {}
        if isinstance(clinical, str):
            clinical = {"description": clinical}
        if not isinstance(clinical, dict):
            clinical = {}
        traits = clinical.get("trait_set") or rec.get("trait_set") or rec.get("traits") or []
        if isinstance(traits, str):
            conditions = [traits]
        else:
            conditions = [x.get("trait_name") or x.get("name") or str(x) for x in traits if isinstance(x, (dict, str))]
        records.append({
            "uid": uid,
            "accession": rec.get("accession") or rec.get("accession_version") or rec.get("variation_set", [{}])[0].get("variation_name"),
            "title": rec.get("title", ""),
            "classification": clinical.get("description") or "",
            "review_status": clinical.get("review_status") or rec.get("review_status") or "",
            "conditions": conditions,
            "updated": rec.get("date_last_updated") or rec.get("date_created") or "",
            "url": f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{uid}/",
            "raw": rec,
        })
    return ("coincidencia_verificada" if len(records) == 1 else "ambigua"), records


def vcf_spdi(chrom: str, pos: int, ref: str, alt: str) -> str:
    ref, alt = ref.upper(), alt.upper()
    start = 0
    while start < min(len(ref), len(alt)) and ref[start] == alt[start]:
        start += 1
    end = 0
    while end < min(len(ref) - start, len(alt) - start) and ref[len(ref) - 1 - end] == alt[len(alt) - 1 - end]:
        end += 1
    deleted = ref[start:len(ref) - end if end else len(ref)]
    inserted = alt[start:len(alt) - end if end else len(alt)]
    return f"{chrom}:{pos - 1 + start}:{deleted}:{inserted}"


async def annotate_missing() -> None:
    client = httpx.AsyncClient(headers={"User-Agent": "Mitoview-local/0.1 (ClinVar annotation)"})
    try:
        while True:
            with connect() as db:
                variants = db.execute("""SELECT DISTINCT o.chrom,o.pos,o.ref,o.alt
                    FROM observations o JOIN sample_versions v ON v.id=o.version_id AND v.active=1
                    WHERE o.chrom IS NOT NULL AND o.pos IS NOT NULL AND o.ref IS NOT NULL AND o.alt IS NOT NULL""").fetchall()
                cache = {r["variant_key"]: r["checked_at"] for r in db.execute("SELECT variant_key,checked_at FROM clinvar_cache")}
            for variant in variants:
                key = f"{variant['chrom']}:{variant['pos']}:{variant['ref']}:{variant['alt']}"
                checked = cache.get(key)
                if checked:
                    try:
                        age = (datetime.now(timezone.utc) - datetime.fromisoformat(checked)).total_seconds()
                        if age < CLINVAR_TTL_SECONDS:
                            continue
                    except ValueError:
                        pass
                try:
                    spdi = vcf_spdi(variant["chrom"], variant["pos"], variant["ref"], variant["alt"])
                    variation_url = "https://api.ncbi.nlm.nih.gov/variation/v0/spdi/" + httpx.URL("https://x/" + spdi).path + "/canonical_representative"
                    norm = await ncbi_json(client, variation_url)
                    canonical_data = norm.get("data", {}).get("spdi", norm.get("data", {}))
                    canonical = ":".join(str(canonical_data[k]) for k in ("seq_id", "position", "deleted_sequence", "inserted_sequence"))
                    status, records = await fetch_clinvar_record(client, canonical)
                    with connect() as db:
                        db.execute("INSERT OR REPLACE INTO clinvar_cache VALUES(?,?,?,?,?,NULL)", (key, status, canonical, json.dumps(records, ensure_ascii=False), now_iso()))
                except Exception as error:
                    with connect() as db:
                        db.execute("INSERT OR REPLACE INTO clinvar_cache VALUES(?, 'error', NULL, '[]', ?, ?)", (key, now_iso(), str(error)[:1000]))
                await asyncio.sleep(0.05)
            try:
                await asyncio.wait_for(clinvar_refresh_event.wait(), timeout=3600)
                clinvar_refresh_event.clear()
            except TimeoutError:
                pass
    finally:
        await client.aclose()


background_tasks: list[asyncio.Task[Any]] = []
clinvar_refresh_event = asyncio.Event()


async def sync_loop() -> None:
    while True:
        await sync_folder()
        await asyncio.sleep(POLL_SECONDS)


@asynccontextmanager
async def lifespan(_: FastAPI):
    initialize()
    await sync_folder()
    background_tasks.extend([asyncio.create_task(sync_loop()), asyncio.create_task(annotate_missing())])
    yield
    for task in background_tasks:
        task.cancel()


app = FastAPI(title="Mitoview API", version="0.1.0", lifespan=lifespan)
app.add_middleware(CORSMiddleware, allow_origins=["http://127.0.0.1:5173", "http://localhost:5173"], allow_methods=["GET", "POST", "PUT", "DELETE"], allow_headers=["Content-Type"])


@app.get("/api/health")
def health() -> dict[str, str]:
    return {"status": "ok", "folder": str(current_tsv_dir())}


@app.get("/api/samples")
def samples() -> dict[str, Any]:
    with connect() as db:
        records = db.execute("""SELECT s.id,s.label,s.path,s.active,v.id AS version_id,v.sha256,v.imported_at,v.row_count,
            (SELECT COUNT(DISTINCT chrom || ':' || pos || ':' || ref || ':' || alt) FROM observations WHERE version_id=v.id) AS distinct_count
            FROM samples s LEFT JOIN sample_versions v ON v.id=s.current_version_id ORDER BY s.label""").fetchall()
        issues = db.execute("SELECT path,line,message,updated_at FROM import_issues ORDER BY path,line").fetchall()
    return {"samples": [dict(r) for r in records], "issues": [dict(r) for r in issues], "folder": str(current_tsv_dir())}


@app.get("/api/summary")
def summary() -> dict[str, Any]:
    with connect() as db:
        row = db.execute("""SELECT COUNT(*) AS observations,COUNT(DISTINCT o.chrom||':'||o.pos||':'||o.ref||':'||o.alt) AS variants,
            COUNT(DISTINCT s.id) AS samples FROM observations o JOIN sample_versions v ON v.id=o.version_id AND v.active=1
            JOIN samples s ON s.id=v.sample_id AND s.active=1""").fetchone()
        pending = db.execute("""SELECT COUNT(*) FROM observations o JOIN sample_versions v ON v.id=o.version_id AND v.active=1
            JOIN samples s ON s.id=v.sample_id AND s.active=1 LEFT JOIN reviews r ON r.observation_id=o.id
            WHERE r.status IS NULL OR r.status='pendiente'""").fetchone()[0]
        statuses = db.execute("SELECT status,COUNT(*) AS count FROM clinvar_cache GROUP BY status").fetchall()
    return {**dict(row), "pending": pending, "clinvar": {r["status"]: r["count"] for r in statuses}}


@app.post("/api/observations/search")
def search_observations(request: SearchRequest) -> dict[str, Any]:
    start = (request.page - 1) * request.pageSize
    params: list[Any] = []
    where = compile_filter(request.filter, params)
    order_field = FIELDS.get(request.sortBy, "o.pos")
    direction = "DESC" if request.sortDirection == "desc" else "ASC"
    with connect() as db:
        total = db.execute(f"SELECT COUNT(*) {BASE_COUNT} WHERE v.active=1 AND s.active=1 AND {where}",params).fetchone()[0]
        rows = db.execute(f"{BASE_SELECT} WHERE v.active=1 AND s.active=1 AND {where} ORDER BY {order_field} {direction},o.pos ASC LIMIT ? OFFSET ?",[*params,request.pageSize,start]).fetchall()
    return {"rows": [serialize(row) for row in rows], "total": total, "page": request.page, "pageSize": request.pageSize}


@app.get("/api/observations/{observation_id}")
def observation_detail(observation_id: str) -> dict[str, Any]:
    with connect() as db:
        row = db.execute(f"{BASE_SELECT} WHERE o.id=?", (observation_id,)).fetchone()
    if not row:
        raise HTTPException(404, "No se encontró la variante")
    return serialize(row)


@app.put("/api/reviews/{observation_id}")
def save_review(observation_id: str, request: ReviewRequest) -> dict[str, Any]:
    with connect() as db:
        exists = db.execute("SELECT 1 FROM observations WHERE id=?", (observation_id,)).fetchone()
        if not exists:
            raise HTTPException(404, "No se encontró la variante")
        db.execute("INSERT OR REPLACE INTO reviews VALUES(?,?,?,?,?)", (observation_id, request.status, request.note, int(request.flagged), now_iso()))
        row = db.execute("SELECT * FROM reviews WHERE observation_id=?", (observation_id,)).fetchone()
    return dict(row)


@app.get("/api/views")
def list_views() -> list[dict[str, Any]]:
    with connect() as db:
        rows = db.execute("SELECT * FROM saved_views ORDER BY name COLLATE NOCASE").fetchall()
    return [{"id": r["id"], "name": r["name"], "definition": json.loads(r["definition_json"]), "updated_at": r["updated_at"]} for r in rows]


@app.post("/api/views")
def upsert_view(request: ViewRequest) -> dict[str, Any]:
    view_id = request.id or str(uuid.uuid4())
    with connect() as db:
        db.execute("INSERT OR REPLACE INTO saved_views VALUES(?,?,?,?)", (view_id, request.name, json.dumps(request.definition), now_iso()))
    return {"id": view_id, "name": request.name, "definition": request.definition}


@app.delete("/api/views/{view_id}")
def delete_view(view_id: str) -> dict[str, bool]:
    with connect() as db:
        db.execute("DELETE FROM saved_views WHERE id=?", (view_id,))
    return {"deleted": True}


@app.post("/api/sync")
async def manual_sync() -> dict[str, Any]:
    return await sync_folder()


@app.put("/api/settings/folder")
async def set_folder(request: FolderRequest) -> dict[str, Any]:
    path = Path(request.path).expanduser().resolve()
    if not path.is_dir():
        raise HTTPException(422, "La carpeta no existe o no está disponible")
    with connect() as db:
        db.execute("INSERT OR REPLACE INTO settings VALUES('tsv_dir',?)", (str(path),))
    return await sync_folder()


@app.get("/api/clinvar/status")
def clinvar_status() -> dict[str, Any]:
    with connect() as db:
        rows = db.execute("SELECT status,COUNT(*) AS count,MAX(checked_at) AS last_checked FROM clinvar_cache GROUP BY status").fetchall()
        total = db.execute("SELECT COUNT(DISTINCT chrom||':'||pos||':'||ref||':'||alt) FROM observations o JOIN sample_versions v ON v.id=o.version_id AND v.active=1").fetchone()[0]
    return {"total": total, "states": [dict(r) for r in rows]}


@app.post("/api/clinvar/refresh")
async def refresh_clinvar() -> dict[str, str]:
    with connect() as db:
        db.execute("DELETE FROM clinvar_cache WHERE status IN ('sin_coincidencia','coincidencia_verificada','ambigua','error')")
    clinvar_refresh_event.set()
    return {"status": "encolado"}


@app.post("/api/charts")
def chart_data(request: FilterRequest) -> dict[str, Any]:
    rows, _ = query_rows(request.filter)
    frequencies = [r["af_alt"] for r in rows if r["af_alt"] is not None]
    depths = [r["coverage"] for r in rows if r["coverage"] is not None]
    effect_counts: dict[str, int] = {}
    for row in rows:
        effect_counts[row["effect"] or "Sin anotación"] = effect_counts.get(row["effect"] or "Sin anotación", 0) + 1
    points = [{"af": r["af_alt"] * 100, "coverage": r["coverage"], "id": r["id"], "pos": r["pos"], "gene": r["gene"]} for r in rows if r["af_alt"] is not None and r["coverage"] is not None]
    return {"points": points, "effects": [{"name": k, "count": v} for k, v in effect_counts.items()], "count": len(rows), "afRange": [min(frequencies) * 100, max(frequencies) * 100] if frequencies else None, "depthRange": [min(depths), max(depths)] if depths else None}


@app.post("/api/compare")
def compare(request: CompareRequest) -> dict[str, Any]:
    if len(set(request.sampleIds)) != len(request.sampleIds):
        raise HTTPException(422, "Hay muestras repetidas")
    params: list[Any] = []
    where = compile_filter(request.filter, params)
    placeholders = ",".join("?" for _ in request.sampleIds)
    params.extend(request.sampleIds)
    with connect() as db:
        rows = db.execute(f"""{BASE_SELECT} WHERE v.active=1 AND s.active=1 AND s.id IN ({placeholders}) AND {where}
            ORDER BY o.pos,s.label""", params).fetchall()
        sample_rows = db.execute(f"SELECT id,label FROM samples WHERE active=1 AND id IN ({placeholders})", request.sampleIds).fetchall()
    items = [serialize(r) for r in rows]
    keys = sorted({f"{r['chrom']}:{r['pos']}:{r['ref']}:{r['alt']}" for r in items}, key=lambda v: int(v.split(':')[1]))
    samples_by_id = {row["id"]: row["label"] for row in sample_rows}
    cells = {key: {sample_id: None for sample_id in request.sampleIds} for key in keys}
    for row in items:
        key = f"{row['chrom']}:{row['pos']}:{row['ref']}:{row['alt']}"
        cells[key][row["sample_id"]] = {"af_alt": row["af_alt"], "coverage": row["coverage"], "observation_id": row["id"]}
    return {"samples": [{"id": sid, "label": samples_by_id.get(sid, sid)} for sid in request.sampleIds], "variants": [{"key": key, "cells": cells[key]} for key in keys]}


@app.post("/api/export")
def export_rows(request: FilterRequest) -> Response:
    rows, total = query_rows(request.filter, limit=None)
    out = []
    for row in rows:
        out.append({**row["raw"], "SAMPLE": row["sample_label"], "REVIEW_STATUS": row["review_status"], "REVIEW_NOTE": row["note"], "FLAGGED": row["flagged"], "CLINVAR_STATUS": row.get("clinvar_status") or "pendiente", "CLINVAR_RECORDS": json.dumps(row["clinvar_records"], ensure_ascii=False)})
    import io
    buffer = io.StringIO(newline="")
    fields = list(out[0].keys()) if out else ["SAMPLE", "CHROM", "POS", "REF", "ALT", "AF_ALT", "COVERAGE", "REVIEW_STATUS", "REVIEW_NOTE", "CLINVAR_STATUS"]
    writer = csv.DictWriter(buffer, fieldnames=fields, delimiter="\t", extrasaction="ignore")
    writer.writeheader()
    writer.writerows(out)
    context = {
        "generatedAt": now_iso(),
        "filter": request.filter,
        "totalObservations": total,
        "samples": sorted({r["sample_label"] for r in rows}),
        "sourceVersions": sorted({f"{r['sample_label']} · {r['version_sha256']} · {r['imported_at']}" for r in rows}),
        "clinvarSource": "NCBI ClinVar; query timestamps and record URLs are included per row",
    }
    import io
    import zipfile
    archive = io.BytesIO()
    with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED) as zipped:
        zipped.writestr("variantes.tsv", "\ufeff" + buffer.getvalue())
        zipped.writestr("contexto.json", json.dumps(context, ensure_ascii=False, indent=2))
    return Response(archive.getvalue(), media_type="application/zip", headers={"Content-Disposition": "attachment; filename=mitoview-exportacion.zip"})


@app.get("/{path:path}")
def frontend(path: str) -> Any:
    dist = ROOT / "dist"
    target = (dist / path).resolve()
    if target.is_relative_to(dist.resolve()) and target.is_file():
        return FileResponse(target)
    index = dist / "index.html"
    if index.exists():
        return FileResponse(index)
    raise HTTPException(404, "Frontend todavía no compilado")
