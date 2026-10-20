export type Sample = {
  id: string; label: string; path: string; active: number; version_id?: string
  imported_at?: string; row_count?: number; distinct_count?: number
}

export type FilterCondition = {
  type: 'condition'; field: string; operator: string; value?: unknown
}

export type FilterGroup = {
  type: 'group'; op: 'and' | 'or' | 'not'; children: FilterNode[]
}

export type FilterNode = FilterCondition | FilterGroup

export type Observation = {
  id: string; sample_id: string; sample_label: string; sample_path: string
  version_id: string; version_sha256: string; imported_at: string
  chrom: string; pos: number; ref: string; alt: string; af_alt: number | null
  coverage: number | null; ad_ref: number | null; ad_alt: number | null
  gene: string | null; effect: string | null; impact: string | null
  hgvs_c: string | null; hgvs_p: string | null; clinvar_id: string | null
  clnsig: string | null; clndn: string | null; af_hom: number | null; af_het: number | null
  review_status: 'pendiente' | 'en_revision' | 'revisada'; note: string
  flagged: boolean; clinvar_status?: string; clinvar_checked_at?: string
  clinvar_error?: string; clinvar_records: ClinVarRecord[]; raw: Record<string, string>
}

export type ClinVarRecord = {
  uid: string; accession: string; title: string; classification: string
  review_status: string; conditions: string[]; updated: string; url: string
}

export async function api<T>(url: string, init?: RequestInit): Promise<T> {
  const response = await fetch(url, { ...init, headers: { 'Content-Type': 'application/json', ...init?.headers } })
  if (!response.ok) {
    const error = await response.json().catch(() => ({}))
    throw new Error(error.detail || `Error ${response.status}`)
  }
  return response.json() as Promise<T>
}
