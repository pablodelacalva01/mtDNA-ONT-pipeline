import { useEffect, useMemo, useState } from 'react'
import { useMutation, useQuery, useQueryClient } from '@tanstack/react-query'
import { createColumnHelper, flexRender, getCoreRowModel, useReactTable, type VisibilityState } from '@tanstack/react-table'
import { Activity, ArrowDownUp, Bookmark, Check, ChevronDown, ChevronLeft, ChevronRight, CircleHelp, Database, Download, ExternalLink, Filter, Flag, FolderSync, GitCompareArrows, HeartPulse, MoreHorizontal, PanelRightClose, Plus, RefreshCw, Search, Settings2, SlidersHorizontal, X, FolderOpen } from 'lucide-react'
import { CartesianGrid, Cell, Legend, ResponsiveContainer, Scatter, ScatterChart, Tooltip, XAxis, YAxis, Bar, BarChart } from 'recharts'
import { api, type FilterNode, type Observation, type Sample } from './types'

type Tab = 'variantes' | 'comparar' | 'revisión'
type Issue = { path: string; line: number; message: string; updated_at: string }
const col = createColumnHelper<Observation>()
const columns = [
  col.accessor('pos', { header: 'Posición', cell: i => <span className="mono position">{i.getValue()}</span> }),
  col.accessor(row => `${row.ref} → ${row.alt}`, { id: 'change', header: 'Cambio', cell: i => <span className="mono allele">{i.getValue()}</span> }),
  col.accessor('gene', { header: 'Gen', cell: i => <span className="gene">{i.getValue() || '—'}</span> }),
  col.accessor('effect', { header: 'Efecto', cell: i => <span className="effect-text">{i.getValue()?.replaceAll('_', ' ') || 'Sin anotación'}</span> }),
  col.accessor('af_alt', { header: 'Fracción alélica', cell: i => i.getValue() == null ? '—' : <span className="mono">{(i.getValue()! * 100).toFixed(1)}%</span> }),
  col.accessor('coverage', { header: 'Profundidad', cell: i => i.getValue() == null ? '—' : <span className="mono">{i.getValue()}×</span> }),
  col.accessor('impact', { header: 'Impacto', cell: i => <span className={`impact ${i.getValue()?.toLowerCase()}`}>{i.getValue() || '—'}</span> }),
  col.accessor('review_status', { header: 'Revisión', cell: i => <span className={`review-pill ${i.getValue()}`}>{i.getValue() === 'revisada' ? 'Revisada' : i.getValue() === 'en_revision' ? 'En revisión' : 'Pendiente'}</span> }),
  ...([
    ['CLINVAR_ID', 'ClinVar ID'], ['CLNSIG', 'Clasificación ClinVar'], ['CLNDN', 'Enfermedad ClinVar'],
    ['AC_hom', 'gnomAD AC hom'], ['AC_het', 'gnomAD AC het'], ['AD', 'AD (ref, alt)'],
    ['VAR_TYPE', 'Tipo de variante'], ['N_SUPPORT', 'N.º combinaciones'], ['CALLED_BY', 'Detectada por'],
    ['CROSS_CALLER', 'Concordancia callers'], ['AF_MINIMAP2', 'AF minimap2'], ['DP_MINIMAP2', 'DP minimap2'],
    ['AD_MINIMAP2', 'AD minimap2'], ['AF_NGMLR', 'AF NGMLR'], ['DP_NGMLR', 'DP NGMLR'],
    ['AD_NGMLR', 'AD NGMLR'], ['AF_SOURCE', 'Fuente AF']
  ] as const).map(([key, label]) => col.accessor(row => row.raw?.[key] ?? '', {
    id: `tsv_${key}`, header: label,
    cell: i => <span className="mono">{String(i.getValue() || '—')}</span>
  })),
]

const fieldOptions = [
  ['gene', 'Gen'], ['pos', 'Posición'], ['afAlt', 'Fracción alélica'], ['coverage', 'Profundidad'],
  ['effect', 'Efecto'], ['impact', 'Impacto'], ['afHom', 'gnomAD homoplasmia'],
  ['afHet', 'gnomAD heteroplasmia'], ['clinvarStatus', 'Estado ClinVar'], ['reviewStatus', 'Estado de revisión'],
  ['note', 'Nota de revisión'], ['flagged', 'Marcada para revisar'],
] as const
const numFields = new Set(['pos', 'afAlt', 'coverage', 'afHom', 'afHet'])

function App() {
  const queryClient = useQueryClient()
  const [tab, setTab] = useState<Tab>('variantes')
  const [selectedSamples, setSelectedSamples] = useState<string[]>([])
  const [selected, setSelected] = useState<Observation | null>(null)
  const [page, setPage] = useState(1)
  const [sortBy, setSortBy] = useState('pos')
  const [sortDirection, setSortDirection] = useState<'asc'|'desc'>('asc')
  const [columnVisibility, setColumnVisibility] = useState<VisibilityState>({})
  const [search, setSearch] = useState('')
  const [gene, setGene] = useState('')
  const [effect, setEffect] = useState('')
  const [review, setReview] = useState('')
  const [minAf, setMinAf] = useState('')
  const [minDepth, setMinDepth] = useState('')
  const [advancedOpen, setAdvancedOpen] = useState(false)
  const [advancedOp, setAdvancedOp] = useState<'and' | 'or' | 'not'>('and')
  const [rules, setRules] = useState<{field:string; operator:string; value:string; value2:string}[]>([])
  const [viewId, setViewId] = useState('')
  const [compareIds, setCompareIds] = useState<string[] | null>(null)
  const [loadedFilter, setLoadedFilter] = useState<FilterNode | null>(null)
  const [folderInput, setFolderInput] = useState('')
  const [settingsOpen, setSettingsOpen] = useState(false)
  const [showCharts, setShowCharts] = useState(false)
  const [chartTab, setChartTab] = useState<'alineamiento'|'variantes'>('variantes')
  const [coverageFiles, setCoverageFiles] = useState<File[]>([])
  const [haplogroupFiles, setHaplogroupFiles] = useState<File[]>([])
  const pageSize = 25

  const samplesQuery = useQuery({ queryKey: ['samples'], queryFn: () => api<{samples: Sample[]; issues: Issue[]; folder: string}>('/api/samples'), refetchInterval: 5000 })
  const samples = samplesQuery.data?.samples.filter(s => s.active && s.version_id) ?? []
  useEffect(()=>{if(!selectedSamples.length&&samples.length)setSelectedSamples([samples[0].id])},[samples.length])
  const issues = samplesQuery.data?.issues ?? []
  const effectiveSamples = selectedSamples.length ? selectedSamples : samples.slice(0,1).map(s => s.id)
  const filter = useMemo<FilterNode>(() => {
    const children: FilterNode[] = []
    if (effectiveSamples.length && effectiveSamples.length < samples.length) children.push({type:'condition', field:'sample', operator:'in', value: effectiveSamples})
    if (search.trim()) children.push({type:'group', op:'or', children:['pos','sampleLabel','gene','hgvsC','hgvsP','ref','alt'].map(field=>({type:'condition',field,operator:'contains',value:search.trim()} as FilterNode))})
    if (gene) children.push({type:'condition',field:'gene',operator:'eq',value:gene})
    if (effect) children.push({type:'condition',field:'effect',operator:'eq',value:effect})
    if (review) children.push({type:'condition',field:'reviewStatus',operator:'eq',value:review})
    if (minAf) children.push({type:'condition',field:'afAlt',operator:'gte',value:Number(minAf)/100})
    if (minDepth) children.push({type:'condition',field:'coverage',operator:'gte',value:Number(minDepth)})
    const advanced = rules.filter(r=>r.field && r.operator && (r.operator === 'isMissing' || r.operator === 'isPresent' || (r.value !== '' && (r.operator!=='between'||r.value2!==''))))
      .map(r=>{const convert=(v:string)=>r.field==='afAlt'?Number(v)/100:numValue(v);return {type:'condition',field:r.field,operator:r.operator,value:r.operator==='between'?[convert(r.value),convert(r.value2)]:convert(r.value)} as FilterNode})
    if (advanced.length) children.push(advancedOp === 'not' ? {type:'group',op:'not',children:[{type:'group',op:'and',children:advanced}]} : {type:'group',op:advancedOp,children:advanced})
    if (loadedFilter) children.push(loadedFilter)
    return {type:'group',op:'and',children}
  }, [effectiveSamples.join(','), samples, search, gene, effect, review, minAf, minDepth, rules, advancedOp, loadedFilter])
  const filterKey = JSON.stringify(filter)
  const resultsQuery = useQuery({
    queryKey: ['observations', filterKey, page, sortBy, sortDirection],
    queryFn: () => api<{rows: Observation[]; total: number}>('/api/observations/search', {method:'POST', body:JSON.stringify({filter,page,pageSize,sortBy,sortDirection})}),
    enabled: samples.length > 0,
    refetchInterval: 5000,
  })
  const chartsQuery = useQuery({ queryKey:['charts',filterKey], queryFn:()=>api<{points:{af:number;coverage:number;id:string;pos:number;gene:string}[];effects:{name:string;count:number}[];count:number}>('/api/charts',{method:'POST',body:JSON.stringify({filter})}), enabled:samples.length>0 && showCharts })
  const viewsQuery = useQuery({queryKey:['views'],queryFn:()=>api<{id:string;name:string;definition:Record<string,unknown>}[]>('/api/views')})
  const summaryQuery = useQuery({queryKey:['summary'],queryFn:()=>api<{observations:number;variants:number;samples:number;pending:number;clinvar:Record<string,number>}>('/api/summary'),refetchInterval:10000})
  const refreshMutation = useMutation({mutationFn:()=>api('/api/sync',{method:'POST'}),onSuccess:()=>{queryClient.invalidateQueries({queryKey:['samples']});queryClient.invalidateQueries({queryKey:['summary']});queryClient.invalidateQueries({queryKey:['observations']})}})
  const clinvarMutation = useMutation({mutationFn:()=>api('/api/clinvar/refresh',{method:'POST'}),onSuccess:()=>queryClient.invalidateQueries({queryKey:['summary']})})
  const reviewMutation = useMutation({mutationFn:(value:{id:string;status:string;note:string;flagged:boolean})=>api(`/api/reviews/${value.id}`,{method:'PUT',body:JSON.stringify(value)}),onSuccess:()=>{queryClient.invalidateQueries({queryKey:['observations']});queryClient.invalidateQueries({queryKey:['observations',filterKey,page]})}})
  const folderMutation = useMutation({mutationFn:(path:string)=>api('/api/settings/folder',{method:'PUT',body:JSON.stringify({path})}),onSuccess:()=>{queryClient.invalidateQueries({queryKey:['samples']});setSettingsOpen(false)}})
  const saveViewMutation = useMutation({mutationFn:(name:string)=>api('/api/views',{method:'POST',body:JSON.stringify({id:viewId||undefined,name,definition:{filter,columns:columns.map(c=>c.id).filter(Boolean)}})}),onSuccess:()=>queryClient.invalidateQueries({queryKey:['views']})})
  const rows = resultsQuery.data?.rows ?? []
  const table = useReactTable({data:rows,columns,getCoreRowModel:getCoreRowModel(),state:{columnVisibility},onColumnVisibilityChange:setColumnVisibility})
  const total = resultsQuery.data?.total ?? 0
  const pages = Math.max(1,Math.ceil(total/pageSize))

  const chooseSamples = (id:string) => {
    setSelectedSamples([id]);setPage(1)
  }
  const clearFilters = () => {setSearch('');setGene('');setEffect('');setReview('');setMinAf('');setMinDepth('');setRules([]);setLoadedFilter(null);setSelectedSamples([]);setPage(1)}
  const saveReview = () => {if(selected)reviewMutation.mutate({id:selected.id,status:selected.review_status,note:selected.note,flagged:selected.flagged})}
  const exportData = async () => {
    const response=await fetch('/api/export',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({filter})})
    if(!response.ok) return
    const blob=await response.blob();const url=URL.createObjectURL(blob);const a=document.createElement('a');a.href=url;a.download='mitoview-variantes.tsv';a.click();URL.revokeObjectURL(url)
  }
  const loadView=(v:{id:string;name:string;definition:Record<string,unknown>})=>{setViewId(v.id);setSearch('');setGene('');setEffect('');setReview('');setMinAf('');setMinDepth('');setRules([]);setSelectedSamples([]);const filterDef=v.definition.filter as FilterNode|undefined;setLoadedFilter(filterDef??null);setPage(1)}
  const issueCount=issues.length
  useEffect(()=>{setPage(1)},[filterKey])
  useEffect(()=>{if(selected){const fresh=rows.find(row=>row.id===selected.id);if(fresh&&fresh!==selected)setSelected(fresh)}},[rows,selected])

  return <div className="app-shell">
    <header className="topbar">
      <div className="brand"><div className="brand-mark"><Activity size={19}/></div><span>mito<span className="brand-light">view</span></span><span className="local-tag"><span className="local-dot"/>LOCAL</span></div>
      <div className="top-actions"><div className="sync-state"><span className="pulse"/>Carpeta sincronizada</div><button className="icon-button" aria-label="Ajustes de carpeta" onClick={()=>{setFolderInput(samplesQuery.data?.folder??'');setSettingsOpen(true)}}><Settings2 size={17}/></button><button className="button button-quiet" onClick={()=>refreshMutation.mutate()} disabled={refreshMutation.isPending}><RefreshCw size={15} className={refreshMutation.isPending?'spin':''}/> Actualizar datos</button></div>
    </header>
    <div className="workspace">
      <aside className="sidebar">
        <div className="side-label">ESPACIO DE TRABAJO</div>
        <button className={`nav-link ${tab==='variantes'?'active':''}`} onClick={()=>setTab('variantes')}><Database size={17}/>Explorar variantes</button>
        <button className={`nav-link ${tab==='comparar'?'active':''}`} onClick={()=>setTab('comparar')}><GitCompareArrows size={17}/>Comparar muestras</button>
        <button className={`nav-link ${tab==='revisión'?'active':''}`} onClick={()=>setTab('revisión')}><Bookmark size={17}/>Revisión guardada</button>
        <div className="side-divider"/>
        <div className="sample-head"><span className="side-label">MUESTRAS</span><span className="count-badge">{samples.length}</span></div>
        <div className="sample-list">{samples.map((sample,i)=><button key={sample.id} className={`sample-item ${effectiveSamples.includes(sample.id)?'selected':''}`} onClick={()=>chooseSamples(sample.id)}><span className={`sample-icon sample-color-${i%5}`}>{sample.label.slice(-2)}</span><span className="sample-name">{sample.label}</span><small>{sample.row_count??0}</small></button>)}</div>
        <div className="side-divider"/>
        <div className="side-label saved-title">FILTROS GUARDADOS</div>
        {viewsQuery.data?.map(view=><button className="saved-filter" key={view.id} onClick={()=>loadView(view)}><Bookmark size={14}/><span>{view.name}</span><MoreHorizontal size={14}/></button>)}
        {!viewsQuery.data?.length&&<div className="sidebar-empty">Guarda una combinación de filtros para encontrarla aquí.</div>}
        <div className="sidebar-foot"><div className="help-mark"><CircleHelp size={16}/></div><div><strong>Datos de origen</strong><span>NC_012920.1 · rCRS</span></div></div>
      </aside>

      <main className="main-content">
        <div className="page-heading"><div><div className="eyebrow">ANÁLISIS MITOCONDRIAL <span>·</span> ONT</div><h1>{tab==='variantes'?'Explorador de variantes':tab==='comparar'?'Comparación de muestras':'Revisión guardada'}</h1><p>{tab==='variantes'?'Busca, filtra y revisa las variantes de la muestra seleccionada.':tab==='comparar'?'Observa cómo cambia una misma variante entre muestras.':'Filtra las variantes marcadas o pendientes de revisión.'}</p></div><div className="heading-metric"><div className="metric-icon"><HeartPulse size={17}/></div><div><strong>{summaryQuery.data?.variants??'—'}</strong><span>variantes únicas</span></div></div></div>

        <div className="metric-row"><Metric label="Observaciones" value={summaryQuery.data?.observations??'—'} icon={<Database size={16}/>} detail="registros en TSV"/><Metric label="Muestras activas" value={samples.length} icon={<Activity size={16}/>} detail="en esta carpeta"/><Metric label="Anotación ClinVar" value={Object.values(summaryQuery.data?.clinvar??{}).reduce((a,b)=>a+b,0)} icon={<CircleHelp size={16}/>} detail="variantes consultadas"/></div>

        {issueCount>0&&<div className="alert-banner"><CircleHelp size={16}/><span><strong>{issueCount} incidencias al importar</strong> · Revisa las líneas indicadas antes de interpretar los resultados.</span><button onClick={()=>setSettingsOpen(true)}>Ver detalles</button></div>}

        <section className="data-card">
          <div className="card-head"><div className="tab-strip"><button className={tab==='variantes'?'tab active':''} onClick={()=>setTab('variantes')}>Todas las variantes</button><button className={tab==='revisión'?'tab active':''} onClick={()=>{setTab('revisión');setReview('pendiente')}}>Pendientes <span className="tab-count">{summaryQuery.data?.pending??'—'}</span></button></div><button className="button button-outline" onClick={exportData}><Download size={15}/> Exportar selección</button></div>
          <div className="filter-bar">
            <label className="search-box"><Search size={16}/><input placeholder="Buscar posición, gen o HGVS..." value={search} onChange={e=>{setSearch(e.target.value);setPage(1)}}/></label>
            <select aria-label="Filtrar por gen" value={gene} onChange={e=>{setGene(e.target.value);setPage(1)}}><option value="">Todos los genes</option>{[...new Set(rows.map(r=>r.gene).filter(Boolean))].sort().map(g=><option key={g!} value={g!}>{g}</option>)}</select>
            <select aria-label="Filtrar por efecto" value={effect} onChange={e=>{setEffect(e.target.value);setPage(1)}}><option value="">Todos los efectos</option>{['missense_variant','synonymous_variant','upstream_gene_variant','frameshift_variant'].map(e=><option key={e} value={e}>{e.replaceAll('_',' ')}</option>)}</select>
            <button className={`button ${advancedOpen?'button-active':'button-outline'} advanced-trigger`} onClick={()=>setAdvancedOpen(!advancedOpen)}><SlidersHorizontal size={15}/> Filtros <span className="filter-count">{rules.length+(minAf?1:0)+(minDepth?1:0)}</span><ChevronDown size={13}/></button>
          </div>
          {advancedOpen&&<div className="advanced-panel"><div className="advanced-title"><div><Filter size={16}/><strong>Filtros avanzados</strong><span>Combina criterios para acotar la búsqueda.</span></div><button className="icon-button" onClick={()=>setAdvancedOpen(false)} aria-label="Cerrar filtros"><X size={16}/></button></div>
            <div className="quick-thresholds"><label>Fracción alélica mín.<div className="input-suffix"><input type="number" min="0" max="100" step="0.1" value={minAf} onChange={e=>setMinAf(e.target.value)} placeholder="0"/><span>%</span></div></label><label>Profundidad mín.<div className="input-suffix"><input type="number" min="0" value={minDepth} onChange={e=>setMinDepth(e.target.value)} placeholder="0"/><span>×</span></div></label><label>Estado revisión<select value={review} onChange={e=>setReview(e.target.value)}><option value="">Cualquiera</option><option value="pendiente">Pendiente</option><option value="en_revision">En revisión</option><option value="revisada">Revisada</option></select></label></div>
            <div className="rule-top"><strong>Condiciones personalizadas</strong><select value={advancedOp} onChange={e=>setAdvancedOp(e.target.value as 'and'|'or'|'not')}><option value="and">Cumplir todas</option><option value="or">Cumplir alguna</option><option value="not">Excluir este grupo</option></select></div>
            {rules.map((rule,index)=><div className="rule-row" key={index}><select aria-label="Campo" value={rule.field} onChange={e=>updateRule(setRules,index,{field:e.target.value,operator:numFields.has(e.target.value)?'gte':'eq'})}><option value="">Campo...</option>{fieldOptions.map(([v,l])=><option key={v} value={v}>{l}</option>)}</select><select aria-label="Operador" value={rule.operator} onChange={e=>updateRule(setRules,index,{operator:e.target.value})}><option value="eq">es</option><option value="neq">no es</option>{numFields.has(rule.field)&&<><option value="gt">mayor que</option><option value="gte">mayor o igual</option><option value="lt">menor que</option><option value="lte">menor o igual</option><option value="between">entre</option></>}<option value="contains">contiene</option><option value="isMissing">sin dato</option><option value="isPresent">con dato</option></select>{!['isMissing','isPresent'].includes(rule.operator)&&<div className="rule-values"><input aria-label="Valor" value={rule.value} onChange={e=>updateRule(setRules,index,{value:e.target.value})} placeholder="Valor"/>{rule.operator==='between'&&<input aria-label="Valor máximo" value={rule.value2} onChange={e=>updateRule(setRules,index,{value2:e.target.value})} placeholder="Hasta"/>}</div>}<button className="icon-button danger-hover" aria-label="Eliminar condición" onClick={()=>setRules(old=>old.filter((_,i)=>i!==index))}><X size={15}/></button></div>)}
            <div className="advanced-footer"><button className="text-button" onClick={()=>setRules(old=>[...old,{field:'',operator:'eq',value:'',value2:''}])}><Plus size={14}/> Añadir condición</button><button className="text-button" onClick={clearFilters}>Limpiar filtros</button></div>
          </div>}
          <div className="saved-view-row"><Bookmark size={14}/><select value={viewId} onChange={e=>{const v=viewsQuery.data?.find(x=>x.id===e.target.value);if(v)loadView(v)}}><option value="">Aplicar filtro guardado...</option>{viewsQuery.data?.map(v=><option value={v.id} key={v.id}>{v.name}</option>)}</select><button className="text-button" onClick={()=>{const name=prompt('Nombre para esta vista');if(name)saveViewMutation.mutate(name)}}><Plus size={14}/> Guardar vista actual</button></div>

          <div className="result-toolbar"><div><strong>{total.toLocaleString('es-ES')}</strong> resultados <span>·</span> {samples.find(s=>effectiveSamples.includes(s.id))?.label??'Sin muestra seleccionada'}</div><div className="result-actions"><details className="column-picker"><summary><Settings2 size={14}/> Columnas</summary><div>{table.getAllLeafColumns().map(column=><label key={column.id}><input type="checkbox" checked={column.getIsVisible()} onChange={column.getToggleVisibilityHandler()}/>{String(column.columnDef.header??column.id)}</label>)}</div></details><button className={`text-button ${showCharts?'selected-text':''}`} onClick={()=>setShowCharts(!showCharts)}><Activity size={15}/>{showCharts?'Ocultar':'Gráficos'}</button><button className="text-button" onClick={()=>clinvarMutation.mutate()} disabled={clinvarMutation.isPending}><ExternalLink size={14}/>Actualizar ClinVar</button><button className="text-button" onClick={()=>setTab('comparar')}><GitCompareArrows size={15}/>Comparar</button></div></div>
          {showCharts&&<div className="charts-section"><div className="chart-tabs"><button className={chartTab==='alineamiento'?'active':''} onClick={()=>setChartTab('alineamiento')}>Datos de alineamiento</button><button className={chartTab==='variantes'?'active':''} onClick={()=>setChartTab('variantes')}>Datos de variantes</button></div>{chartTab==='alineamiento'?<div className="chart-card alignment-card"><div className="chart-head"><strong>Profundidad de lecturas alineadas · NC_012920</strong><FolderPicker onFiles={setCoverageFiles} extensions={['.tsv']}/></div><p className="chart-note">Selecciona la carpeta Coverage con los TSV de minimap2 y ngmlr. {coverageFiles.length?`${coverageFiles.length} archivos cargados.`:'Los datos se procesan localmente.'}</p><div className="chart-body alignment-body">{coverageFiles.length?<CoverageCharts files={coverageFiles} selectedSample={samples.find(s=>effectiveSamples.includes(s.id))?.label??''}/>:<div className="chart-empty">Selecciona la carpeta Coverage para visualizar las profundidades.</div>}</div></div>:<div className="charts-panel"><div className="chart-card"><div className="chart-head"><strong>Efectos anotados</strong><span>{chartsQuery.data?.count??0} observaciones</span></div><div className="chart-body">{chartsQuery.data?.effects.length?<ResponsiveContainer width="100%" height="100%"><BarChart data={chartsQuery.data.effects} layout="vertical" margin={{left:12,right:12,top:4,bottom:2}}><XAxis type="number" hide/><YAxis type="category" dataKey="name" width={145} tick={{fontSize:11,fill:'#536078'}}/><Tooltip/><Bar dataKey="count" radius={[0,5,5,0]}>{chartsQuery.data.effects.map((_,i)=><Cell key={i} fill={['#3867e8','#53a8a1','#b58bdb','#e0a34b'][i%4]}/>)}</Bar></BarChart></ResponsiveContainer>:<div className="chart-empty">Sin datos.</div>}</div></div><div className="chart-card"><div className="chart-head"><strong>Fracción alélica y profundidad</strong><span>cada punto es una observación</span></div><div className="chart-body">{chartsQuery.data?.points.length?<ResponsiveContainer width="100%" height="100%"><ScatterChart margin={{top:8,right:12,bottom:4,left:2}}><CartesianGrid strokeDasharray="3 3" stroke="#e8edf3"/><XAxis type="number" dataKey="af" name="Fracción alélica" unit="%" domain={[0,100]} tick={{fontSize:11,fill:'#748198'}}/><YAxis type="number" dataKey="coverage" name="Profundidad" unit="×" tick={{fontSize:11,fill:'#748198'}}/><Tooltip cursor={{strokeDasharray:'3 3'}} formatter={(value:any,name:any)=>[value,name==='af'?'Fracción alélica (%)':'Profundidad']} labelFormatter={()=>''}/><Scatter data={chartsQuery.data.points} fill="#3867e8" onClick={(data:any)=>{const row=rows.find(r=>r.id===data.id);if(row)setSelected(row)}}/></ScatterChart></ResponsiveContainer>:<div className="chart-empty">Sin observaciones para graficar.</div>}</div></div><div className="chart-card haplo-card"><div className="chart-head"><strong>Haplogrupos</strong><FolderPicker onFiles={setHaplogroupFiles} extensions={['.txt']}/></div><HaplogroupInfo files={haplogroupFiles} selectedSample={samples.find(s=>effectiveSamples.includes(s.id))?.label??''}/></div></div>}</div>}

          {tab==='comparar'?<ComparePanel samples={samples} compareIds={compareIds??samples.slice(0,4).map(s=>s.id)} setCompareIds={setCompareIds} filter={filter} onSelect={id=>{setTab('variantes');setSelectedSamples([]);setSearch('');setTimeout(()=>document.getElementById(`row-${id}`)?.scrollIntoView({block:'center'}),50)}}/>:tab==='revisión'?<div className="review-shortcuts"><div><Bookmark size={18}/><strong>Tu revisión</strong><span>Usa los filtros de arriba para mostrar pendientes, notas o variantes marcadas.</span></div><div className="review-toggle"><button className={review==='pendiente'?'on':''} onClick={()=>setReview(review==='pendiente'?'':'pendiente')}>Pendientes</button><button onClick={()=>setReview(review==='revisada'?'':'revisada')}>Revisadas</button><button onClick={clearFilters}>Ver todas</button></div></div>:<div className="table-wrap"><table className="variant-table"><thead><tr>{table.getHeaderGroups().map(hg=>hg.headers.map(h=><th key={h.id}><button className="sort-header" onClick={()=>{const id=h.column.id;const apiField=({af_alt:"afAlt",review_status:"reviewStatus"} as Record<string,string>)[id]||id;if(sortBy===apiField)setSortDirection(sortDirection==='asc'?'desc':'asc');else{setSortBy(apiField);setSortDirection('asc')}}}>{flexRender(h.column.columnDef.header,h.getContext())}<ArrowDownUp size={12}/></button></th>))}</tr></thead><tbody>{resultsQuery.isLoading?<tr><td className="empty-cell" colSpan={8}>Cargando variantes…</td></tr>:rows.length?table.getRowModel().rows.map(row=><tr id={`row-${row.original.id}`} key={row.id} className={selected?.id===row.original.id?'row-selected':''} onClick={()=>setSelected(row.original)}>{row.getVisibleCells().map(cell=><td key={cell.id}>{flexRender(cell.column.columnDef.cell,cell.getContext())}</td>)}</tr>):<tr><td className="empty-cell" colSpan={8}>{samples.length?'No hay resultados con estos filtros.':'Esperando archivos TSV en la carpeta configurada…'}</td></tr>}</tbody></table></div>}
          {tab!=='comparar'&&<div className="table-footer"><span>Mostrando {total?((page-1)*pageSize+1):0}–{Math.min(page*pageSize,total)} de {total.toLocaleString('es-ES')}</span><div className="pagination"><button aria-label="Página anterior" disabled={page<=1} onClick={()=>setPage(page-1)}><ChevronLeft size={16}/></button><span>Página <strong>{page}</strong> de {pages}</span><button aria-label="Página siguiente" disabled={page>=pages} onClick={()=>setPage(page+1)}><ChevronRight size={16}/></button></div></div>}
        </section>
        <footer className="page-footer"><span>Los datos provienen de los TSV originales y de las fuentes citadas en cada anotación.</span><span>Revisión manual · rCRS NC_012920.1</span></footer>
      </main>
    </div>
    {selected&&<aside className="detail-backdrop" onClick={()=>setSelected(null)}><section className="detail-drawer" onClick={e=>e.stopPropagation()}><div className="drawer-head"><div><div className="eyebrow">DETALLE DE VARIANTE</div><h2 className="mono">{selected.chrom}:{selected.pos}</h2><span className="drawer-sub">{selected.ref} → {selected.alt} <span>·</span> {selected.sample_label}</span></div><button className="icon-button" onClick={()=>setSelected(null)} aria-label="Cerrar detalle"><PanelRightClose size={18}/></button></div>
      <div className="drawer-scroll"><div className="detail-tags"><span className={`impact ${selected.impact?.toLowerCase()}`}>{selected.impact||'Sin impacto'}</span><span className="detail-effect">{selected.effect?.replaceAll('_',' ')||'Sin efecto anotado'}</span></div>
      <div className="detail-grid"><Detail label="Gen" value={selected.gene}/><Detail label="HGVS c." value={selected.hgvs_c}/><Detail label="HGVS p." value={selected.hgvs_p}/><Detail label="Fracción alélica" value={selected.af_alt==null?'—':`${(selected.af_alt*100).toFixed(2)}%`}/><Detail label="Profundidad" value={selected.coverage==null?'—':`${selected.coverage} lecturas`}/><Detail label="AD (ref, alt)" value={`${selected.ad_ref??'—'}, ${selected.ad_alt??'—'}`}/><Detail label="gnomAD homoplasmia" value={selected.af_hom==null?'Sin dato':selected.af_hom.toExponential(2)}/><Detail label="gnomAD heteroplasmia" value={selected.af_het==null?'Sin dato':selected.af_het.toExponential(2)}/></div>
      <div className="section-divider"/><div className="section-heading"><div><strong>ClinVar</strong><span>Fuente externa · actualización automática</span></div><button className="text-button" onClick={()=>clinvarMutation.mutate()}><RefreshCw size={13}/>Consultar</button></div>
      {selected.clinvar_records?.length?selected.clinvar_records.map((record,i)=><article className="clinvar-record" key={`${record.uid}-${i}`}><div className="clinvar-top"><strong>{record.accession||`Registro ${record.uid}`}</strong><a href={record.url} target="_blank" rel="noreferrer" aria-label="Abrir registro ClinVar"><ExternalLink size={14}/></a></div><div className="clinvar-class">{record.classification||'Clasificación no indicada'}</div>{record.review_status&&<div className="clinvar-review">{record.review_status}</div>}{record.conditions?.length>0&&<div className="condition-list">{record.conditions.map(condition=><span key={condition}>{condition}</span>)}</div>}<small>Actualizado {record.updated||selected.clinvar_checked_at||'—'}</small></article>):<div className="clinvar-empty"><div className="clinvar-state">{clinvarLabel(selected.clinvar_status)}</div>{selected.clinvar_error&&<small>{selected.clinvar_error}</small>}<span>Las anotaciones de ClinVar son una referencia para la revisión profesional.</span></div>}
      <div className="section-divider"/><div className="section-heading"><div><strong>Revisión médica</strong><span>Notas locales asociadas a esta versión</span></div><Flag size={16} className={selected.flagged?'flag-active':''}/></div>
      <label className="detail-field">Estado<select value={selected.review_status} onChange={e=>setSelected({...selected,review_status:e.target.value as Observation['review_status']})}><option value="pendiente">Pendiente</option><option value="en_revision">En revisión</option><option value="revisada">Revisada</option></select></label>
      <label className="detail-field">Nota<textarea rows={4} value={selected.note} onChange={e=>setSelected({...selected,note:e.target.value})} placeholder="Añade observaciones para esta variante…"/></label>
      <label className="flag-checkbox"><input type="checkbox" checked={selected.flagged} onChange={e=>setSelected({...selected,flagged:e.target.checked})}/> <Flag size={14}/>Marcar para revisar después</label>
      <div className="source-box"><Database size={14}/><div><strong>Archivo de origen</strong><span>{selected.sample_path.split(/[\\/]/).pop()}</span><span>Importado {new Date(selected.imported_at).toLocaleString('es-ES')}</span></div></div>
      <details className="raw-details"><summary>Ver todos los campos del TSV</summary><div className="raw-fields">{Object.entries(selected.raw).map(([k,v])=><div key={k}><span>{k}</span><code>{v||'—'}</code></div>)}</div></details>
      </div><div className="drawer-footer"><button className="button button-outline" onClick={()=>setSelected(null)}>Cerrar</button><button className="button button-primary" onClick={saveReview} disabled={reviewMutation.isPending}><Check size={15}/>Guardar revisión</button></div>
    </section></aside>}
    {settingsOpen&&<div className="modal-backdrop" onClick={()=>setSettingsOpen(false)}><section className="settings-modal" onClick={e=>e.stopPropagation()}><div className="modal-title"><div className="metric-icon"><FolderSync size={18}/></div><div><h2>Carpeta de resultados</h2><p>Se comprueba automáticamente cada cinco segundos.</p></div><button className="icon-button" onClick={()=>setSettingsOpen(false)} aria-label="Cerrar"><X size={17}/></button></div><label className="detail-field">Ruta local<input className="folder-path" value={folderInput} onChange={e=>setFolderInput(e.target.value)} placeholder="C:\\ruta\\TSV results"/></label><button className="button button-primary save-folder" onClick={()=>folderMutation.mutate(folderInput)} disabled={folderMutation.isPending}><Check size={15}/>Guardar y sincronizar</button>{issues.length>0&&<div className="issue-list"><strong>Incidencias de importación</strong>{issues.map((issue,i)=><div key={i}><code>{issue.path.split(/[\\/]/).pop()}{issue.line?` · línea ${issue.line}`:''}</code><span>{issue.message}</span></div>)}</div>}<div className="privacy-note"><span className="local-dot"/>La web solo se abre en este equipo. ClinVar recibe únicamente la descripción de la variante.</div></section></div>}
  </div>
}

function numValue(value:string): unknown { if(value==='')return ''; const n=Number(value); return Number.isFinite(n)?n:value }

function FolderPicker({onFiles,extensions}:{onFiles:(files:File[])=>void;extensions:string[]}) {
  return <label className="folder-picker"><FolderOpen size={14}/> Seleccionar carpeta<input type="file" multiple {...({webkitdirectory:'',directory:''} as any)} accept={extensions.join(',')} onChange={e=>onFiles(Array.from(e.target.files??[]).filter(f=>extensions.some(ext=>f.name.toLowerCase().endsWith(ext))))}/></label>
}

function CoverageCharts({files,selectedSample}:{files:File[];selectedSample:string}) {
  const [series,setSeries]=useState<{name:string;data:{position:number;depth:number}[]}[]>([])
  useEffect(()=>{let alive=true;Promise.all(files.map(async file=>({file,text:await file.text()}))).then(items=>{
    const relevant=items.filter(({file})=>!selectedSample||file.webkitRelativePath.toLowerCase().includes(selectedSample.toLowerCase())||file.name.toLowerCase().includes(selectedSample.toLowerCase()))
    const source=relevant.length?relevant:items
    const parsed=source.map(({file,text})=>{const lines=text.split(/\r?\n/).map(line=>line.trim()).filter(Boolean);let positionCol=0,depthCol=1;const first=lines[0]?.split(/[\t,; ]+/).map(x=>x.toLowerCase())??[];const posHeader=first.findIndex(x=>/^(pos|position|start|start_position|coordinate)$/.test(x));const depthHeader=first.findIndex(x=>/(depth|coverage|read_depth|meandepth)/.test(x));if(posHeader>=0)positionCol=posHeader;if(depthHeader>=0)depthCol=depthHeader;const data=lines.slice(posHeader>=0||depthHeader>=0?1:0).map(line=>{const cols=line.split(/[\t,; ]+/);const nums=cols.map(Number);let p=nums[positionCol],d=nums[depthCol];if(!Number.isFinite(p)||!Number.isFinite(d)){const numeric=nums.map((v,i)=>({v,i})).filter(x=>Number.isFinite(x.v));if(numeric.length>=2){p=numeric[0].v;d=numeric[numeric.length-1].v}}if(!Number.isFinite(p)||!Number.isFinite(d))return null;return {position:p,depth:d}}).filter((x):x is {position:number;depth:number}=>x!==null);const path=file.webkitRelativePath.toLowerCase();const aligner=path.includes('ngmlr')?'ngmlr':path.includes('minimap')?'minimap2':file.name.replace(/\.tsv$/i,'');return {name:aligner,data}}).filter(x=>x.data.length)
    if(alive)setSeries(parsed)
  });return()=>{alive=false}},[files,selectedSample])
  return series.length?<ResponsiveContainer width="100%" height="100%"><ScatterChart margin={{top:12,right:18,bottom:12,left:12}}><CartesianGrid strokeDasharray="3 3" stroke="#e8edf3"/><XAxis type="number" dataKey="position" name="Posición mitocondrial" domain={[1,16569]} tick={{fontSize:10,fill:'#748198'}} label={{value:'Posición en NC_012920 (pb)',position:'insideBottom',offset:-7,fontSize:10}}/><YAxis type="number" dataKey="depth" name="Profundidad" tick={{fontSize:10,fill:'#748198'}} label={{value:'Profundidad (×)',angle:-90,position:'insideLeft',fontSize:10}}/><Tooltip/><Legend/><Scatter name={series.length===1?series[0].name:'minimap2'} data={series[0].data} fill="#3867e8" line/><Scatter name={series.length>1?series[1].name:'ngmlr'} data={series[1]?.data??[]} fill="#53a8a1" line/></ScatterChart></ResponsiveContainer>:<div className="chart-empty">No se reconocieron pares posición/profundidad en los TSV. Se esperan columnas numéricas posición y profundidad.</div>
}

function HaplogroupInfo({files,selectedSample}:{files:File[];selectedSample:string}) {
  const [text,setText]=useState('')
  useEffect(()=>{let alive=true;const match=files.find(f=>!selectedSample||f.webkitRelativePath.toLowerCase().includes(selectedSample.toLowerCase())||f.name.toLowerCase().includes(selectedSample.toLowerCase()))??files[0];if(!match){setText('');return}match.text().then(value=>{if(alive)setText(value)});return()=>{alive=false}},[files,selectedSample])
  return <div className="haplo-content">{files.length?<><div className="haplo-file">{files.find(f=>f.webkitRelativePath.toLowerCase().includes(selectedSample.toLowerCase())||f.name.toLowerCase().includes(selectedSample.toLowerCase()))?.name??files[0].name}</div><pre>{text||'El archivo está vacío.'}</pre></>:<div className="chart-empty">Selecciona la carpeta Haplogroups con los TXT de las muestras.</div>}</div>
}
function updateRule(setter:React.Dispatch<React.SetStateAction<{field:string;operator:string;value:string;value2:string}[]>>,index:number,patch:Partial<{field:string;operator:string;value:string;value2:string}>){setter(old=>old.map((rule,i)=>i===index?{...rule,...patch}:rule))}
function Metric({label,value,icon,detail}:{label:string;value:string|number;icon:React.ReactNode;detail:string}){return <div className="metric-card"><div className="metric-icon">{icon}</div><div className="metric-copy"><span>{label}</span><div><strong>{typeof value==='number'?value.toLocaleString('es-ES'):value}</strong><small>{detail}</small></div></div></div>}
function Detail({label,value}:{label:string;value:React.ReactNode}){return <div className="detail-item"><span>{label}</span><strong>{value||'—'}</strong></div>}
function clinvarLabel(status?:string){return ({pendiente:'Consulta pendiente',coincidencia_verificada:'Coincidencia verificada',sin_coincidencia:'Sin registro coincidente',ambigua:'Más de una coincidencia',error:'No se pudo consultar ClinVar'} as Record<string,string>)[status||'']||'Consultando ClinVar…'}

function ComparePanel({samples,compareIds,setCompareIds,filter,onSelect}:{samples:Sample[];compareIds:string[];setCompareIds:(v:string[])=>void;filter:FilterNode;onSelect:(id:string)=>void}) {
  const ids=compareIds.filter(id=>samples.some(s=>s.id===id))
  const query=useQuery({queryKey:['compare',ids.join(','),JSON.stringify(filter)],queryFn:()=>api<{samples:{id:string;label:string}[];variants:{key:string;cells:Record<string,{af_alt:number|null;coverage:number|null;observation_id:string}|null>}[]}>('/api/compare',{method:'POST',body:JSON.stringify({sampleIds:ids,filter})}),enabled:ids.length>0})
  return <div className="compare-panel"><div className="compare-controls"><div><strong>Selecciona hasta 12 muestras</strong><span>Las variantes se agrupan por posición y alelos.</span></div><div className="compare-checks">{samples.map(s=><label key={s.id}><input type="checkbox" checked={ids.includes(s.id)} onChange={e=>{const next=e.target.checked?[...ids,s.id].slice(-12):ids.filter(x=>x!==s.id);setCompareIds(next)}}/>{s.label}</label>)}</div></div>{query.isLoading?<div className="compare-empty">Preparando comparación…</div>:ids.length===0?<div className="compare-empty">Selecciona al menos una muestra para comparar.</div>:<div className="compare-table-wrap"><table className="compare-table"><thead><tr><th>Variante</th>{query.data?.samples.map(s=><th key={s.id}>{s.label}</th>)}</tr></thead><tbody>{query.data?.variants.length?query.data.variants.map(v=><tr key={v.key}><td className="mono">{v.key}</td>{query.data?.samples.map(s=>{const cell=v.cells[s.id];return <td key={s.id} className={cell?'compare-found':'compare-missing'}>{cell?<button onClick={()=>onSelect(cell.observation_id)}>{cell.af_alt==null?'—':`${(cell.af_alt*100).toFixed(1)}%`} <small>· {cell.coverage??'—'}×</small></button>:<span>No consta en el TSV</span>}</td>})}</tr>):<tr><td colSpan={ids.length+1} className="empty-cell">No hay variantes en las muestras seleccionadas con estos filtros.</td></tr>}</tbody></table></div>}<div className="compare-legend"><span className="legend-dot found"/> Variedad registrada <span className="legend-dot missing"/> Sin registro en el TSV</div></div>
}

export default App
