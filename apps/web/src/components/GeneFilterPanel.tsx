import { useState } from "react";
import type { ProjectPayload } from "../types";

export type GeneFilterSettings = {
  exclude_mitochondrial: boolean; exclude_pseudogenes: boolean; exclude_rrna: boolean;
  protein_coding_only: boolean; excluded_genes: string[]; min_expression: number; min_samples: number;
  annotation_path: string; expression_matrix_path: string;
};
const defaults: GeneFilterSettings = {exclude_mitochondrial:false, exclude_pseudogenes:false,
  exclude_rrna:false, protein_coding_only:false, excluded_genes:[], min_expression:0, min_samples:1,
  annotation_path:"", expression_matrix_path:""};

export function GeneFilterPanel({project, busy, error, onApply}: {
  project: ProjectPayload; busy: boolean; error?: string;
  onApply: (settings: GeneFilterSettings) => void;
}) {
  const saved = project.manifest.metadata.gene_filter as {
    config?: Partial<GeneFilterSettings>; annotation_tsv?: string; gtf?: string; expression_matrix?: string;
  } | undefined;
  const [settings, setSettings] = useState<GeneFilterSettings>(() => ({...defaults, ...saved?.config,
    annotation_path:saved?.annotation_tsv ?? saved?.gtf ?? "",
    expression_matrix_path:saved?.expression_matrix ?? ""}));
  const [list, setList] = useState(settings.excluded_genes.join("\n"));
  const [fileError, setFileError] = useState("");
  const enabled = settings.exclude_mitochondrial || settings.exclude_pseudogenes || settings.exclude_rrna
    || settings.protein_coding_only || list.trim() || settings.min_expression > 0;
  return <details className="run-plan-panel gene-filter-panel">
    <summary><strong>Gene exclusion / 基因排除</strong> · {saved ? "Enabled / 已启用" : "All genes / 未排除"}</summary>
    <fieldset disabled={busy}>
      <legend>Bulk result scope / bulk 结果范围</legend>
      {([
        ["exclude_mitochondrial", "Exclude mitochondrial-encoded genes / 排除线粒体编码基因"],
        ["exclude_pseudogenes", "Exclude pseudogenes / 排除假基因"],
        ["exclude_rrna", "Exclude rRNA / 排除核糖体 RNA"],
        ["protein_coding_only", "Protein-coding only (also excludes lncRNA and unannotated genes) / 仅保留蛋白编码基因"],
      ] as const).map(([key,label]) => <label key={key} className="gene-filter-check"><input type="checkbox"
        checked={settings[key]} onChange={event=>setSettings({...settings,[key]:event.target.checked})}/>{label}</label>)}
      <label>Annotation TSV / GTF path (optional override; uses the project GTF otherwise)
        <input value={settings.annotation_path} onChange={e=>setSettings({...settings,annotation_path:e.target.value})} placeholder="Local annotation file with gene type and chromosome"/></label>
      <label>Additional gene IDs / symbols to exclude / 自定义排除名单
        <textarea value={list} onChange={e=>setList(e.target.value)} rows={3}/></label>
      <label>Load a text gene list / 导入名单<input type="file" accept=".txt,.tsv,.csv" onChange={async e=>{
        const file=e.target.files?.[0]; if (!file) return;
        if (file.size>2*1024*1024) {setFileError("Gene list must be smaller than 2 MiB."); return;}
        try {const text=await file.text();setList(previous=>previous+"\n"+text);setFileError("");}
        catch {setFileError("Could not read the gene list.");}
      }}/></label>
      <label>Minimum expression (0 disables) / 最低表达量<input type="number" min="0" step="any"
        value={settings.min_expression} onChange={e=>setSettings({...settings,min_expression:Number(e.target.value)})}/></label>
      <label>Minimum samples / 最少样本数<input type="number" min="1" step="1" value={settings.min_samples}
        onChange={e=>setSettings({...settings,min_samples:Number(e.target.value)})}/></label>
      {settings.min_expression>0 && <label>Expression matrix path (gene_id + sample columns; units follow this file)
        <input value={settings.expression_matrix_path} onChange={e=>setSettings({...settings,expression_matrix_path:e.target.value})}/></label>}
      <p>Creates a separate filtered project and preserves existing results. Rebuild downstream summaries, PCA/heatmaps and enrichment in that scope. Original upstream P values and FDR stay unchanged. No FASTQ/BAM reprocessing.</p>
      <p>按注释识别基因类别；核编码的线粒体相关基因不会因名称被删除。报告自动记录排除规则和数量。</p>
      {(error || fileError) && <p role="alert">{error || fileError}</p>}
      <button type="button" className="secondary-action" disabled={busy || !Number.isFinite(settings.min_expression)
        || settings.min_expression<0 || !Number.isInteger(settings.min_samples) || settings.min_samples<1}
        onClick={()=>onApply({...settings,excluded_genes:[...new Set(list.split(/[\s,;]+/).filter(Boolean))]})}>
        {busy ? "Applying…" : enabled ? "Apply exclusions / 应用排除" : "Use original unfiltered project / 使用未过滤项目"}</button>
    </fieldset>
  </details>;
}
