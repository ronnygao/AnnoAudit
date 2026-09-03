# -*- coding: utf-8 -*-
# v9c M2 敏感性分析：
# ① 官方谷氨酸能的 marker 表达谱（正面证据：神经元 marker 低、小胶质 marker 高）
# ② score_genes（AddModuleScore 风格）第二打分器复核污染率
import os, gzip, glob, time
import numpy as np, pandas as pd, scipy.io, scipy.sparse as sp

t0 = time.time()
RAW = "/Users/apple/WorkBuddy/2026-07-22-16-57-18/cerebri_data/raw"
ANNOT = "/Users/apple/Desktop/颅脑损伤专病数据库/CEREBRI_KCNC3生信分析/v3分析结果/cell_annotations.csv"
OUT = "/Users/apple/Desktop/颅脑损伤专病数据库/CEREBRI_KCNC3生信分析/v9_unbiased_clustering"
os.makedirs(OUT, exist_ok=True)
def log(m):
    print(f"[{time.time()-t0:6.1f}s] {m}", flush=True)

MARKERS = {
    "Excitatory": ["Slc17a7", "Rbfox3", "Neurod6", "Tbr1", "Camk2a"],
    "Inhibitory": ["Gad1", "Gad2", "Sst", "Pvalb", "Npy"],
    "Astrocyte": ["Gfap", "Aqp4", "Slc1a3", "Aldh1l1"],
    "Microglia": ["C1qb", "Tyrobp", "Cx3cr1", "Hexb", "Aif1"],
    "Oligodendrocyte": ["Mbp", "Plp1", "Mog", "Mag"],
    "OPC": ["Pdgfra", "Olig1", "Cspg4"],
    "Endothelial": ["Pecam1", "Cldn5", "Flt1"],
    "Pericyte": ["Pdgfrb", "Rgs5", "Vtn"],
}
ann = pd.read_csv(ANNOT)
ann["bc_plain"] = ann["barcode"].str.replace(r"_GSM\d+$", "", regex=True)
marker_genes = sorted({g for v in MARKERS.values() for g in v})
# 扩充基因集（marker + 离子通道）供 score_genes 构建对照基因
ion_genes = [g.strip() for g in open("/tmp/ionchannels_mouse.txt") if g.strip()]
need_all = sorted(set(marker_genes) | set(ion_genes))
feat_files = sorted(os.listdir(RAW))
feat_path = next(f for f in feat_files if f.endswith("_features.tsv.gz"))
with gzip.open(os.path.join(RAW, feat_path), "rt") as f:
    features = [l.strip().split("\t") for l in f if not l.startswith("#")]
gene_names = [x[1] if len(x) > 1 else x[0] for x in features]
gene_to_idx = {g: i for i, g in enumerate(gene_names)}
need_idx = {g: gene_to_idx[g] for g in need_all if g in gene_to_idx}
need_cols = [need_idx[g] for g in need_all if g in need_idx]
loaded_genes = [g for g in need_all if g in need_idx]
log(f"基因集: marker {len(marker_genes)} + ion {len(ion_genes)} -> 命中 {len(need_cols)}")

X_parts, meta_parts = [], []
for gsm, grp in ann.groupby("gsm"):
    mtx_f = glob.glob(os.path.join(RAW, f"{gsm}_*_matrix.mtx.gz"))
    bar_f = glob.glob(os.path.join(RAW, f"{gsm}_*_barcodes.tsv.gz"))
    if not mtx_f or not bar_f: continue
    with gzip.open(mtx_f[0], "rb") as f:
        mat = scipy.io.mmread(f)
    mat = sp.csc_matrix(mat).T.tocsr()
    with gzip.open(bar_f[0], "rt") as f:
        barcodes = [l.strip() for l in f]
    bc2row = {bc: i for i, bc in enumerate(barcodes)}
    keep = grp["bc_plain"].isin(barcodes)
    grp2 = grp.loc[grp.index[keep]]
    rows = [bc2row[b] for b in grp2["bc_plain"]]
    if not rows: continue
    X_parts.append(mat[rows, :][:, need_cols].tocsr()); meta_parts.append(grp2)
X = sp.vstack(X_parts).tocsr()
meta = pd.concat(meta_parts).reset_index(drop=True)
meta["total_umi"] = meta["total_umi"].astype(float)
qc = (meta["n_genes"] > 200) & (meta["mito_pct"] < 20)
meta = meta[qc].reset_index(drop=True); X = X[qc.values]
umi = meta["total_umi"].values
X_cpm = X.multiply(1e6 / umi[:, None]).tocsr()
X_log = X_cpm.copy(); X_log.data = np.log1p(X_log.data)
g2c = {g: i for i, g in enumerate(loaded_genes)}
official = (meta["cell_type"] == "Glutamatergic_neuron").values
log(f"矩阵 {X.shape} | 官方谷氨酸能 {official.sum()}")

# ============ ① 官方谷氨酸能 marker 表达谱（正面证据） ============
log("\n=== ① 官方谷氨酸能 marker 表达谱 ===")
Xo = X_log[official].toarray()
rows = []
for ct, genes in MARKERS.items():
    idxs = [g2c[g] for g in genes if g in g2c]
    m_ct = Xo[:, idxs].mean()
    det = (Xo[:, idxs] > 0).mean() * 100
    rows.append({"marker_set": ct, "mean_expr_in_official_glut": round(m_ct, 3), "detection_pct": round(det, 1)})
marker_profile = pd.DataFrame(rows)
log(marker_profile.to_string(index=False))
marker_profile.to_csv(f"{OUT}/official_glut_marker_profile.csv", index=False)

# 对比：真谷氨酸能（max-score Excitatory）的 marker 谱
score_df = pd.DataFrame(index=meta.index)
for ct, genes in MARKERS.items():
    idxs = [g2c[g] for g in genes if g in g2c]
    score_df[ct] = np.asarray(X_log[:, idxs].mean(axis=1)).ravel()
max_ct = score_df.idxmax(axis=1); max_val = score_df.max(axis=1)
true_glut = ((max_ct == "Excitatory") & (max_val > 0)).values
Xt = X_log[true_glut].toarray()
rows2 = []
for ct, genes in MARKERS.items():
    idxs = [g2c[g] for g in genes if g in g2c]
    rows2.append({"marker_set": ct, "mean_expr_in_true_glut": round(Xt[:, idxs].mean(), 3)})
log(f"\n真谷氨酸能 n={true_glut.sum()}，其 marker 谱:")
log(pd.DataFrame(rows2).to_string(index=False))

# ============ ② score_genes（AddModuleScore 风格）第二打分器 ============
log("\n=== ② score_genes 第二打分器（AddModuleScore 风格） ===")
import scanpy as sc
adata = sc.AnnData(X.astype(np.float32))
adata.obs_names = meta["barcode"].values
adata.var_names = loaded_genes
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
for ct, genes in MARKERS.items():
    gs = [g for g in genes if g in adata.var_names]
    sc.tl.score_genes(adata, gene_list=gs, score_name=f"score_{ct}")
scores = pd.DataFrame({ct: adata.obs[f"score_{ct}"].values for ct in MARKERS})
sg_max = scores.idxmax(axis=1)
sg_val = scores.max(axis=1)
sg_ct = np.where(sg_val > 0, sg_max, "Unknown")
log("score_genes 全库注释分布:\n" + pd.Series(sg_ct).value_counts().to_string())

# 官方谷氨酸能的 score_genes 构成
off_comp = pd.Series(sg_ct[official]).value_counts()
off_comp_pct = (off_comp / official.sum() * 100).round(1)
log(f"\n官方谷氨酸能 ({official.sum()}) 的 score_genes 注释构成:")
for ct, n in off_comp.items():
    log(f"  {ct:>14}: {n:>7} ({off_comp_pct[ct]:.1f}%)")
n_neuron_sg = off_comp.get("Excitatory", 0) + off_comp.get("Inhibitory", 0)
log(f"\nscore_genes 真神经元占比: {n_neuron_sg}/{official.sum()} = {n_neuron_sg/official.sum()*100:.1f}% | 污染率: {(1-n_neuron_sg/official.sum())*100:.1f}%")

# 与 max-score 对比（加 Unknown 过滤，与 v5 口径一致）
m_ct_all = score_df.loc[official].idxmax(axis=1)
m_val_all = score_df.loc[official].max(axis=1)
m_ct_f = np.where(m_val_all > 0, m_ct_all, "Unknown")
m_neuron_n = ((m_ct_f == "Excitatory") | (m_ct_f == "Inhibitory")).sum()
m_unknown = (m_ct_f == "Unknown").sum()
log(f"max-score（含 Unknown 过滤）: 真神经元 {m_neuron_n}/{official.sum()} = {m_neuron_n/official.sum()*100:.1f}% | Unknown {m_unknown} | 污染率: {(1-m_neuron_n/official.sum())*100:.1f}%")

comp = pd.DataFrame({
    "method": ["max-score (v5)", "score_genes (AddModuleScore)"],
    "neuron_pct": [round(m_neuron_n/official.sum()*100, 1), round(n_neuron_sg/official.sum()*100, 1)],
    "pollution_pct": [round((1-m_neuron_n/official.sum())*100, 1), round((1-n_neuron_sg/official.sum())*100, 1)],
})
comp.to_csv(f"{OUT}/scorer_comparison.csv", index=False)
log("\n两打分器对比:\n" + comp.to_string(index=False))
log("\nDONE")
