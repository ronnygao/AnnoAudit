# -*- coding: utf-8 -*-
# v9b marker 稳健性分析：marker 列表留一法扰动，重算"官方谷氨酸能污染率"
# 目的：证明 97.9% 污染率不是 marker 列表选择的人为产物
import os, gzip, glob, time, itertools
import numpy as np, pandas as pd, scipy.io, scipy.sparse as sp

t0 = time.time()
RAW = "/Users/apple/WorkBuddy/2026-07-22-16-57-18/cerebri_data/raw"
ANNOT = "/Users/apple/Desktop/颅脑损伤专病数据库/CEREBRI_KCNC3生信分析/v3分析结果/cell_annotations.csv"
OUT = "/Users/apple/Desktop/颅脑损伤专病数据库/CEREBRI_KCNC3生信分析/v9_unbiased_clustering"
os.makedirs(OUT, exist_ok=True)

MARKERS_FULL = {
    "Excitatory": ["Slc17a7", "Rbfox3", "Neurod6", "Tbr1", "Camk2a"],
    "Inhibitory": ["Gad1", "Gad2", "Sst", "Pvalb", "Npy"],
    "Astrocyte": ["Gfap", "Aqp4", "Slc1a3", "Aldh1l1"],
    "Microglia": ["C1qb", "Tyrobp", "Cx3cr1", "Hexb", "Aif1"],
    "Oligodendrocyte": ["Mbp", "Plp1", "Mog", "Mag"],
    "OPC": ["Pdgfra", "Olig1", "Cspg4"],
    "Endothelial": ["Pecam1", "Cldn5", "Flt1"],
    "Pericyte": ["Pdgfrb", "Rgs5", "Vtn"],
}
def log(m):
    print(f"[{time.time()-t0:6.1f}s] {m}", flush=True)

# 读矩阵（443 基因子集：marker 基因足够）
ann = pd.read_csv(ANNOT)
ann["bc_plain"] = ann["barcode"].str.replace(r"_GSM\d+$", "", regex=True)
marker_genes = sorted({g for v in MARKERS_FULL.values() for g in v})
feat_files = sorted(os.listdir(RAW))
feat_path = next(f for f in feat_files if f.endswith("_features.tsv.gz"))
with gzip.open(os.path.join(RAW, feat_path), "rt") as f:
    features = [l.strip().split("\t") for l in f if not l.startswith("#")]
gene_names = [x[1] if len(x) > 1 else x[0] for x in features]
gene_to_idx = {g: i for i, g in enumerate(gene_names)}
need_idx = {g: gene_to_idx[g] for g in marker_genes if g in gene_to_idx}
need_cols = [need_idx[g] for g in marker_genes]
log(f"marker genes: {len(need_cols)}")

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
g2c = {g: i for i, g in enumerate(marker_genes)}
official = (meta["cell_type"] == "Glutamatergic_neuron").values
log(f"矩阵 {X.shape}, 官方谷氨酸能 {official.sum()}")

def annotate(markers):
    scores = {}
    for ct, genes in markers.items():
        idx = [g2c[g] for g in genes if g in g2c]
        scores[ct] = np.asarray(X_log[:, idx].mean(axis=1)).ravel() if idx else np.zeros(X_log.shape[0])
    sdf = pd.DataFrame(scores)
    mx = sdf.idxmax(axis=1); mv = sdf.max(axis=1)
    return np.where(mv > 0, mx, "Unknown")

# 基线（全 marker）
base = annotate(MARKERS_FULL)
base_poll = 1 - (np.isin(base[official], ["Excitatory", "Inhibitory"]).sum() / official.sum())
log(f"基线污染率: {base_poll*100:.2f}%")

# 留一扰动：每类去掉 1 个 marker（每类至少留 2 个）
results = []
for ct, genes in MARKERS_FULL.items():
    for drop_g in genes:
        if len(genes) <= 2: continue
        mk = {c: (gs if c != ct else [g for g in gs if g != drop_g]) for c, gs in MARKERS_FULL.items()}
        ann2 = annotate(mk)
        poll = 1 - (np.isin(ann2[official], ["Excitatory", "Inhibitory"]).sum() / official.sum())
        results.append({"perturbation": f"{ct}-{drop_g}", "pollution": poll*100})
        log(f"  {ct}-{drop_g}: 污染率 {poll*100:.2f}%")
res = pd.DataFrame(results)
res.loc[len(res)] = ["baseline", base_poll*100]
res.to_csv(os.path.join(OUT, "marker_robustness_leaveoneout.csv"), index=False)
log(f"\n污染率: 基线 {base_poll*100:.2f}% | 扰动均值 {res[res.perturbation!='baseline'].pollution.mean():.2f}% "
    f"± {res[res.perturbation!='baseline'].pollution.std():.2f}% | 范围 [{res[res.perturbation!='baseline'].pollution.min():.2f}, {res[res.perturbation!='baseline'].pollution.max():.2f}]")
log("DONE")
