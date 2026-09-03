#!/usr/bin/env python3
# ============================================================
# v15_tau_sensitivity: C2 margin gate tau 敏感性分析
# 复用 v14 终版协议的 C2 逻辑（scanpy score_genes + margin = glia_max - neu_max）
# tau 只作用于最终阈值 (margin > tau).mean()，因此 margin 分布只需算一次；
# C1/C3 直接读取 v14 终版 final_audit.json（不重跑 KMeans/CellTypist）。
# 输出: tau_sensitivity.json（tau=1.0/1.25/1.5/1.75/2.0 的 C2 与 ACS）
#        margin_distribution.csv（margin 值，供复现/绘图）
# 用法: python3 v15_tau_sensitivity.py --dataset cerebri|gse330130
# ============================================================
import os, sys, glob, gzip, time, json, warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import scipy.io, scipy.sparse as sp
import anndata as ad
import scanpy as sc

T0 = time.time()
def log(m): print(f"[{time.time()-T0:7.1f}s] {m}", flush=True)

BASE = "/Users/apple/Desktop/颅脑损伤专病数据库/CEREBRI_KCNC3生信分析"
TAUS = [1.0, 1.25, 1.5, 1.75, 2.0]

MARKERS_MOUSE = {
    "Excitatory": ["Slc17a7","Rbfox3","Neurod6","Tbr1","Camk2a"],
    "Inhibitory": ["Gad1","Gad2","Sst","Pvalb","Npy"],
    "Astrocyte": ["Gfap","Aqp4","Slc1a3","Aldh1l1"],
    "Microglia": ["C1qb","Tyrobp","Cx3cr1","Hexb","Aif1"],
    "Oligodendrocyte": ["Mbp","Plp1","Mog","Mag"],
    "OPC": ["Pdgfra","Olig1","Cspg4"],
    "Endothelial": ["Pecam1","Cldn5","Flt1"],
    "Pericyte": ["Pdgfrb","Rgs5","Vtn"],
}
MARKERS_HUMAN = {
    "Excitatory": ["SLC17A7","RBFOX3","NEUROD6","TBR1","CAMK2A"],
    "Inhibitory": ["GAD1","GAD2","SST","PVALB","NPY"],
    "Astrocyte": ["GFAP","AQP4","SLC1A3","ALDH1L1"],
    "Microglia": ["C1QB","TYROBP","CX3CR1","HEXB","AIF1"],
    "Oligodendrocyte": ["MBP","PLP1","MOG","MAG"],
    "OPC": ["PDGFRA","OLIG1","CSPG4"],
    "Endothelial": ["PECAM1","CLDN5","FLT1"],
    "Pericyte": ["PDGFRB","RGS5","VTN"],
}
NEU_TYPES = ["Excitatory", "Inhibitory"]

def margin_from_scores(S, glial_types):
    neu_max = S[[f"mod_{t}" for t in NEU_TYPES]].max(axis=1).values
    glia_max = S[[f"mod_{t}" for t in glial_types]].max(axis=1).values
    return glia_max - neu_max

def tau_table(margin):
    return {str(t): float((margin > t).mean() * 100) for t in TAUS}

# ------------------------------------------------------------
def run_cerebri():
    RAW_DIR = "/Users/apple/Workbuddy/2026-07-22-16-57-18/cerebri_data/raw"
    ANNOT_CSV = os.path.join(BASE, "v3分析结果/cell_annotations.csv")
    V14_JSON = os.path.join(BASE, "v14_cerebri_final/final_audit.json")
    OUT_DIR = os.path.join(BASE, "v15_tau_sensitivity")
    os.makedirs(OUT_DIR, exist_ok=True)

    ann = pd.read_csv(ANNOT_CSV)
    ann["bc_plain"] = ann["barcode"].str.replace(r"_GSM\d+$", "", regex=True)
    sel = ann[ann["cell_type"].astype(str) == "Glutamatergic_neuron"].copy()
    qc = (sel["n_genes"] > 200) & (sel["mito_pct"] < 20)
    sel = sel[qc].reset_index(drop=True)
    log(f"官方 Glutamatergic_neuron QC 后: {len(sel)}")

    feat_files = sorted(os.listdir(RAW_DIR))
    feat_path = next(f for f in feat_files if f.endswith("_features.tsv.gz"))
    with gzip.open(os.path.join(RAW_DIR, feat_path), "rt") as f:
        features = [l.strip().split("\t") for l in f if not l.startswith("#")]
    gene_names = [x[1] if len(x) > 1 else x[0] for x in features]
    n_g = len(gene_names)

    X_parts = []
    for gi, (gsm, grp) in enumerate(sel.groupby("gsm")):
        mtx_f = glob.glob(os.path.join(RAW_DIR, f"{gsm}_*_matrix.mtx.gz"))
        bar_f = glob.glob(os.path.join(RAW_DIR, f"{gsm}_*_barcodes.tsv.gz"))
        with gzip.open(mtx_f[0], "rb") as f:
            mat = scipy.io.mmread(f)
        mat = sp.csc_matrix(mat).T.tocsr()
        with gzip.open(bar_f[0], "rt") as f:
            barcodes = [l.strip() for l in f]
        bc2row = {bc: i for i, bc in enumerate(barcodes)}
        rows = [bc2row[b] for b in grp["bc_plain"]]
        X_parts.append(mat[rows, :].astype(np.float32).tocsr())
        log(f"[{gi+1}] {gsm}: {len(rows)} 细胞")
    X = sp.vstack(X_parts).tocsr()
    del X_parts

    sym = pd.Index(gene_names)
    dup = sym.duplicated(keep=False)
    if dup.any():
        dup_pos = np.where(dup)[0]
        sums = np.asarray(X[:, dup].sum(axis=0)).ravel()
        order = pd.Series(sums, index=dup_pos).groupby(sym[dup].values).idxmax()
        keep_cols = np.zeros(n_g, dtype=bool)
        keep_cols[np.asarray(order.values, dtype=int)] = True
        keep_cols[~dup] = True
        X = X[:, keep_cols].tocsr()
        genes = sym[keep_cols].values
    else:
        genes = sym.values
    log(f"全基因矩阵 {X.shape}")

    umi = np.asarray(X.sum(axis=1)).ravel()
    X_cpm = X.multiply(1e6 / np.maximum(umi, 1)[:, None]).tocsr()
    X_cpm.data = np.log1p(X_cpm.data)

    A = ad.AnnData(X=X_cpm)
    A.obs_names = [f"c{i}" for i in range(X.shape[0])]
    A.var_names = list(genes)
    np.random.seed(0)
    for ct, gs in MARKERS_MOUSE.items():
        use = [g for g in gs if g in A.var_names]
        sc.tl.score_genes(A, gene_list=use, score_name=f"mod_{ct}")
    S = A.obs[[c for c in A.obs.columns if c.startswith("mod_")]]
    margin = margin_from_scores(S, [t for t in MARKERS_MOUSE if t not in NEU_TYPES])
    pd.DataFrame({"margin_glia_minus_neu": margin}).to_csv(
        os.path.join(OUT_DIR, "cerebri_margin_distribution.csv"), index=False)
    c2_tau = tau_table(margin)
    log(f"CEREBRI C2 by tau: {c2_tau}")

    v14 = json.load(open(V14_JSON))
    c1, c3 = v14["checks"]["C1"], v14["checks"]["C3"]
    out = {"dataset": "CEREBRI mouse TBI atlas (GSE269748)",
           "official_subset": "Glutamatergic_neuron (n=45051, QC-passed)",
           "n_cells": int(X.shape[0]),
           "C1_pct": c1, "C3_pct": c3,
           "C2_pct_by_tau": {k: round(v, 1) for k, v in c2_tau.items()},
           "ACS_pct_by_tau": {k: round(np.mean([c1, v, c3]), 1) for k, v in c2_tau.items()},
           "baseline_tau": 1.5,
           "baseline_ACS_pct": v14["ACS_pct"]}
    with open(os.path.join(OUT_DIR, "cerebri_tau_sensitivity.json"), "w") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    log("CEREBRI tau 敏感性: " + json.dumps(out["C2_pct_by_tau"]) + " ACS: " + json.dumps(out["ACS_pct_by_tau"]))

# ------------------------------------------------------------
def run_gse330130():
    H5 = "/Users/apple/Desktop/TBI分子时钟/数据/GSE330130_LSC_SNRNASEQ_COUNTS.h5ad"
    V14_JSON = os.path.join(BASE, "v14_gse330130_final/final_audit.json")
    OUT_DIR = os.path.join(BASE, "v15_tau_sensitivity")
    os.makedirs(OUT_DIR, exist_ok=True)
    NEURON_LABELS = ["Inhibitory Neuron", "Excitatory  Neuron", "Motor Neurons"]

    a = ad.read_h5ad(H5)
    log(f"loaded {a.shape}")
    a.var["symbol"] = a.var["common_name"].astype(str)
    a.var_names = a.var["symbol"].values
    a = a[:, ~pd.Index(a.var_names).duplicated(keep="first")].copy()
    X = a.X.tocsr()
    umi = np.asarray(X.sum(axis=1)).ravel()
    Xn = X.multiply(1e6 / np.maximum(umi, 1)[:, None]).tocsr()
    Xn.data = np.log1p(Xn.data)
    a.X = Xn
    del X, Xn

    off = a.obs["Cell_Type"].isin(NEURON_LABELS).values
    log(f"official neurons: {off.sum()}")
    sub = a[off].copy()
    np.random.seed(0)
    for ct, genes in MARKERS_HUMAN.items():
        use = [g for g in genes if g in sub.var_names]
        sc.tl.score_genes(sub, gene_list=use, score_name=f"mod_{ct}")
    S = sub.obs[[c for c in sub.obs.columns if c.startswith("mod_")]]
    margin = margin_from_scores(S, [t for t in MARKERS_HUMAN if t not in NEU_TYPES])
    pd.DataFrame({"margin_glia_minus_neu": margin}).to_csv(
        os.path.join(OUT_DIR, "gse330130_margin_distribution.csv"), index=False)
    c2_tau = tau_table(margin)
    log(f"GSE330130 C2 by tau: {c2_tau}")

    v14 = json.load(open(V14_JSON))
    c1, c3 = v14["checks"]["C1"], v14["checks"]["C3"]
    out = {"dataset": "GSE330130 human ALS lumbar spinal cord snRNA-seq",
           "official_subset": "Neuron (n=2450)",
           "n_cells": int(sub.n_obs),
           "C1_pct": c1, "C3_pct": c3,
           "C2_pct_by_tau": {k: round(v, 1) for k, v in c2_tau.items()},
           "ACS_pct_by_tau": {k: round(np.mean([c1, v, c3]), 1) for k, v in c2_tau.items()},
           "baseline_tau": 1.5,
           "baseline_ACS_pct": v14["ACS_pct"]}
    with open(os.path.join(OUT_DIR, "gse330130_tau_sensitivity.json"), "w") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    log("GSE330130 tau 敏感性: " + json.dumps(out["C2_pct_by_tau"]) + " ACS: " + json.dumps(out["ACS_pct_by_tau"]))

if __name__ == "__main__":
    ds = sys.argv[1] if len(sys.argv) > 1 else "both"
    if ds in ("cerebri", "both"):
        run_cerebri()
    if ds in ("gse330130", "both"):
        run_gse330130()
    log("V15 DONE")
