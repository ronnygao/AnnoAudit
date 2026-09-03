#!/usr/bin/env python3
# ============================================================
# v10: 全基因 signal-provenance 指纹矩阵
# 1) 307 离子通道基因 + marker 基因在 9 个群体中的条件轨迹
#    （query = official_glut；候选 = 8 个 marker-defined 类型）
# 2) 逐基因 Spearman 指纹 → 基因 × 细胞类型 相关矩阵
# 3) 对照 query = marker_glut（真谷氨酸能）
# 4) 溯源统计：argmax 类型分布、家族聚合、Wilcoxon 检验
# ============================================================
import gzip
import glob
import json
import os
import time
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
from scipy.stats import spearmanr, wilcoxon

T0 = time.time()
def log(msg):
    print(f"[{time.time()-T0:7.1f}s] {msg}", flush=True)

RAW_DIR = "/Users/apple/Workbuddy/2026-07-22-16-57-18/cerebri_data/raw"
ANNOT_CSV = "/Users/apple/Desktop/颅脑损伤专病数据库/CEREBRI_KCNC3生信分析/v3分析结果/cell_annotations.csv"
ION_CSV = "/Users/apple/Desktop/颅脑损伤专病数据库/A1_离子通道组景观/ion_channel_landscape_summary.csv"
OUT_DIR = "/Users/apple/Desktop/颅脑损伤专病数据库/CEREBRI_KCNC3生信分析/v10_provenance"
os.makedirs(OUT_DIR, exist_ok=True)

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
CONDS = ["Naive_24h", "rCHI_24h", "CCI_ipsi_24h", "CCI_contra_24h",
         "CCI+HS_24h", "CCI_ipsi_7d", "CCI_ipsi_6mo"]

# ---------------- 加载 ----------------
ann = pd.read_csv(ANNOT_CSV)
ann["bc_plain"] = ann["barcode"].str.replace(r"_GSM\d+$", "", regex=True)

feat_files = sorted(os.listdir(RAW_DIR))
feat_path = next(f for f in feat_files if f.endswith("_features.tsv.gz"))
with gzip.open(os.path.join(RAW_DIR, feat_path), "rt") as f:
    features = [l.strip().split("\t") for l in f if not l.startswith("#")]
gene_names = [x[1] if len(x) > 1 else x[0] for x in features]
gene_to_idx = {g: i for i, g in enumerate(gene_names)}
log(f"features 基因数: {len(gene_names)}")

ion_sum = pd.read_csv(ION_CSV)
ion_genes = ion_sum["gene"].tolist()
need = sorted(set(m for v in MARKERS.values() for m in v) | set(ion_genes))
need_found = [g for g in need if g in gene_to_idx]
log(f"需用基因 {len(need)}，检出 {len(need_found)}")
need_cols = [gene_to_idx[g] for g in need_found]
g2c = {g: i for i, g in enumerate(need_found)}

X_parts, meta_parts = [], []
for gsm, grp in ann.groupby("gsm"):
    mtx_f = glob.glob(os.path.join(RAW_DIR, f"{gsm}_*_matrix.mtx.gz"))
    bar_f = glob.glob(os.path.join(RAW_DIR, f"{gsm}_*_barcodes.tsv.gz"))
    if not mtx_f or not bar_f:
        continue
    with gzip.open(mtx_f[0], "rb") as f:
        mat = scipy.io.mmread(f)
    mat = sp.csc_matrix(mat).T.tocsr()
    with gzip.open(bar_f[0], "rt") as f:
        barcodes = [l.strip() for l in f]
    bc2row = {bc: i for i, bc in enumerate(barcodes)}
    grp_rows = [bc2row[b] for b in grp["bc_plain"] if b in bc2row]
    if not grp_rows:
        continue
    sub = mat[grp_rows, :][:, need_cols].tocsr()
    keep = grp["bc_plain"].isin(barcodes)
    X_parts.append(sub)
    meta_parts.append(grp.loc[grp.index[keep]].copy())
    del mat, sub

X = sp.vstack(X_parts).tocsr()
meta = pd.concat(meta_parts).reset_index(drop=True)
log(f"合并矩阵: {X.shape}")

meta["total_umi"] = meta["total_umi"].astype(float)
qc = (meta["n_genes"] > 200) & (meta["mito_pct"] < 20)
meta = meta[qc].reset_index(drop=True)
X = X[qc.values]
umi = meta["total_umi"].values
X_cpm = X.multiply(1e6 / umi[:, None]).tocsr()
X_log = X_cpm.copy(); X_log.data = np.log1p(X_log.data)
log(f"QC 后: {X.shape[0]} 细胞")

# marker 注释（同 v5/v6）
scores = {}
for ct, genes in MARKERS.items():
    idx = [g2c[g] for g in genes if g in g2c]
    scores[ct] = np.asarray(X_log[:, idx].mean(axis=1)).ravel() if idx else np.zeros(X_log.shape[0])
score_df = pd.DataFrame(scores)
max_ct = score_df.idxmax(axis=1)
meta["celltype_marker"] = np.where(score_df.max(axis=1) > 0, max_ct, "Unknown")

official_glut = (meta["cell_type"] == "Glutamatergic_neuron").values
POPS = {}
for ct in MARKERS:
    POPS[ct] = (meta["celltype_marker"] == ct).values
POPS["official_glut"] = official_glut
log(f"marker 谷氨酸能: {int(POPS['Excitatory'].sum())} | 官方: {int(official_glut.sum())}")
for ct, m in POPS.items():
    log(f"  {ct}: {int(m.sum())}")

# 每群体每条件的细胞数（决定哪些条件-群体组合有效）
pop_cond_n = {}
for pname, mask in POPS.items():
    sub = meta[mask]
    pop_cond_n[pname] = {c: int((sub["condition"] == c).sum()) for c in CONDS}
pcn_df = pd.DataFrame(pop_cond_n).T
pcn_df.to_csv(os.path.join(OUT_DIR, "pop_condition_ncells.csv"))
log("\n群体×条件细胞数:\n" + pcn_df.to_string())

# ---------------- 轨迹矩阵 ----------------
# traj[(gene, pop)] = 7 条件的 mean log1p CPM
gene_list = [g for g in ion_genes if g in g2c] + [g for v in MARKERS.values() for g in v]
gene_list = sorted(set(gene_list))
log(f"\n轨迹基因数（离子通道+marker）: {len(gene_list)}")

traj_rows = []
for gi, g in enumerate(gene_list):
    xg = X_log[:, g2c[g]].toarray().ravel()
    row = {"gene": g}
    for pname, mask in POPS.items():
        sub_meta = meta[mask]
        sub_x = xg[mask]
        for c in CONDS:
            m = (sub_meta["condition"] == c).values
            row[f"{pname}::{c}"] = float(sub_x[m].mean()) if m.sum() > 0 else np.nan
    traj_rows.append(row)
    if (gi + 1) % 100 == 0:
        log(f"  轨迹进度 {gi+1}/{len(gene_list)}")
traj = pd.DataFrame(traj_rows)
traj.to_csv(os.path.join(OUT_DIR, "trajectories_all_genes.csv"), index=False)
log("轨迹矩阵已保存")

# ---------------- 指纹相关矩阵 ----------------
def fingerprint(query_pop, cand_pops, min_n=30):
    """对每基因：query 轨迹 vs 各候选类型轨迹的 Spearman（仅用双方都有数据的条件）"""
    rows = []
    for _, r in traj.iterrows():
        q = np.array([r[f"{query_pop}::{c}"] for c in CONDS], dtype=float)
        row = {"gene": r["gene"]}
        ok_q = ~np.isnan(q)
        for cp in cand_pops:
            t = np.array([r[f"{cp}::{c}"] for c in CONDS], dtype=float)
            ok = ok_q & ~np.isnan(t)
            if ok.sum() >= 5 and np.nanstd(q[ok]) > 0 and np.nanstd(t[ok]) > 0:
                rho, p = spearmanr(q[ok], t[ok])
                row[cp] = round(float(rho), 4)
                row[cp + "_p"] = round(float(p), 4)
                row[cp + "_n"] = int(ok.sum())
            else:
                row[cp] = np.nan
        rows.append(row)
    return pd.DataFrame(rows)

cand = list(MARKERS.keys())  # 8 个 marker-defined 类型
fp_official = fingerprint("official_glut", cand)
fp_official.to_csv(os.path.join(OUT_DIR, "fingerprint_query_official_glut.csv"), index=False)
fp_true = fingerprint("Excitatory", cand)
fp_true.to_csv(os.path.join(OUT_DIR, "fingerprint_query_true_glut.csv"), index=False)
log(f"指纹矩阵: official {fp_official.notna().mean(axis=0).round(2).to_dict()}")

# ---------------- 溯源统计 ----------------
def provenance_stats(fp, label):
    have = fp.dropna(subset=cand)
    argmax_type = have[cand].idxmax(axis=1)
    dist = argmax_type.value_counts()
    # 检出基因：离子通道且轨迹非恒定
    stats = {
        "query": label,
        "n_genes_total": len(fp),
        "n_genes_valid": len(have),
        "argmax_distribution": dist.to_dict(),
        "median_r_microglia": float(have["Microglia"].median()),
        "median_r_excitatory": float(have["Excitatory"].median()),
    }
    return stats, argmax_type

st_off, am_off = provenance_stats(fp_official, "official_glut")
st_true, am_true = provenance_stats(fp_true, "true_glut")
log(f"\nofficial 溯源: {st_off}")
log(f"true glut 溯源: {st_true}")

# Wilcoxon：official vs microglia 的 r 分布 vs true glut vs microglia
common = fp_official.dropna(subset=["Microglia"]).merge(
    fp_true[["gene", "Microglia"]], on="gene", suffixes=("_off", "_true")).dropna(subset=["Microglia_true"])
r_off = common["Microglia_off"].values
r_true = common["Microglia_true"].values
w_p = wilcoxon(r_off, r_true) if len(common) > 10 else (np.nan, np.nan)
w0_off = wilcoxon(r_off - 0) if len(r_off) > 10 else (np.nan, np.nan)
stats_out = {
    "official_provenance": st_off,
    "true_provenance": st_true,
    "wilcoxon_official_vs_true_microglia_r": {"n": len(common), "W": float(w_p[0]) if not np.isnan(w_p[0]) else None, "p": float(w_p[1])},
    "wilcoxon_official_r_vs_zero": {"n": len(r_off), "p": float(w0_off[1])},
    "mean_r_official_vs_microglia": float(np.mean(r_off)),
    "mean_r_true_vs_microglia": float(np.mean(r_true)),
}
with open(os.path.join(OUT_DIR, "provenance_stats.json"), "w") as f:
    json.dump(stats_out, f, indent=2, default=str)
log(f"\nWilcoxon official vs true (microglia r): n={len(common)}, p={w_p[1]:.3e}")
log(f"mean r official-vs-microglia = {np.mean(r_off):.3f} | true-vs-microglia = {np.mean(r_true):.3f}")

# ---------------- 家族聚合 ----------------
fam = ion_sum[["gene", "family"]].dropna(subset=["gene"])
def fam_prefix(g):
    import re
    rules = [
        (r"^Kctd", "KCTD"), (r"^Kcne", "Kcne"), (r"^Kcnip", "Kcnip"), (r"^Kcnab", "Kcnab"),
        (r"^Kcn[ns]", "Kcn voltage-gated"), (r"^Kcn[mt]", "Ca-activated K"), (r"^Kcnj", "Inward-rectifier K"),
        (r"^Kcnk", "Two-pore K"), (r"^Scn", "Scn Na"), (r"^Cacna", "Cacna"), (r"^Cacnb", "Cacnb"), (r"^Cacng", "Cacng"),
        (r"^Clcn", "ClC"), (r"^Clca", "ClCA"), (r"^Best", "Bestrophin"), (r"^Ano", "Anoctamin"),
        (r"^Trp", "TRP"), (r"^Mcoln", "TRP-ML"), (r"^Pkd2", "TRPP"), (r"^Hcn", "HCN"), (r"^Cng", "CNG"),
        (r"^Ry", "Ryanodine"), (r"^Itpr", "IP3R"), (r"^Gabr", "GABA-A"), (r"^Glra", "Glycine"), (r"^Glrb", "Glycine"),
        (r"^Gria", "iGluR AMPA"), (r"^Grin", "iGluR NMDA"), (r"^Grik", "iGluR kainate"),
        (r"^Chrn", "Nicotinic"), (r"^Htr3", "5-HT3"), (r"^P2rx", "P2X"),
        (r"^Tmem38", "TRIC"), (r"^Tpcn", "Two-pore TPC"), (r"^Asic", "ASIC"), (r"^Scnn", "ENaC"),
    ]
    for pat, name in rules:
        if re.match(pat, g):
            return name
    return "Other"
fam["family"] = fam["family"].fillna(fam["gene"].map(fam_prefix))
fam_map = dict(zip(fam["gene"], fam["family"]))

for fp, name in [(fp_official, "official"), (fp_true, "true")]:
    fp2 = fp.copy()
    fp2["family"] = fp2["gene"].map(fam_map)
    fam_agg = fp2.dropna(subset=cand).groupby("family")[cand].agg(["median", "count"])
    fam_agg.to_csv(os.path.join(OUT_DIR, f"family_fingerprint_{name}.csv"))
    log(f"\n{name} 家族×类型 指纹（median r）:\n{fam_agg['median'].round(2).to_string()}")

log("\nV10 DONE")
