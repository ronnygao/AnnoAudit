#!/usr/bin/env python3
# ============================================================
# v14_cerebri_final: CEREBRI 官方 glutamatergic 子集终版审计
# 协议与 v12 模拟基准完全一致：
#   C1 marker max-score（log1p CPM, 8 类 argmax, 阈值>0）
#   C2 模块评分 margin 门控（scanpy score_genes + τ=1.5）
#   C3 无监督聚类（top-2000 变异基因 → PCA30 → KMeans k=10 → 簇 majority 映射）
#   C4 适用性门控 CellTypist（适用性表沿用 v12_simulation：
#      DMB 12.4% / MTG 18.1% / PFC 1.6%，全部 <80% → N/A）
#   ACS = 可用检查均值（本数据 = mean(C1, C2, C3)）
# ============================================================
import os, glob, gzip, pickle, time, json, warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import scipy.io, scipy.sparse as sp
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
import scanpy as sc
import anndata as ad

T0 = time.time()
def log(m): print(f"[{time.time()-T0:7.1f}s] {m}", flush=True)

RAW_DIR = "/Users/apple/Workbuddy/2026-07-22-16-57-18/cerebri_data/raw"
ANNOT_CSV = "/Users/apple/Desktop/颅脑损伤专病数据库/CEREBRI_KCNC3生信分析/v3分析结果/cell_annotations.csv"
OUT_DIR = "/Users/apple/Desktop/颅脑损伤专病数据库/CEREBRI_KCNC3生信分析/v14_cerebri_final"
MODEL_DIR = "/Users/apple/.celltypist/data/models"
APPL_CSV = "/Users/apple/Desktop/颅脑损伤专病数据库/CEREBRI_KCNC3生信分析/v12_simulation/applicability_table.csv"
TAU = 1.5
os.makedirs(OUT_DIR, exist_ok=True)

MARKERS = {
    "Excitatory": ["Slc17a7","Rbfox3","Neurod6","Tbr1","Camk2a"],
    "Inhibitory": ["Gad1","Gad2","Sst","Pvalb","Npy"],
    "Astrocyte": ["Gfap","Aqp4","Slc1a3","Aldh1l1"],
    "Microglia": ["C1qb","Tyrobp","Cx3cr1","Hexb","Aif1"],
    "Oligodendrocyte": ["Mbp","Plp1","Mog","Mag"],
    "OPC": ["Pdgfra","Olig1","Cspg4"],
    "Endothelial": ["Pecam1","Cldn5","Flt1"],
    "Pericyte": ["Pdgfrb","Rgs5","Vtn"],
}
NEU_TYPES = ["Excitatory", "Inhibitory"]
GLIAL_TYPES = [t for t in MARKERS if t not in NEU_TYPES]

# ---------- 选细胞：官方 Glutamatergic_neuron + QC ----------
ann = pd.read_csv(ANNOT_CSV)
ann["bc_plain"] = ann["barcode"].str.replace(r"_GSM\d+$", "", regex=True)
off_all = ann["cell_type"].astype(str) == "Glutamatergic_neuron"
log(f"官方 Glutamatergic_neuron 总数: {off_all.sum()}")
sel = ann[off_all].copy()
qc = (sel["n_genes"] > 200) & (sel["mito_pct"] < 20)
sel = sel[qc].reset_index(drop=True)
log(f"QC 后审计子集: {len(sel)}")

# ---------- 全基因加载 ----------
feat_files = sorted(os.listdir(RAW_DIR))
feat_path = next(f for f in feat_files if f.endswith("_features.tsv.gz"))
with gzip.open(os.path.join(RAW_DIR, feat_path), "rt") as f:
    features = [l.strip().split("\t") for l in f if not l.startswith("#")]
gene_names = [x[1] if len(x) > 1 else x[0] for x in features]
n_g = len(gene_names)

X_parts = []
ngsm = sel["gsm"].nunique()
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
    del mat
    log(f"[{gi+1}/{ngsm}] {gsm}: {len(rows)} 细胞")
X = sp.vstack(X_parts).tocsr()
del X_parts
log(f"全基因矩阵: {X.shape}")

# symbol 去重（保留总表达最高）
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
    log(f"去重: {X.shape[1]} 基因")
else:
    genes = sym.values
gene_pos = {g: i for i, g in enumerate(genes)}

# ---------- 归一化 ----------
umi = np.asarray(X.sum(axis=1)).ravel()
X_cpm = X.multiply(1e6 / np.maximum(umi, 1)[:, None]).tocsr()
X_cpm.data = np.log1p(X_cpm.data)
X_10k = X.multiply(1e4 / np.maximum(umi, 1)[:, None]).tocsr()
X_10k.data = np.log1p(X_10k.data)
log("归一化完成（log1p CPM + log1p CP10K）")

mk_idx = {ct: np.array([gene_pos[g] for g in gs if g in gene_pos]) for ct, gs in MARKERS.items()}
missing = [(ct, g) for ct, gs in MARKERS.items() for g in gs if g not in gene_pos]
log(f"marker 缺失: {missing if missing else '无'}")

# ---------- C1: marker max-score ----------
n = X.shape[0]
s1 = {}
for ct, gi in mk_idx.items():
    s1[ct] = np.asarray(X_cpm[:, gi].mean(axis=1)).ravel() if len(gi) else np.zeros(n)
s1_df = pd.DataFrame(s1)
bval = s1_df.max(axis=1).values
best = s1_df.values.argmax(axis=1)
types = list(MARKERS.keys())
mtype = np.where(bval > 0, np.array(types)[best], "Unknown")
c1 = 1 - pd.Series(mtype).isin(NEU_TYPES).mean()
log(f"C1 marker max-score 非神经元率 = {c1*100:.1f}%")
log("marker 定型分布:\n" + pd.Series(mtype).value_counts().to_string())

# ---------- C2: margin 门控模块评分 ----------
A = ad.AnnData(X=X_cpm)
A.obs_names = [f"c{i}" for i in range(n)]
A.var_names = list(genes)
np.random.seed(0)
for ct, gs in MARKERS.items():
    use = [g for g in gs if g in A.var_names]
    sc.tl.score_genes(A, gene_list=use, score_name=f"mod_{ct}")
mods = [c for c in A.obs.columns if c.startswith("mod_")]
S = A.obs[mods]
neu_max = S[[f"mod_{t}" for t in NEU_TYPES]].max(axis=1).values
glia_max = S[[f"mod_{t}" for t in GLIAL_TYPES]].max(axis=1).values
margin = glia_max - neu_max
c2 = float((margin > TAU).mean())
amax = S.values.argmax(axis=1)
amod = np.array([m.replace("mod_", "") for m in S.columns])[amax]
c2_argmax = 1 - pd.Series(amod).isin(NEU_TYPES).mean()
log(f"C2 margin 门控 = {c2*100:.1f}% (无门控 argmax = {c2_argmax*100:.1f}%)")

# ---------- C3: 无监督聚类 ----------
gm = np.asarray(X_cpm.mean(axis=0)).ravel()
gv = np.asarray(X_cpm.multiply(X_cpm).mean(axis=0)).ravel() - gm ** 2
keepvar = gm > 0.01
top = np.argsort(np.where(keepvar, gv, -1))[-min(2000, keepvar.sum()):]
Xh = np.asarray(X_cpm[:, top].todense(), dtype=np.float32)
mu, sd = Xh.mean(0), Xh.std(0) + 1e-9
Xs_ = np.clip((Xh - mu) / sd, -10, 10)
del Xh
Xp_ = PCA(n_components=30, random_state=0, svd_solver="randomized").fit_transform(Xs_)
del Xs_
cl = KMeans(n_clusters=10, random_state=0, n_init=10).fit_predict(Xp_)
cl_type = pd.Series(mtype).groupby(cl).agg(lambda s: s.value_counts().index[0])
log("簇 majority 类型:\n" + cl_type.to_string())
ctype = pd.Series(cl).map(cl_type).values
c3 = 1 - pd.Series(ctype).isin(NEU_TYPES).mean()
log(f"C3 KMeans 非神经元率 = {c3*100:.1f}%")
pd.DataFrame({"cluster": cl, "cluster_type": ctype, "marker_type": mtype}).to_csv(
    os.path.join(OUT_DIR, "c3_clusters.csv"), index=False)

# ---------- C4: 适用性门控（原始值仅供参考） ----------
APPLIC = pd.read_csv(APPL_CSV)
log("适用性表 (v12):\n" + APPLIC.to_string(index=False))

def is_neuronal_dmb(lab):
    return "Neuron" in str(lab) or "Excitatory" in str(lab) or "Inhibitory" in str(lab)

with open(os.path.join(MODEL_DIR, "Developing_Mouse_Brain.pkl"), "rb") as f:
    d = pickle.load(f)
lr, scl = d["Model"], d["Scaler_"]
feats = [str(x) for x in lr.features]
idx = [gene_pos.get(g) for g in feats]
F = np.zeros((n, len(feats)), dtype=np.float64)
X10 = X_10k.tocsr()
for j, i in enumerate(idx):
    if i is not None:
        F[:, j] = np.asarray(X10[:, i].todense()).ravel()
lab = lr.predict(scl.transform(F))
pd.Series(lab).value_counts().to_csv(os.path.join(OUT_DIR, "celltypist_DMB_official.csv"))
c4_raw = 1 - np.mean([is_neuronal_dmb(l) for l in lab])
appl_dmb = bool(APPLIC.loc[APPLIC.model == "Developing_Mouse_Brain", "applicable_exc"].iloc[0])
nmatch = sum(i is not None for i in idx)
log(f"C4 DMB: 特征匹配 {nmatch}/{len(feats)}, 原始非神经元率 {c4_raw*100:.1f}%, 适用性 {'PASS' if appl_dmb else 'FAIL'}")
c4 = c4_raw if appl_dmb else np.nan

# ---------- ACS ----------
checks = {"C1": round(c1*100,1), "C2": round(c2*100,1), "C3": round(c3*100,1),
          "C4": round(c4*100,1) if not np.isnan(c4) else None}
vals = [c1*100, c2*100, c3*100] + ([c4*100] if not np.isnan(c4) else [])
acs = float(np.mean(vals))
final = {
    "dataset": "CEREBRI mouse TBI atlas (GSE269748)",
    "official_label": f"Glutamatergic_neuron (n={int(off_all.sum())}, QC-passed n={n})",
    "checks": checks,
    "C2_ungated_argmax_pct": round(c2_argmax*100,1),
    "C4_raw_DMB_pct": round(c4_raw*100,1),
    "C4_applicable_models": [],
    "C4_excluded_models": ["Developing_Mouse_Brain (applicability 12.4%)",
                            "Adult_Human_MTG (cross-species, 18.1%)",
                            "Adult_Human_PrefrontalCortex (cross-species, 1.6%)"],
    "ACS_pct": round(acs,1),
    "protocol": "C1 argmax(>0); C2 scanpy score_genes + margin tau=1.5; C3 top2000HVG-PCA30-KMeans10; C4 applicability-gated CellTypist (all models failed gate)",
}
with open(os.path.join(OUT_DIR, "final_audit.json"), "w") as f:
    json.dump(final, f, indent=2, ensure_ascii=False)
log("最终审计: " + json.dumps(final, ensure_ascii=False))
log("V14 CEREBRI DONE")
