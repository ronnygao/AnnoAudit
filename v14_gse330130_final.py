#!/usr/bin/env python3
# v14: GSE330130 终版审计（与 v12 模拟协议完全一致）
#   [2026-08-31 修订] panel 口径统一：20 基因小 panel → 32 基因严格人类直系同源 panel
#     （与 v14_cerebri_final 小鼠 32 基因 panel 及手稿 "32 marker genes" 一致；不含 SLC17A6）
#     C1 阈值 0.1 → >0（与 v12/v14_cerebri 协议一致）
#     C3 不再沿用 v11 硬编码 100.0（v11 majority 映射基于旧 panel），用新 panel 重算
#   C1 marker max-score（log1p CPM, argmax, 阈值 >0）
#   C2 模块评分 margin 门控（scanpy score_genes + τ=1.5）
#   C3 无监督聚类（60k 子采样 → top-2000 HVG → PCA50 → KMeans k=20 → 簇 majority 映射）
#   C4 适用性门控 CellTypist（MTG/PFC 用 common_name 匹配；DMB 跨物种排除）
#       适用性: marker 定型参考池上目标类型(神经元)准确率 ≥80%
#   ACS = 可用检查的均值
import os, time, pickle, warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import scanpy as sc
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA

T0 = time.time()
def log(m): print(f"[{time.time()-T0:7.1f}s] {m}", flush=True)

H5 = "/Users/apple/Desktop/TBI分子时钟/数据/GSE330130_LSC_SNRNASEQ_COUNTS.h5ad"
OUT_DIR = "/Users/apple/Desktop/颅脑损伤专病数据库/CEREBRI_KCNC3生信分析/v14_gse330130_final"
MODEL_DIR = "/Users/apple/.celltypist/data/models"
TAU = 1.5
os.makedirs(OUT_DIR, exist_ok=True)

# 32 基因严格人类直系同源 panel（= 小鼠 32 基因大 panel 的直系同源，不含 SLC17A6）
MARKERS_H = {
    "Excitatory": ["SLC17A7", "RBFOX3", "NEUROD6", "TBR1", "CAMK2A"],
    "Inhibitory": ["GAD1", "GAD2", "SST", "PVALB", "NPY"],
    "Astrocyte": ["GFAP", "AQP4", "SLC1A3", "ALDH1L1"],
    "Microglia": ["C1QB", "TYROBP", "CX3CR1", "HEXB", "AIF1"],
    "Oligodendrocyte": ["MBP", "PLP1", "MOG", "MAG"],
    "OPC": ["PDGFRA", "OLIG1", "CSPG4"],
    "Endothelial": ["PECAM1", "CLDN5", "FLT1"],
    "Pericyte": ["PDGFRB", "RGS5", "VTN"],
}
NGENE = sum(len(v) for v in MARKERS_H.values())
assert NGENE == 32, NGENE
NEU_TYPES = ["Excitatory", "Inhibitory"]
NEURON_LABELS = ["Inhibitory Neuron", "Excitatory  Neuron", "Motor Neurons"]

# ---------- 加载 ----------
a = ad.read_h5ad(H5)
log(f"loaded {a.shape}")
a.var["symbol"] = a.var["common_name"].astype(str)
a.var_names = a.var["symbol"].values
a = a[:, ~pd.Index(a.var_names).duplicated(keep="first")].copy()
log(f"dedup 后 {a.shape}")
missing = [g for gs in MARKERS_H.values() for g in gs if g not in a.var_names]
log(f"panel 基因缺失: {missing if missing else '无（32/32 全部命中）'}")

# ---------- 归一化 log1p CPM ----------
X = a.X.tocsr()
umi = np.asarray(X.sum(axis=1)).ravel()
Xn = X.multiply(1e6 / np.maximum(umi, 1)[:, None]).tocsr()
Xn.data = np.log1p(Xn.data)
a.X = Xn
del X, Xn
log("归一化完成")

# ---------- marker 定型（阈值 >0，与 v12/CEREBRI 一致） ----------
gidx = {g: a.var_names.get_loc(g) for ms in MARKERS_H.values() for g in ms if g in a.var_names}
scores = {}
for ct, genes in MARKERS_H.items():
    idx = [gidx[g] for g in genes if g in gidx]
    scores[ct] = np.asarray(a.X[:, idx].mean(axis=1)).ravel() if idx else np.zeros(a.n_obs)
sdf = pd.DataFrame(scores)
best = sdf.idxmax(axis=1)
bval = sdf.max(axis=1)
a.obs["marker_type"] = np.where(bval > 0, best.values, "Unknown")
log("marker 定型分布:\n" + a.obs["marker_type"].value_counts().to_string())

off = a.obs["Cell_Type"].isin(NEURON_LABELS).values
log(f"official neurons: {off.sum()}")
ct1 = pd.crosstab(a.obs.loc[off, "Cell_Type"], a.obs.loc[off, "marker_type"])
ct1.to_csv(os.path.join(OUT_DIR, "official_neuron_marker_crosstab.csv"))

# ---------- C1 ----------
m_neu = a.obs["marker_type"].isin(NEU_TYPES).values
c1 = 1 - (off & m_neu).sum() / off.sum()
log(f"C1 marker max-score 非神经元率 = {c1*100:.1f}%")

# ---------- C2 (margin 门控) ----------
sub = a[off].copy()
np.random.seed(0)
for ct, genes in MARKERS_H.items():
    use = [g for g in genes if g in sub.var_names]
    sc.tl.score_genes(sub, gene_list=use, score_name=f"mod_{ct}")
mods = [c for c in sub.obs.columns if c.startswith("mod_")]
S = sub.obs[mods]
neu_max = S[[f"mod_{t}" for t in NEU_TYPES]].max(axis=1).values
glial_types = [t for t in MARKERS_H if t not in NEU_TYPES]
glia_max = S[[f"mod_{t}" for t in glial_types]].max(axis=1).values
margin = glia_max - neu_max
c2 = float((margin > TAU).mean())
amod = S.values.argmax(axis=1)
c2_argmax = 1 - pd.Series(np.array([m.replace("mod_","") for m in S.columns])[amod]).isin(NEU_TYPES).mean()
log(f"C2 margin 门控 = {c2*100:.1f}% (无门控 argmax = {c2_argmax*100:.1f}%)")
del sub

# ---------- C3: 无监督 KMeans（60k 子采样，与 v11 同管线，新 panel majority 映射） ----------
rng = np.random.RandomState(0)
n_sub = min(60000, a.n_obs)
sub_idx = rng.choice(a.n_obs, n_sub, replace=False)
b = a[sub_idx].copy()
Xb = b.X.tocsr()
gene_mean = np.asarray(Xb.mean(axis=0)).ravel()
gene_var = np.asarray(Xb.multiply(Xb).mean(axis=0)).ravel() - gene_mean ** 2
keep = gene_mean > 0.01
top_idx = np.argsort(np.where(keep, gene_var, -1))[-min(2000, keep.sum()):]
Xh = np.asarray(Xb[:, top_idx].todense(), dtype=np.float32)
del Xb
mu = Xh.mean(axis=0); sd = Xh.std(axis=0) + 1e-9
Xs = (Xh - mu) / np.where(sd > 0, sd, 1)
Xs = np.clip(Xs, -10, 10)
del Xh
pca = PCA(n_components=50, random_state=0, svd_solver="randomized")
Xp = pca.fit_transform(Xs)
del Xs
km = KMeans(n_clusters=20, random_state=0, n_init=10)
cl = km.fit_predict(Xp).astype(str)
del Xp
cl_map = pd.Series(b.obs["marker_type"].values).groupby(cl).agg(lambda s: s.value_counts().index[0])
cl_type = pd.Series(cl).map(cl_map).values
neu_b = pd.Series(cl_type).isin(NEU_TYPES).values
off_sub = off[sub_idx]
c3 = 1 - (off_sub & neu_b).sum() / max(off_sub.sum(), 1)
log(f"C3 KMeans 非神经元率 = {c3*100:.1f}% (子采样 n={n_sub}, 官方神经元 {off_sub.sum()})")
ct3 = pd.crosstab(b.obs.loc[off_sub, "Cell_Type"], pd.Series(cl_type[off_sub], index=b.obs.index[off_sub]))
ct3.to_csv(os.path.join(OUT_DIR, "official_neuron_kmeans_crosstab.csv"))
del b

# ---------- C4: 适用性 + 官方子集注释 ----------
sym_pos = {g: i for i, g in enumerate(a.var_names)}

def predict_rows(rows, mf):
    """从磁盘重新加载 raw counts 指定行, 构建 CP10K log1p 特征并预测"""
    with open(os.path.join(MODEL_DIR, mf + ".pkl"), "rb") as f:
        d = pickle.load(f)
    lr, scl = d["Model"], d["Scaler_"]
    feats = [str(x) for x in lr.features]
    raw = ad.read_h5ad(H5)
    raw.var["symbol"] = raw.var["common_name"].astype(str)
    raw.var_names = raw.var["symbol"].values
    raw = raw[:, ~pd.Index(raw.var_names).duplicated(keep="first")].copy()
    Xr = raw.X[rows].tocsr()
    u = np.asarray(Xr.sum(axis=1)).ravel()
    X10 = Xr.multiply(1e4 / np.maximum(u, 1)[:, None]).tocsr()
    X10.data = np.log1p(X10.data)
    idx = [sym_pos.get(g) for g in feats]
    F = np.zeros((X10.shape[0], len(feats)), dtype=np.float64)
    for j, i in enumerate(idx):
        if i is not None:
            F[:, j] = np.asarray(X10[:, i].todense()).ravel()
    lab = lr.predict(scl.transform(F))
    return lab, sum(i is not None for i in idx), len(feats)

# 参考池（marker 定型）
rr = np.random.RandomState(0)
refs = {}
for ct in ["Excitatory", "Inhibitory", "Oligodendrocyte", "Astrocyte", "Microglia"]:
    pool = np.where((a.obs["marker_type"] == ct).values)[0]
    if len(pool) >= 50:
        refs[ct] = rr.choice(pool, min(1000, len(pool)), replace=False)
log("参考池: " + ", ".join(f"{k}={len(v)}" for k, v in refs.items()))

is_neuronal_mtg = lambda l: (not any(k in str(l) for k in
    ("Oligo","Astro","OPC","Micro","Endo","VLMC","Peri","SMC","Pax6"))) and \
    any(k in str(l) for k in ("IT","ET","CT","NP","L6b","Sst","Pvalb","Vip","Lamp5","Sncg","Chandelier"))
is_neuronal_pfc = lambda l: str(l).startswith("InN") or str(l).startswith("ExN")

off_rows = np.where(off)[0]
appl_rows, c4_rates = [], []
for mf, rule in [("Adult_Human_MTG", is_neuronal_mtg), ("Adult_Human_PrefrontalCortex", is_neuronal_pfc)]:
    # 适用性: 神经元参考池的目标类型准确率（≥80% 才可用）
    neu_ref = np.concatenate([refs.get("Excitatory", []), refs.get("Inhibitory", [])]).astype(int)
    lab_ref, nmatch, nfeat = predict_rows(neu_ref, mf)
    acc_neu = np.mean([rule(l) for l in lab_ref])
    log(f"{mf}: 特征匹配 {nmatch}/{nfeat}, 神经元参考准确率 {acc_neu*100:.1f}%")
    applicable = acc_neu >= 0.80
    # 官方子集注释
    lab_off, _, _ = predict_rows(off_rows, mf)
    pd.Series(lab_off).value_counts().to_csv(os.path.join(OUT_DIR, f"celltypist_{mf}_official.csv"))
    nonneu = 1 - np.mean([rule(l) for l in lab_off])
    log(f"{mf}: 官方 neuron 子集非神经元率 = {nonneu*100:.1f}%  适用性={'PASS' if applicable else 'FAIL'}")
    appl_rows.append({"model": mf, "features_matched": nmatch, "features_total": nfeat,
                      "neuron_ref_accuracy_pct": round(acc_neu*100,1),
                      "applicable": applicable,
                      "official_nonneuronal_pct": round(nonneu*100,1)})
    if applicable:
        c4_rates.append(nonneu)
pd.DataFrame(appl_rows).to_csv(os.path.join(OUT_DIR, "applicability_table.csv"), index=False)

c4 = float(np.mean(c4_rates)) if c4_rates else np.nan
checks = {"C1": round(c1*100,1), "C2": round(c2*100,1), "C3": round(c3*100,1),
          "C4": round(c4*100,1) if not np.isnan(c4) else None}
acs = float(np.mean([c1*100, c2*100, c3*100] + ([c4*100] if not np.isnan(c4) else [])))
final = {
    "dataset": "GSE330130 human ALS lumbar spinal cord snRNA-seq",
    "official_label": "Neuron (n=%d)" % off.sum(),
    "panel": f"{NGENE} human ortholog marker genes (mouse 32-gene panel orthologs, excl. SLC17A6)",
    "checks": checks,
    "C2_ungated_argmax_pct": round(c2_argmax*100,1),
    "C3_protocol": "60k subsample, top-2000 HVG, PCA50, KMeans k=20, cluster-majority mapping (32-gene panel)",
    "C4_applicable_models": [r["model"] for r in appl_rows if r["applicable"]],
    "C4_excluded_models": [r["model"] for r in appl_rows if not r["applicable"]] + ["Developing_Mouse_Brain (cross-species: 4/7416 features)"],
    "ACS_pct": round(acs,1),
    "protocol": "C1 argmax(threshold>0); C2 scanpy score_genes + margin tau=1.5; C3 KMeans(60k, k=20); C4 applicability-gated CellTypist (>=80%)",
}
import json
with open(os.path.join(OUT_DIR, "final_audit.json"), "w") as f:
    json.dump(final, f, indent=2, ensure_ascii=False)
log("最终审计: " + json.dumps(final, ensure_ascii=False))
log("V14 DONE")
