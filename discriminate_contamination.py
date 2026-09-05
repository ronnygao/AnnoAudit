#!/usr/bin/env python3
# ============================================================
# discriminate_contamination.py — 真错标 vs ambient 双阳判别
# 对官方神经元/实质标签中"被判非本型"的细胞，用独立基因得分分类：
#   真错标 = 异型基因高 + 本型/神经元基因低（如 CEREBRI 微胶质流入）
#   ambient = 本型基因仍高（双阳，如 GSE330130 少突流入）
# 输出每官方标签的"判别矩阵"与校正后污染率
# ============================================================
import os, time, json, warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import anndata as ad

T0 = time.time()
def log(m): print(f"[{time.time()-T0:7.1f}s] {m}", flush=True)
OUT = "/Users/apple/Desktop/颅脑损伤专病数据库/P2_离子通道组景观论文/AnnoAudit_survey/discrimination"
os.makedirs(OUT, exist_ok=True)

def score_mean(Xn, var_names, genes):
    cols = [var_names.get_loc(g) for g in genes if g in var_names]
    if not cols:
        return np.zeros(Xn.shape[0])
    return np.asarray(Xn[:, cols].mean(axis=1)).ravel()

# ---------------- GSE330130 ----------------
log("== GSE330130 ==")
H5 = os.path.join(os.path.dirname(OUT), "data/GSE330130_LSC_SNRNASEQ_COUNTS.h5ad")
a = ad.read_h5ad(H5)
a.var["symbol"] = a.var["common_name"].astype(str)
a.var_names = a.var["symbol"].values
a = a[:, ~pd.Index(a.var_names).duplicated(keep="first")].copy()
X = a.X.tocsr().astype(np.float32)
umi = np.asarray(X.sum(axis=1)).ravel()
Xn = X.multiply(1e4 / np.maximum(umi, 1)[:, None]).tocsr()   # CP10K（与 CellTypist 一致）
Xn.data = np.log1p(Xn.data)
VN = a.var_names
NEU_LAB = ["Inhibitory Neuron", "Excitatory  Neuron", "Motor Neurons"]
off = a.obs["Cell_Type"].astype(str).isin(NEU_LAB).values
log(f"官方神经元 n = {off.sum()}")

# 得分
s_neur = score_mean(Xn, VN, ["SNAP25", "STMN2", "SYT1"])
s_olig = score_mean(Xn, VN, ["MBP", "PLP1", "MOG", "MAG"])
s_astr = score_mean(Xn, VN, ["GFAP", "AQP4", "SLC1A3", "ALDH1L1"])
s_mg   = score_mean(Xn, VN, ["C1QB", "TYROBP", "CX3CR1", "HEXB", "AIF1"])

# 参考值（真细胞类型的得分中位数，用官方纯净标签）
ref = {}
for lab, m in [("Oligo", s_olig), ("Astro", s_astr), ("Mg", s_mg), ("Neuron", s_neur)]:
    for o, name in [("Oligodendrocytes", "Oligo"), ("Astrocyte", "Astro"), ("Microglia", "Mg")]:
        pass
# 简化：用官方 Oligodendrocytes/Astrocyte/Microglia 标签（纯度高）做参考
olig_pure = a.obs["Cell_Type"].astype(str) == "Oligodendrocytes"
astr_pure = a.obs["Cell_Type"].astype(str) == "Astrocyte"
mg_pure   = a.obs["Cell_Type"].astype(str) == "Microglia"
ref_olig_neur = np.median(s_neur[olig_pure.values]); ref_olig_olig = np.median(s_olig[olig_pure.values])
ref_astr_neur = np.median(s_neur[astr_pure.values])
ref_mg_neur   = np.median(s_neur[mg_pure.values])
log(f"参考: 真少突 SNAP25组中位 {ref_olig_neur:.2f}, MBP组中位 {ref_olig_olig:.2f}; 真星形神经基因中位 {ref_astr_neur:.2f}")
# 判别阈值：神经元身份 = s_neur > (ref_olig_neur + ref_astr_neur)/2 附近 → 用 1.0 (log1p CP10K ~ SNAP25 检测)

df = pd.DataFrame({"official": a.obs["Cell_Type"].astype(str).values,
                   "neur": s_neur, "olig": s_olig, "astr": s_astr, "mg": s_mg})
df["identity"] = np.select(
    [df["mg"] > df["neur"] * 1.5, df["olig"] > df["neur"] * 1.5, df["astr"] > df["neur"] * 1.5],
    ["Microglia_like", "Oligo_like", "Astro_like"], default="Neuron_like")
for lab in NEU_LAB:
    sub = df[df["official"] == lab]
    vc = sub["identity"].value_counts()
    log(f"官方 {lab} (n={len(sub)}): " + "; ".join(f"{k}={v} ({v/len(sub)*100:.1f}%)" for k, v in vc.items()))
    # 判别: Neuron_like = 本质神经元(含 ambient 双阳)
    true_neu = (sub["identity"] == "Neuron_like").mean()
    log(f"  → 校正后真神经元占比 {true_neu*100:.1f}%, 真非神经元(胶质样) {(1-true_neu)*100:.1f}%")
out330 = df[df["official"].isin(NEU_LAB)]
out330.to_csv(os.path.join(OUT, "gse330130_neuron_identity.csv"), index=False)

# ---------------- CEREBRI ----------------
log("\n== CEREBRI ==")
RAW_DIR = "/Users/apple/Workbuddy/2026-07-22-16-57-18/cerebri_data/raw"
ANNOT_CSV = "/Users/apple/Desktop/颅脑损伤专病数据库/CEREBRI_KCNC3生信分析/v3分析结果/cell_annotations.csv"
import glob, gzip, scipy.io, scipy.sparse as sp
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
VERIFY_G = ["Trem2","Cd68","Lyz2","Gja1","Slc1a2","Snap25","Stmn2","Syt1","Mbp","Plp1"]
ann = pd.read_csv(ANNOT_CSV)
ann["bc_plain"] = ann["barcode"].str.replace(r"_GSM\d+$", "", regex=True)
feat_path = sorted(glob.glob(os.path.join(RAW_DIR, "*_features.tsv.gz")))[0]
with gzip.open(feat_path, "rt") as f:
    feats = [l.strip().split("\t") for l in f if not l.startswith("#")]
gene_names = [x[1] if len(x) > 1 else x[0] for x in feats]
gene_idx = {g: i for i, g in enumerate(gene_names)}
want = sorted(set(g for gs in MARKERS.values() for g in gs) | set(VERIFY_G))
widx = {g: gene_idx[g] for g in want if g in gene_idx}
Xp, bcp = [], []
for gsm in sorted(ann["gsm"].unique()):
    mtx_f = glob.glob(os.path.join(RAW_DIR, f"{gsm}_*_matrix.mtx.gz"))
    bar_f = glob.glob(os.path.join(RAW_DIR, f"{gsm}_*_barcodes.tsv.gz"))
    if not mtx_f: continue
    with gzip.open(mtx_f[0], "rb") as f: mat = scipy.io.mmread(f)
    mat = mat.tocsr()
    sub = mat[sorted(set(widx.values())), :].T.tocsr().astype(np.float32)
    del mat
    with gzip.open(bar_f[0], "rt") as f: barcodes = [l.strip() for l in f]
    Xp.append(sub); bcp.append(pd.DataFrame({"gsm": gsm, "barcode": barcodes}))
Xc = sp.vstack(Xp).tocsr()
bcm = pd.concat(bcp, ignore_index=True)
umi = np.asarray(Xc.sum(axis=1)).ravel()
Xcn = Xc.multiply(1e4 / np.maximum(umi, 1)[:, None]).tocsr()
Xcn.data = np.log1p(Xcn.data)
inv = {i: g for g, i in widx.items()}
colg = [inv[i] for i in sorted(set(widx.values()))]
gpos = {g: j for j, g in enumerate(colg)}
def s_mean(genes):
    cols = [gpos[g] for g in genes if g in gpos]
    return np.asarray(Xcn[:, cols].mean(axis=1)).ravel() if cols else np.zeros(Xcn.shape[0])
s_neur = s_mean(["Snap25","Stmn2","Syt1"]); s_olig = s_mean(["Mbp","Plp1","Mog","Mag"])
s_astr = s_mean(["Gfap","Aqp4","Slc1a3","Aldh1l1"]); s_mg = s_mean(["C1qb","Tyrobp","Cx3cr1","Hexb","Aif1"])
# 对齐注释
mrg = ann.merge(bcm, left_on=["gsm","bc_plain"], right_on=["gsm","barcode"], how="left")
rows = mrg.index  # ann 行序 = Xcn 行序（bcm 拼接序）？
# bcm 拼接序与 Xcn 行一致；mrg 保留 ann 行序但丢 X 行号 → 重建
mrg2 = bcm.reset_index().rename(columns={"index": "xrow"})
m2 = ann.merge(mrg2, left_on=["gsm","bc_plain"], right_on=["gsm","barcode"], how="left")
xr = m2["xrow"].values
NEU_LAB_C = ["Glutamatergic_neuron", "GABAergic_neuron"]
dc = pd.DataFrame({"official": ann["cell_type"].astype(str).values,
                   "neur": s_neur[xr], "olig": s_olig[xr], "astr": s_astr[xr], "mg": s_mg[xr]})
dc["identity"] = np.select(
    [dc["mg"] > dc["neur"] * 1.5, dc["olig"] > dc["neur"] * 1.5, dc["astr"] > dc["neur"] * 1.5],
    ["Microglia_like", "Oligo_like", "Astro_like"], default="Neuron_like")
for lab in NEU_LAB_C + ["Astrocyte", "OPC"]:
    sub = dc[dc["official"] == lab]
    vc = sub["identity"].value_counts()
    log(f"官方 {lab} (n={len(sub)}): " + "; ".join(f"{k}={v} ({v/len(sub)*100:.1f}%)" for k, v in vc.items()))
dc.to_csv(os.path.join(OUT, "cerebri_labels_identity.csv"), index=False)
log("\nDONE — 输出 " + OUT)
