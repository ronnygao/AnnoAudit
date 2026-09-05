#!/usr/bin/env python3
# ============================================================
# classify_all_labels.py — 两图谱全官方标签 × 判别身份 汇总（出版用表）
# 对每个官方标签的每个细胞，用独立基因判别"最可能真实身份"：
#   Neuron_like / Microglia_like / Oligo_like / Astro_like / (Endo/OPC 由 32 panel argmax 兜底)
# 输出：每官方标签的判别身份构成 + 判定（clean / contaminated / ambiguous）
# ============================================================
import os, glob, gzip, json, time, warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import anndata as ad
import scipy.io, scipy.sparse as sp

T0 = time.time()
def log(m): print(f"[{time.time()-T0:7.1f}s] {m}", flush=True)
OUT = "/Users/apple/Desktop/颅脑损伤专病数据库/P2_离子通道组景观论文/AnnoAudit_survey/classification"
os.makedirs(OUT, exist_ok=True)

# ================= GSE330130（全部 9 类） =================
log("== GSE330130 ==")
H5 = "/Users/apple/Desktop/颅脑损伤专病数据库/P2_离子通道组景观论文/AnnoAudit_survey/data/GSE330130_LSC_SNRNASEQ_COUNTS.h5ad"
a = ad.read_h5ad(H5)
a.var["symbol"] = a.var["common_name"].astype(str)
a.var_names = a.var["symbol"].values
a = a[:, ~pd.Index(a.var_names).duplicated(keep="first")].copy()
X = a.X.tocsr().astype(np.float32)
umi = np.asarray(X.sum(axis=1)).ravel()
Xn = X.multiply(1e6 / np.maximum(umi, 1)[:, None]).tocsr()
Xn.data = np.log1p(Xn.data)
VN = a.var_names
def sm(genes):
    cols = [VN.get_loc(g) for g in genes if g in VN]
    return np.asarray(Xn[:, cols].mean(axis=1)).ravel() if cols else np.zeros(a.n_obs)
s_neur = sm(["SNAP25","STMN2","SYT1"]); s_olig = sm(["MBP","PLP1","MOG","MAG"])
s_astr = sm(["GFAP","AQP4","SLC1A3","ALDH1L1"]); s_mg = sm(["TREM2","CD68","CX3CR1","HEXB","AIF1"])
s_endo = sm(["PECAM1","CLDN5","FLT1"]); s_opc = sm(["PDGFRA","OLIG1","CSPG4"])
# argmax 型（8 类 + Unknown）
score_map = {"Excitatory": sm(["SLC17A7","RBFOX3","NEUROD6","TBR1","CAMK2A"]),
             "Inhibitory": sm(["GAD1","GAD2","SST","PVALB","NPY"]),
             "Astrocyte": s_astr, "Microglia": s_mg, "Oligodendrocyte": s_olig,
             "OPC": s_opc, "Endothelial": s_endo, "Pericyte": sm(["PDGFRB","RGS5","VTN"])}
sdf = pd.DataFrame(score_map)
argmax_t = sdf.idxmax(axis=1).values
argmax_v = sdf.max(axis=1).values
argmax_t = np.where(argmax_v > 0, argmax_t, "Unknown")
# 独立判别（4 主身份）
ident = np.select(
    [s_mg >= np.maximum(s_neur, s_olig) * 1.2, s_olig >= np.maximum(s_neur, s_mg) * 1.2,
     s_astr >= np.maximum(s_neur, s_mg) * 1.2],
    ["Microglia", "Oligodendrocyte", "Astrocyte"], default="Neuron")
# 神经元身份再用 argmax 细拆 Ex/Inh
ident = np.where((ident == "Neuron") & (argmax_t == "Inhibitory"), "Inhibitory", ident)
ident = np.where((ident == "Neuron") & (argmax_t == "Excitatory"), "Excitatory", ident)
df = pd.DataFrame({"official": a.obs["Cell_Type"].astype(str).values, "ident": ident,
                   "s_neur": s_neur.round(3), "s_olig": s_olig.round(3), "s_astr": s_astr.round(3), "s_mg": s_mg.round(3)})
SELF330 = {"Oligodendrocytes": "Oligodendrocyte", "Astrocyte": "Astrocyte", "Microglia": "Microglia",
           "OPCs": "OPC", "Endothelial": "Endothelial",
           "Inhibitory Neuron": "Inhibitory", "Excitatory  Neuron": "Excitatory",
           "Motor Neurons": "Neuron", "Meninges": None}
rows = []
for lab in df["official"].unique():
    sub = df[df["official"] == lab]
    vc = sub["ident"].value_counts()
    self_t = SELF330[lab]
    match = vc.get(self_t, 0) / len(sub) if self_t else np.nan
    rows.append({"atlas": "GSE330130", "official_label": lab, "n": len(sub),
                 "matched_identity_pct": round(match * 100, 1) if not np.isnan(match) else None,
                 "top_ident": vc.index[0], "top_pct": round(vc.iloc[0] / len(sub) * 100, 1),
                 "compo": "; ".join(f"{k}={round(v/len(sub)*100,1)}%" for k, v in vc.items())})
df.to_csv(os.path.join(OUT, "gse330130_all_labels_identity.csv"), index=False)

# ================= CEREBRI（全部 11 类） =================
log("\n== CEREBRI ==")
RAW_DIR = "/Users/apple/Workbuddy/2026-07-22-16-57-18/cerebri_data/raw"
ANNOT_CSV = "/Users/apple/Desktop/颅脑损伤专病数据库/CEREBRI_KCNC3生信分析/v3分析结果/cell_annotations.csv"
ann = pd.read_csv(ANNOT_CSV)
ann["bc_plain"] = ann["barcode"].str.replace(r"_GSM\d+$", "", regex=True)
feat_path = sorted(glob.glob(os.path.join(RAW_DIR, "*_features.tsv.gz")))[0]
with gzip.open(feat_path, "rt") as f:
    feats = [l.strip().split("\t") for l in f if not l.startswith("#")]
gene_names = [x[1] if len(x) > 1 else x[0] for x in feats]
gene_idx = {g: i for i, g in enumerate(gene_names)}
MK = {"Excitatory": ["Slc17a7","Rbfox3","Neurod6","Tbr1","Camk2a"], "Inhibitory": ["Gad1","Gad2","Sst","Pvalb","Npy"],
      "Astrocyte": ["Gfap","Aqp4","Slc1a3","Aldh1l1"], "Microglia": ["C1qb","Tyrobp","Cx3cr1","Hexb","Aif1"],
      "Oligodendrocyte": ["Mbp","Plp1","Mog","Mag"], "OPC": ["Pdgfra","Olig1","Cspg4"],
      "Endothelial": ["Pecam1","Cldn5","Flt1"], "Pericyte": ["Pdgfrb","Rgs5","Vtn"]}
want = sorted(set(g for gs in MK.values() for g in gs) | {"Trem2","Cd68","Lyz2","Gja1","Slc1a2","S100b","Snap25","Stmn2","Syt1"})
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
Xc = sp.vstack(Xp).tocsr(); bcm = pd.concat(bcp, ignore_index=True)
umi = np.asarray(Xc.sum(axis=1)).ravel()
Xcn = Xc.multiply(1e6 / np.maximum(umi, 1)[:, None]).tocsr()
Xcn.data = np.log1p(Xcn.data)
inv = {i: g for g, i in widx.items()}
colg = [inv[i] for i in sorted(set(widx.values()))]; gpos = {g: j for j, g in enumerate(colg)}
def smc(genes):
    cols = [gpos[g] for g in genes if g in gpos]
    return np.asarray(Xcn[:, cols].mean(axis=1)).ravel() if cols else np.zeros(Xcn.shape[0])
scores_c = {t: smc(gs) for t, gs in MK.items()}
s_neur = smc(["Snap25","Stmn2","Syt1"]); s_olig = scores_c["Oligodendrocyte"]
s_astr = scores_c["Astrocyte"]; s_mg = scores_c["Microglia"]
sdf = pd.DataFrame(scores_c); amt = sdf.idxmax(axis=1).values; amv = sdf.max(axis=1).values
amt = np.where(amv > 0, amt, "Unknown")
ident = np.select(
    [s_mg >= np.maximum(s_neur, s_olig) * 1.2, s_olig >= np.maximum(s_neur, s_mg) * 1.2,
     s_astr >= np.maximum(s_neur, s_mg) * 1.2],
    ["Microglia", "Oligodendrocyte", "Astrocyte"], default="Neuron")
ident = np.where((ident == "Neuron") & (amt == "Inhibitory"), "Inhibitory", ident)
ident = np.where((ident == "Neuron") & (amt == "Excitatory"), "Excitatory", ident)
mrg = bcm.reset_index().rename(columns={"index": "xrow"})
m2 = ann.merge(mrg, left_on=["gsm", "bc_plain"], right_on=["gsm", "barcode"], how="left")
xr = m2["xrow"].values
dc = pd.DataFrame({"official": ann["cell_type"].astype(str).values, "ident": ident[xr]})
SELFC = {"Microglia": "Microglia", "Glutamatergic_neuron": "Excitatory", "GABAergic_neuron": "Inhibitory",
         "Astrocyte": "Astrocyte", "Oligodendrocyte": "Oligodendrocyte", "OPC": "OPC",
         "Endothelial": "Endothelial", "Pericyte": "Pericyte", "VSM": "Pericyte",
         "Fibroblast": None, "Erythrocyte": None, "Unknown": None}
for lab in dc["official"].unique():
    sub = dc[dc["official"] == lab]
    vc = sub["ident"].value_counts()
    self_t = SELFC.get(lab)
    match = vc.get(self_t, 0) / len(sub) if self_t else np.nan
    rows.append({"atlas": "CEREBRI", "official_label": lab, "n": len(sub),
                 "matched_identity_pct": round(match * 100, 1) if not np.isnan(match) else None,
                 "top_ident": vc.index[0], "top_pct": round(vc.iloc[0] / len(sub) * 100, 1),
                 "compo": "; ".join(f"{k}={round(v/len(sub)*100,1)}%" for k, v in vc.items())})
dc.to_csv(os.path.join(OUT, "cerebri_all_labels_identity.csv"), index=False)

res = pd.DataFrame(rows).sort_values(["atlas", "n"], ascending=[True, False])
res.to_csv(os.path.join(OUT, "all_labels_classification_summary.csv"), index=False)
print("\n===== 汇总（判别身份匹配率 = 官方标签与独立基因判别一致的细胞比例）=====")
print(res.to_string(index=False))
log("DONE")
