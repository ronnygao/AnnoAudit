#!/usr/bin/env python3
"""时间结构验证：CEREBRI 官方标签的污染是否随损伤条件（时间）变化？
对每官方标签 × condition，计算 marker 组成（微胶质占比等）→ 检验"状态依赖错标"假说
"""
import os, glob, gzip, time, warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import scipy.io, scipy.sparse as sp

T0 = time.time()
def log(m): print(f"[{time.time()-T0:7.1f}s] {m}", flush=True)
OUT = "/Users/apple/Desktop/颅脑损伤专病数据库/P2_离子通道组景观论文/AnnoAudit_survey/time_structure"
os.makedirs(OUT, exist_ok=True)

RAW_DIR = "/Users/apple/Workbuddy/2026-07-22-16-57-18/cerebri_data/raw"
ANNOT_CSV = "/Users/apple/Desktop/颅脑损伤专病数据库/CEREBRI_KCNC3生信分析/v3分析结果/cell_annotations.csv"
MK = {"Excitatory": ["Slc17a7","Rbfox3","Neurod6","Tbr1","Camk2a"],
      "Inhibitory": ["Gad1","Gad2","Sst","Pvalb","Npy"],
      "Astrocyte": ["Gfap","Aqp4","Slc1a3","Aldh1l1"],
      "Microglia": ["C1qb","Tyrobp","Cx3cr1","Hexb","Aif1"],
      "Oligodendrocyte": ["Mbp","Plp1","Mog","Mag"],
      "OPC": ["Pdgfra","Olig1","Cspg4"],
      "Endothelial": ["Pecam1","Cldn5","Flt1"],
      "Pericyte": ["Pdgfrb","Rgs5","Vtn"]}
ann = pd.read_csv(ANNOT_CSV)
ann["bc_plain"] = ann["barcode"].str.replace(r"_GSM\d+$", "", regex=True)
# condition 列时间分组（Naive=对照; 24h=急性; 6mo=慢性; 7d=?）
cond_map = {}
for c in ann["condition"].unique():
    cs = str(c)
    if "Naive" in cs or "Sham" in cs or "naive" in cs: cond_map[c] = "Naive/control"
    elif "6mo" in cs or "6m" in cs: cond_map[c] = "6mo_chronic"
    elif "7d" in cs or "7_day" in cs: cond_map[c] = "7d"
    elif "24h" in cs: cond_map[c] = "24h_acute"
    else: cond_map[c] = cs
ann["time_group"] = ann["condition"].map(cond_map)
log("时间组分布: " + ann["time_group"].value_counts().to_string().replace("\n", " | "))

feat_path = sorted(glob.glob(os.path.join(RAW_DIR, "*_features.tsv.gz")))[0]
with gzip.open(feat_path, "rt") as f:
    feats = [l.strip().split("\t") for l in f if not l.startswith("#")]
gene_names = [x[1] if len(x) > 1 else x[0] for x in feats]
gene_idx = {g: i for i, g in enumerate(gene_names)}
want = sorted(set(g for gs in MK.values() for g in gs))
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
colg = [inv[i] for i in sorted(set(widx.values()))]
gpos = {g: j for j, g in enumerate(colg)}
scs = {}
for ct, gs in MK.items():
    cols = [gpos[g] for g in gs if g in gpos]
    scs[ct] = np.asarray(Xcn[:, cols].mean(axis=1)).ravel()
sdf = pd.DataFrame(scs); amt = sdf.idxmax(axis=1).values; amv = sdf.max(axis=1).values
amt = np.where(amv > 0, amt, "Unknown")
mrg = bcm.reset_index().rename(columns={"index": "xrow"})
m2 = ann.merge(mrg, left_on=["gsm", "bc_plain"], right_on=["gsm", "barcode"], how="left")
xr = m2["xrow"].values
mt_x = amt[xr]  # 与 ann 行序一致
ann["marker_type"] = mt_x

# 关键标签 × 时间组：微胶质占比 + 自身 marker 占比
FOCUS = ["Glutamatergic_neuron", "GABAergic_neuron", "Astrocyte", "OPC", "Microglia"]
rows = []
for lab in FOCUS:
    sub = ann[ann["cell_type"].astype(str) == lab]
    for tg, grp in sub.groupby("time_group"):
        n = len(grp)
        if n < 50:
            continue
        mg = (grp["marker_type"] == "Microglia").mean() * 100
        self_t = {"Glutamatergic_neuron": "Excitatory", "GABAergic_neuron": "Inhibitory",
                  "Astrocyte": "Astrocyte", "OPC": "OPC", "Microglia": "Microglia"}[lab]
        self_frac = (grp["marker_type"] == self_t).mean() * 100
        rows.append({"official_label": lab, "time": tg, "n": n,
                     "microglia_pct": round(mg, 1), "self_marker_pct": round(self_frac, 1)})
res = pd.DataFrame(rows)
order_t = ["Naive/control", "24h_acute", "7d", "6mo_chronic"]
res["t_ord"] = res["time"].map({t: i for i, t in enumerate(order_t)})
res = res.sort_values(["official_label", "t_ord"])
res.to_csv(os.path.join(OUT, "label_by_time.csv"), index=False)
print(res.to_string(index=False))
log("DONE")
