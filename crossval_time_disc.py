#!/usr/bin/env python3
"""交叉验证 3：判别口径的时间分层——急性期污染开关在 D1（τ=1.5）下是否成立
若 marker 口径（对照干净/急性崩溃）与判别口径一致 → 状态依赖声明双重验证
"""
import os, glob, gzip, warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import scipy.io, scipy.sparse as sp

PANEL = {"Excitatory": ["Slc17a7","Rbfox3","Neurod6","Tbr1","Camk2a"],
         "Inhibitory": ["Gad1","Gad2","Sst","Pvalb","Npy"],
         "Astrocyte": ["Gfap","Aqp4","Slc1a3","Aldh1l1"],
         "Microglia": ["C1qb","Tyrobp","Cx3cr1","Hexb","Aif1"],
         "Oligodendrocyte": ["Mbp","Plp1","Mog","Mag"],
         "OPC": ["Pdgfra","Olig1","Cspg4"],
         "Endothelial": ["Pecam1","Cldn5","Flt1"],
         "Pericyte": ["Pdgfrb","Rgs5","Vtn"]}
MK8 = list(PANEL.keys())
INDEP = {"Excitatory": ["Snap25","Stmn2","Syt1"], "Inhibitory": ["Snap25","Stmn2","Syt1"],
         "Microglia": ["Trem2","Cd68","Lyz2"], "Astrocyte": ["Gja1","Slc1a2","S100b"],
         "Oligodendrocyte": ["Mbp","Plp1","Mog","Mag"], "OPC": ["Pdgfra","Olig1","Cspg4"],
         "Endothelial": ["Pecam1","Cldn5","Flt1"], "Pericyte": ["Pdgfrb","Rgs5","Vtn"]}
TAU = 1.5
RAW_DIR = "/Users/apple/Workbuddy/2026-07-22-16-57-18/cerebri_data/raw"
ANNOT_CSV = "/Users/apple/Desktop/颅脑损伤专病数据库/CEREBRI_KCNC3生信分析/v3分析结果/cell_annotations.csv"

ann = pd.read_csv(ANNOT_CSV)
ann["bc_plain"] = ann["barcode"].str.replace(r"_GSM\d+$", "", regex=True)
def tg(c):
    cs = str(c)
    if "Naive" in cs or "Sham" in cs or "naive" in cs: return "control"
    if "6mo" in cs or "6m" in cs: return "6mo"
    if "7d" in cs or "7_day" in cs: return "7d"
    if "24h" in cs: return "24h"
    return cs
ann["tgroup"] = ann["condition"].map(tg)

feat_path = sorted(glob.glob(os.path.join(RAW_DIR, "*_features.tsv.gz")))[0]
with gzip.open(feat_path, "rt") as f:
    feats = [l.strip().split("\t") for l in f if not l.startswith("#")]
gene_names = [x[1] if len(x) > 1 else x[0] for x in feats]
gene_idx = {g: i for i, g in enumerate(gene_names)}
want = sorted(set(g for gs in PANEL.values() for g in gs) | set(g for gs in INDEP.values() for g in gs))
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
def sm(genes):
    cols = [gpos[g] for g in genes if g in gpos]
    return np.asarray(Xcn[:, cols].mean(axis=1)).ravel() if cols else np.zeros(Xcn.shape[0])
scs = {ct: sm(PANEL[ct]) for ct in MK8}
sdf = pd.DataFrame(scs); amt = sdf.idxmax(axis=1).values; amv = sdf.max(axis=1).values
amt = np.where(amv > 0, amt, "Unknown")
inds = {ct: sm(INDEP[ct]) for ct in MK8}
mrg = bcm.reset_index().rename(columns={"index": "xrow"})
m2 = ann.merge(mrg, left_on=["gsm", "bc_plain"], right_on=["gsm", "barcode"], how="left")
xr = m2["xrow"].values
SELF = {"Glutamatergic_neuron": "Excitatory", "GABAergic_neuron": "Inhibitory",
        "Astrocyte": "Astrocyte", "OPC": "OPC"}

rows = []
for lab, self_t in SELF.items():
    lab_ann = ann[ann["cell_type"].astype(str) == lab]
    for tgname, grp in lab_ann.groupby("tgroup"):
        idx = grp.index
        x_idx = xr[idx]
        sub_amt = amt[x_idx]
        n = len(idx)
        if n < 100: continue
        non_self = sub_amt != self_t
        neu_other = ((sub_amt == "Inhibitory") & (self_t == "Excitatory")) | ((sub_amt == "Excitatory") & (self_t == "Inhibitory"))
        cand = non_self & ~neu_other
        s_self = inds[self_t][x_idx]
        conf = np.zeros(n, dtype=bool)
        for o in MK8:
            if o == self_t: continue
            m_o = (sub_amt == o) & cand
            if m_o.any():
                conf |= m_o & (inds[o][x_idx] > TAU * s_self)
        rows.append({"label": lab, "time": tgname, "n": n,
                     "marker_self_pct": round((sub_amt == self_t).mean() * 100, 1),
                     "disc_contamination_pct": round(conf.mean() * 100, 1),
                     "microglia_argmax_pct": round((sub_amt == "Microglia").mean() * 100, 1)})
res = pd.DataFrame(rows)
order = ["control", "24h", "7d", "6mo"]
res["ord"] = res["time"].map({t: i for i, t in enumerate(order)})
res = res.sort_values(["label", "ord"])
res.to_csv("/Users/apple/Desktop/颅脑损伤专病数据库/P2_离子通道组景观论文/AnnoAudit_survey/crossval_time_discrimination.csv", index=False)
print(res.to_string(index=False))
