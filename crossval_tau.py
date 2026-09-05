#!/usr/bin/env python3
"""交叉验证 1：判别污染率对 τ 的敏感性（τ ∈ {1.0,1.25,1.5,2.0,3.0}）
验证重构手稿核心声明"污染率"不依赖人为选定的 τ=1.5。
CEREBRI 关键标签 + GSE330130 神经元标签。
"""
import os, glob, gzip, warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import scipy.io, scipy.sparse as sp
import anndata as ad

PANEL = {"Excitatory": ["Slc17a7","Rbfox3","Neurod6","Tbr1","Camk2a"],
         "Inhibitory": ["Gad1","Gad2","Sst","Pvalb","Npy"],
         "Astrocyte": ["Gfap","Aqp4","Slc1a3","Aldh1l1"],
         "Microglia": ["C1qb","Tyrobp","Cx3cr1","Hexb","Aif1"],
         "Oligodendrocyte": ["Mbp","Plp1","Mog","Mag"],
         "OPC": ["Pdgfra","Olig1","Cspg4"],
         "Endothelial": ["Pecam1","Cldn5","Flt1"],
         "Pericyte": ["Pdgfrb","Rgs5","Vtn"]}
MK8 = list(PANEL.keys())
NEU_IND = ["Snap25","Stmn2","Syt1"]
INDEP = {"Excitatory": NEU_IND, "Inhibitory": NEU_IND,
         "Microglia": ["Trem2","Cd68","Lyz2"], "Astrocyte": ["Gja1","Slc1a2","S100b"],
         "Oligodendrocyte": ["Mbp","Plp1","Mog","Mag"], "OPC": ["Pdgfra","Olig1","Cspg4"],
         "Endothelial": ["Pecam1","Cldn5","Flt1"], "Pericyte": ["Pdgfrb","Rgs5","Vtn"]}
TAUS = [1.0, 1.25, 1.5, 2.0, 3.0]

def run(expr, var_names, obs_lab, self_map, labels_of_interest, dataset):
    gpos = {g: j for j, g in enumerate(var_names)}
    def sm(genes):
        cols = [gpos[g] for g in genes if g in gpos]
        return np.asarray(expr[:, cols].mean(axis=1)).ravel() if cols else np.zeros(expr.shape[0])
    scs = {ct: sm(PANEL[ct]) for ct in MK8}
    sdf = pd.DataFrame(scs); amt = sdf.idxmax(axis=1).values; amv = sdf.max(axis=1).values
    amt = np.where(amv > 0, amt, "Unknown")
    inds = {ct: sm(INDEP[ct]) for ct in MK8}
    rows = []
    for lab in labels_of_interest:
        self_t = self_map[lab]
        mask = obs_lab == lab
        n = int(mask.sum())
        sub_amt = amt[mask]
        non_self = sub_amt != self_t
        neu_other = ((sub_amt == "Inhibitory") & (self_t == "Excitatory")) | ((sub_amt == "Excitatory") & (self_t == "Inhibitory"))
        cand = non_self & ~neu_other
        s_self = inds[self_t][mask]
        row = {"dataset": dataset, "label": lab, "n": n}
        for tau in TAUS:
            conf = np.zeros(n, dtype=bool)
            for o in MK8:
                if o == self_t:
                    continue
                m_o = (sub_amt == o) & cand
                if m_o.any():
                    conf |= m_o & (inds[o][mask] > tau * s_self)
            row[f"tau{tau}"] = round(conf.mean() * 100, 1)
        rows.append(row)
    return rows

# ---- CEREBRI ----
RAW_DIR = "/Users/apple/Workbuddy/2026-07-22-16-57-18/cerebri_data/raw"
ANNOT_CSV = "/Users/apple/Desktop/颅脑损伤专病数据库/CEREBRI_KCNC3生信分析/v3分析结果/cell_annotations.csv"
ann = pd.read_csv(ANNOT_CSV)
ann["bc_plain"] = ann["barcode"].str.replace(r"_GSM\d+$", "", regex=True)
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
mrg = bcm.reset_index().rename(columns={"index": "xrow"})
m2 = ann.merge(mrg, left_on=["gsm", "bc_plain"], right_on=["gsm", "barcode"], how="left")
xr = m2["xrow"].values
obs_x = np.full(Xcn.shape[0], "Unmatched", dtype=object)
obs_x[xr] = ann["cell_type"].astype(str).values
SELF_C = {"Glutamatergic_neuron": "Excitatory", "GABAergic_neuron": "Inhibitory", "Astrocyte": "Astrocyte",
          "OPC": "OPC", "Microglia": "Microglia"}
r1 = run(Xcn, colg, obs_x, SELF_C, ["Glutamatergic_neuron", "GABAergic_neuron", "Astrocyte", "OPC"], "CEREBRI")
print("CEREBRI τ 敏感性:")
print(pd.DataFrame(r1).to_string(index=False))

# ---- GSE330130 ----
H5 = "/Users/apple/Desktop/颅脑损伤专病数据库/P2_离子通道组景观论文/AnnoAudit_survey/data/GSE330130_LSC_SNRNASEQ_COUNTS.h5ad"
a = ad.read_h5ad(H5)
a.var["symbol"] = a.var["common_name"].astype(str)
a.var_names = a.var["symbol"].values
a = a[:, ~pd.Index(a.var_names).duplicated(keep="first")].copy()
X = a.X.tocsr().astype(np.float32)
umi = np.asarray(X.sum(axis=1)).ravel()
Xn = X.multiply(1e6 / np.maximum(umi, 1)[:, None]).tocsr()
Xn.data = np.log1p(Xn.data)
globals()["PANEL"] = {k: [g.upper() for g in v] for k, v in PANEL.items()}
globals()["INDEP"] = {"Excitatory": ["SNAP25","STMN2","SYT1"], "Inhibitory": ["SNAP25","STMN2","SYT1"],
                      "Microglia": ["TREM2","CD68"], "Astrocyte": ["GJA1","SLC1A2","S100B"],
                      "Oligodendrocyte": ["MBP","PLP1","MOG","MAG"], "OPC": ["PDGFRA","OLIG1","CSPG4"],
                      "Endothelial": ["PECAM1","CLDN5","FLT1"], "Pericyte": ["PDGFRB","RGS5","VTN"]}
SELF_H = {"Inhibitory Neuron": "Inhibitory", "Excitatory  Neuron": "Excitatory"}
r2 = run(Xn, a.var_names, a.obs["Cell_Type"].astype(str).values, SELF_H,
         ["Inhibitory Neuron", "Excitatory  Neuron"], "GSE330130")
print("\nGSE330130 τ 敏感性:")
print(pd.DataFrame(r2).to_string(index=False))
out = pd.concat([pd.DataFrame(r1), pd.DataFrame(r2)], ignore_index=True)
out.to_csv("/Users/apple/Desktop/颅脑损伤专病数据库/P2_离子通道组景观论文/AnnoAudit_survey/crossval_tau_sensitivity.csv", index=False)
