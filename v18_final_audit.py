#!/usr/bin/env python3
# ============================================================
# v18_final_audit.py — 正式版全标签审计 + 本型基准判别（定稿数字）
# 判别协议 D1（正式定义）：
#   对官方标签 L 的每细胞：t0 = argmax 型（32 panel 8 类, >0 阈值）
#   若 t0 ≠ 本型(L)：以独立基因做身份确认——比较"候选异型独立基因得分 s_o"
#   与"本型独立基因得分 s_s"：s_o > τ·s_s (τ=1.5) → 真错标(本质=t0型)
#   否则 → 身份不明（潜在 ambient/低质量，不计污染）
#   本型独立基因：Exc/Snap25-Stmn2-Syt1; Mg/Trem2-Cd68-Lyz2; Astro/Gja1-Slc1a2-S100b;
#                  Oligo/Mbp-Plp1-Mog-Mag(panel 已含,用 panel); Endo-Pericyte-OPC(panel 4marker 均值)
# 输出：每标签 [n, marker-self%, 判别污染%, 判别分解] —— 正式定稿表
# ============================================================
import os, glob, gzip, json, time, warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import scipy.io, scipy.sparse as sp
import anndata as ad

T0 = time.time()
def log(m): print(f"[{time.time()-T0:7.1f}s] {m}", flush=True)
OUT = "/Users/apple/Desktop/颅脑损伤专病数据库/P2_离子通道组景观论文/AnnoAudit_survey/final_audit"
os.makedirs(OUT, exist_ok=True)
TAU = 1.5

MK8 = ["Excitatory", "Inhibitory", "Astrocyte", "Microglia", "Oligodendrocyte", "OPC", "Endothelial", "Pericyte"]
PANEL = {
    "Excitatory": ["Slc17a7","Rbfox3","Neurod6","Tbr1","Camk2a"],
    "Inhibitory": ["Gad1","Gad2","Sst","Pvalb","Npy"],
    "Astrocyte": ["Gfap","Aqp4","Slc1a3","Aldh1l1"],
    "Microglia": ["C1qb","Tyrobp","Cx3cr1","Hexb","Aif1"],
    "Oligodendrocyte": ["Mbp","Plp1","Mog","Mag"],
    "OPC": ["Pdgfra","Olig1","Cspg4"],
    "Endothelial": ["Pecam1","Cldn5","Flt1"],
    "Pericyte": ["Pdgfrb","Rgs5","Vtn"],
}
# 本型独立基因（判别用，避开 panel 基因）
# 神经元标签统一用泛神经元独立基因（SNAP25/STMN2/SYT1）——审计问"是否神经元"(类级)，
# 非"兴奋/抑制亚型"；亚型 marker（GAD/SLC17A7）在 argmax_comp 描述层保留
NEU_IND_M = ["Snap25", "Stmn2", "Syt1"]; NEU_IND_H = ["SNAP25", "STMN2", "SYT1"]
INDEP = {"Excitatory": NEU_IND_M, "Inhibitory": NEU_IND_M,
         "Microglia": ["Trem2","Cd68","Lyz2"], "Astrocyte": ["Gja1","Slc1a2","S100b"],
         "Oligodendrocyte": ["Mbp","Plp1","Mog","Mag"], "OPC": ["Pdgfra","Olig1","Cspg4"],
         "Endothelial": ["Pecam1","Cldn5","Flt1"], "Pericyte": ["Pdgfrb","Rgs5","Vtn"]}
# 官方标签 → 本型 marker 类
SELF_MAP = {"Glutamatergic_neuron": "Excitatory", "GABAergic_neuron": "Inhibitory",
            "Astrocyte": "Astrocyte", "Microglia": "Microglia", "Oligodendrocyte": "Oligodendrocyte",
            "OPC": "OPC", "Endothelial": "Endothelial", "Pericyte": "Pericyte",
            "VSM": "Pericyte", "Fibroblast": "Pericyte", "Erythrocyte": None, "Unknown": None}
SELF_MAP_H = {"Oligodendrocytes": "Oligodendrocyte", "Astrocyte": "Astrocyte", "Microglia": "Microglia",
              "OPCs": "OPC", "Endothelial": "Endothelial", "Inhibitory Neuron": "Inhibitory",
              "Excitatory  Neuron": "Excitatory", "Motor Neurons": None, "Meninges": None}

def run_discrimination(obs_official, expr, var_names, self_map, dataset, extra_mn_chat=False):
    """expr: cells × genes(仅关注基因列, 列名=var_names); obs_official: 官方标签数组"""
    gpos = {g: j for j, g in enumerate(var_names)}
    def sm(genes):
        cols = [gpos[g] for g in genes if g in gpos]
        return np.asarray(expr[:, cols].mean(axis=1)).ravel() if cols else np.zeros(expr.shape[0])
    # 8 类 panel 得分（用完整 panel 基因的均值；已在 var_names 中）
    scs = {ct: sm(PANEL[ct]) for ct in MK8}
    sdf = pd.DataFrame(scs); amt = sdf.idxmax(axis=1).values; amv = sdf.max(axis=1).values
    amt = np.where(amv > 0, amt, "Unknown")
    ind_scores = {ct: sm(INDEP[ct]) for ct in MK8}
    rows = []
    for lab in pd.unique(obs_official):
        if lab not in self_map or self_map[lab] is None:
            continue
        self_t = self_map[lab]
        mask = obs_official == lab
        n = int(mask.sum())
        sub_amt = amt[mask]
        # marker-self 率
        self_frac = float((sub_amt == self_t).mean())
        # 判别: argmax 非本型 & s_o > τ·s_self_ind
        non_self = sub_amt != self_t
        # 特殊: 本型为 Excitatory 且 argmax=Inhibitory → 仍是神经元，不判污染（亚型混淆）
        neu_other = (sub_amt == "Inhibitory") & (self_t == "Excitatory")
        neu_other |= (sub_amt == "Excitatory") & (self_t == "Inhibitory")
        cand = non_self & ~neu_other
        s_self = ind_scores[self_t][mask]
        # 各候选异型的独立分（对 argmax 型）
        conf = np.zeros(n, dtype=bool)
        for o in MK8:
            if o == self_t:
                continue
            m_o = (sub_amt == o) & cand
            if m_o.any():
                s_o = ind_scores[o][mask]
                conf |= m_o & (s_o > TAU * s_self)
        # 神经元内亚型混淆：argmax 互换（Exc↔Inh）不计污染，但计入"神经元内"
        neu_mix = neu_other & (ind_scores[self_t][mask] > 0)
        contamin = conf.sum() / n
        comp = {}
        for o in MK8 + ["Unknown"]:
            v = int(((sub_amt == o) & cand & ~conf).sum())
            if o != self_t and v:
                comp[f"unresolved_{o}"] = v
        rows.append({"dataset": dataset, "official_label": lab, "n": n,
                     "marker_self_pct": round(self_frac * 100, 1),
                     "disc_contamination_pct": round(contamin * 100, 1),
                     "confirmed_self_pct": round((1 - contamin) * 100, 1),
                     "neu_mix_pct": round(neu_mix.sum() / n * 100, 1),
                     "argmax_comp": {k: round(float(v / n * 100), 1) for k, v in
                                     pd.Series(sub_amt).value_counts().items()}})
        log(f"{dataset} | {lab:<22} n={n:>7} | marker-self {self_frac*100:5.1f}% | 判别污染 {contamin*100:5.1f}%")
    return pd.DataFrame(rows)

# ================= CEREBRI =================
log("== CEREBRI ==")
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
res_c = run_discrimination(obs_x, Xcn, colg, SELF_MAP, "CEREBRI")
res_c.to_csv(os.path.join(OUT, "cerebri_final_audit.csv"), index=False)

# ================= GSE330130 =================
log("\n== GSE330130 ==")
H5 = "/Users/apple/Desktop/颅脑损伤专病数据库/P2_离子通道组景观论文/AnnoAudit_survey/data/GSE330130_LSC_SNRNASEQ_COUNTS.h5ad"
a = ad.read_h5ad(H5)
a.var["symbol"] = a.var["common_name"].astype(str)
a.var_names = a.var["symbol"].values
a = a[:, ~pd.Index(a.var_names).duplicated(keep="first")].copy()
X = a.X.tocsr().astype(np.float32)
umi = np.asarray(X.sum(axis=1)).ravel()
Xn = X.multiply(1e6 / np.maximum(umi, 1)[:, None]).tocsr()
Xn.data = np.log1p(Xn.data)
# 人类 panel 大写
PANEL_H = {k: [g.upper() for g in v] for k, v in PANEL.items()}
INDEP_H = {"Excitatory": NEU_IND_H, "Inhibitory": NEU_IND_H,
           "Microglia": ["TREM2","CD68"], "Astrocyte": ["GJA1","SLC1A2","S100B"],
           "Oligodendrocyte": ["MBP","PLP1","MOG","MAG"], "OPC": ["PDGFRA","OLIG1","CSPG4"],
           "Endothelial": ["PECAM1","CLDN5","FLT1"], "Pericyte": ["PDGFRB","RGS5","VTN"]}
# 复用 run_discrimination（覆盖全局 PANEL/INDEP）
globals()["PANEL"] = PANEL_H
globals()["INDEP"] = INDEP_H
res_h = run_discrimination(a.obs["Cell_Type"].astype(str).values, Xn, a.var_names, SELF_MAP_H, "GSE330130")
res_h.to_csv(os.path.join(OUT, "gse330130_final_audit.csv"), index=False)
print("\n===== 正式审计表 =====\n")
print(pd.concat([res_c, res_h], ignore_index=True).to_string(index=False))
log("DONE")
