#!/usr/bin/env python3
# ============================================================
# survey_gse330130_all_labels.py — GSE330130 (人 ALS 脊髓) 官方全标签 × C1 审计
# 横断面调查（AnnoAudit survey）第二端点；协议与 survey_cerebri_all_labels.py 完全一致
# 输入：GEO 附件 h5ad（含官方 Cell_Type 注释）；人类 32 基因直系同源 panel
# 独立验证基因：微胶质 Trem2/CD68/LYZ2/APOE、星形 GJA1/SLC1A2/S100B、
#              神经元 SNAP25/STMN2/SYT1、胆碱能 CHAT/SLC5A7（检验 Motor Neurons 身份）
# ============================================================
import os, time, warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import anndata as ad

T0 = time.time()
def log(m): print(f"[{time.time()-T0:7.1f}s] {m}", flush=True)

H5 = "/Users/apple/Desktop/颅脑损伤专病数据库/P2_离子通道组景观论文/AnnoAudit_survey/data/GSE330130_LSC_SNRNASEQ_COUNTS.h5ad"
OUT_DIR = "/Users/apple/Desktop/颅脑损伤专病数据库/P2_离子通道组景观论文/AnnoAudit_survey/gse330130_all_labels"
os.makedirs(OUT_DIR, exist_ok=True)

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
NEU = ["Excitatory", "Inhibitory"]
assert sum(len(v) for v in MARKERS_H.values()) == 32
VERIFY = {"Mg_extra": ["TREM2", "CD68", "LYZ2", "APOE"],
          "Astro_extra": ["GJA1", "SLC1A2", "S100B"],
          "Neuron_extra": ["SNAP25", "STMN2", "SYT1"],
          "Cholinergic": ["CHAT", "SLC5A7"]}

a = ad.read_h5ad(H5)
a.var["symbol"] = a.var["common_name"].astype(str)
a.var_names = a.var["symbol"].values
a = a[:, ~pd.Index(a.var_names).duplicated(keep="first")].copy()
log(f"loaded {a.shape}")

missing = [g for gs in MARKERS_H.values() for g in gs if g not in a.var_names]
log(f"panel 缺失: {missing if missing else '无 (32/32)'}")
vmiss = [g for gs in VERIFY.values() for g in gs if g not in a.var_names]
log(f"验证基因缺失: {vmiss if vmiss else '无'}")

# ---------- 归一化 log1p CPM ----------
X = a.X.tocsr().astype(np.float32)
umi = np.asarray(X.sum(axis=1)).ravel()
Xn = X.multiply(1e6 / np.maximum(umi, 1)[:, None]).tocsr()
Xn.data = np.log1p(Xn.data)

gidx = {g: a.var_names.get_loc(g) for gs in MARKERS_H.values() for g in gs if g in a.var_names}
vidx = {g: a.var_names.get_loc(g) for gs in VERIFY.values() for g in gs if g in a.var_names}

scores = {}
for ct, genes in MARKERS_H.items():
    cols = [gidx[g] for g in genes if g in gidx]
    scores[ct] = np.asarray(Xn[:, cols].mean(axis=1)).ravel()
sdf = pd.DataFrame(scores)
best = sdf.idxmax(axis=1)
bval = sdf.max(axis=1)
marker_type = np.where(bval > 0, best.values, "Unknown")
log("marker 定型分布: " + pd.Series(marker_type).value_counts().to_string().replace("\n", " | "))

a.obs["marker_type"] = marker_type
ann = a.obs[["Cell_Type"]].astype(str).reset_index()
ann["marker_type"] = marker_type
ann = ann.rename(columns={"Cell_Type": "official_label"})

# ---------- 混淆矩阵 ----------
conf = pd.crosstab(ann["official_label"], ann["marker_type"])
conf.to_csv(os.path.join(OUT_DIR, "confusion_matrix.csv"))
log("\n混淆矩阵:\n" + conf.to_string())

SELF = {"Oligodendrocytes": "Oligodendrocyte", "Astrocyte": "Astrocyte",
        "Microglia": "Microglia", "OPCs": "OPC", "Endothelial": "Endothelial",
        "Inhibitory Neuron": "Inhibitory", "Excitatory  Neuron": "Excitatory",
        "Motor Neurons": None, "Meninges": None}
rows = []
for lab in ann["official_label"].unique():
    sub = ann[ann["official_label"] == lab]
    n = len(sub)
    vc = sub["marker_type"].value_counts()
    self_t = SELF.get(lab)
    self_frac = vc.get(self_t, 0) / n if self_t else np.nan
    neu_frac = (sub["marker_type"].isin(NEU)).mean()
    rows.append({"official_label": lab, "n": n,
                 "self_marker_frac": round(self_frac * 100, 1) if not np.isnan(self_frac) else None,
                 "neuron_frac": round(neu_frac * 100, 1),
                 "top_marker": vc.index[0], "top_frac": round(vc.iloc[0] / n * 100, 1),
                 "second_marker": vc.index[1] if len(vc) > 1 else "", "second_frac": round(vc.iloc[1] / n * 100, 1) if len(vc) > 1 else None})
res = pd.DataFrame(rows).sort_values("n", ascending=False)
res.to_csv(os.path.join(OUT_DIR, "label_audit.csv"), index=False)
log("\n每官方标签审计:\n" + res.to_string(index=False))

# ---------- 独立验证基因表达 ----------
vcols = {g: j for j, g in enumerate(a.var_names) if g in vidx}
Xsub = Xn
for grpname, glist in VERIFY.items():
    present = [g for g in glist if g in vcols]
    if not present:
        log(f"!! {grpname} 无基因"); continue
    col_pos = [vcols[g] for g in present]
    expr = pd.DataFrame(np.asarray(Xsub[:, col_pos].todense()), columns=present)
    tab = pd.concat([ann[["official_label", "marker_type"]], expr], axis=1)
    piv = tab.groupby(["official_label", "marker_type"], observed=True)[present].mean().round(3)
    piv.to_csv(os.path.join(OUT_DIR, f"verify_{grpname}_expression.csv"))
log("独立验证基因表已输出")

import json
summary = {"dataset": "GSE330130 human ALS lumbar spinal cord", "total_cells": int(len(ann)),
           "official_labels": [{"label": r["official_label"], "n": int(r["n"]),
                                "self_marker_pct": r["self_marker_frac"],
                                "neuron_marker_pct": r["neuron_frac"],
                                "top_marker": r["top_marker"], "top_pct": r["top_frac"]} for r in rows]}
with open(os.path.join(OUT_DIR, "summary.json"), "w") as f:
    json.dump(summary, f, indent=2, ensure_ascii=False)
log("\nDONE — " + OUT_DIR)
