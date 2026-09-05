#!/usr/bin/env python3
# ============================================================
# survey_cerebri_all_labels.py — CEREBRI 官方全标签 × C1 marker 审计
# 横断面调查（AnnoAudit survey）的"深度第一端点"：
#   对官方注释全部 11 个非 Unknown 标签逐一审计 marker 组成，
#   回答审稿人问题："污染是只集中在谷氨酸能标签，还是系统性？"
# 协议 = v14 C1（32 基因小鼠 panel, 8 类 argmax, 阈值 mean log1p(CPM)>0）
# 优化：仅提取 32 marker 基因的列，340K 细胞全量可一次载入
# ============================================================
import os, glob, gzip, time, warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import scipy.io, scipy.sparse as sp

T0 = time.time()
def log(m): print(f"[{time.time()-T0:7.1f}s] {m}", flush=True)

RAW_DIR = "/Users/apple/Workbuddy/2026-07-22-16-57-18/cerebri_data/raw"
ANNOT_CSV = "/Users/apple/Desktop/颅脑损伤专病数据库/CEREBRI_KCNC3生信分析/v3分析结果/cell_annotations.csv"
OUT_DIR = "/Users/apple/Desktop/颅脑损伤专病数据库/P2_离子通道组景观论文/AnnoAudit_survey/cerebri_all_labels"
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
NEU = ["Excitatory", "Inhibitory"]
assert sum(len(v) for v in MARKERS.values()) == 32

# ---------- 注释 ----------
ann = pd.read_csv(ANNOT_CSV)
ann["bc_plain"] = ann["barcode"].str.replace(r"_GSM\d+$", "", regex=True)
ann = ann.reset_index(drop=True)
log(f"注释 {len(ann)} 行; 官方标签分布:")
log(ann["cell_type"].astype(str).value_counts().to_string().replace("\n", " | "))

# ---------- features 全基因表 ----------
feat_path = sorted(glob.glob(os.path.join(RAW_DIR, "*_features.tsv.gz")))[0]
with gzip.open(feat_path, "rt") as f:
    feats = [l.strip().split("\t") for l in f if not l.startswith("#")]
gene_names = [x[1] if len(x) > 1 else x[0] for x in feats]
gene_idx = {g: i for i, g in enumerate(gene_names)}
mk_idx = {}
missing = []
for ct, genes in MARKERS.items():
    for g in genes:
        if g in gene_idx:
            mk_idx[g] = gene_idx[g]
        else:
            missing.append(g)
log(f"panel 基因命中 {len(mk_idx)}/32; 缺失: {missing if missing else '无'}")

# 独立验证基因（不在 32 panel 内）：微胶质/巨噬、星形、神经元
VERIFY = {"Mg_extra": ["Trem2", "Cd68", "Lyz2", "Apoe"],
          "Astro_extra": ["Gja1", "Slc1a2", "S100b"],
          "Neuron_extra": ["Snap25", "Stmn2", "Syt1"]}
vfy_idx = {}
for grp, genes in VERIFY.items():
    for g in genes:
        if g in gene_idx:
            vfy_idx[g] = gene_idx[g]
        else:
            log(f"!! 验证基因缺失: {g}")
log(f"验证基因命中 {len(vfy_idx)}/{sum(len(v) for v in VERIFY.values())}")

# ---------- 逐 GSM 提取 32 列 ----------
gsm_list = sorted(ann["gsm"].unique())
X_parts, bc_parts = [], []
for gi, gsm in enumerate(gsm_list):
    mtx_f = glob.glob(os.path.join(RAW_DIR, f"{gsm}_*_matrix.mtx.gz"))
    bar_f = glob.glob(os.path.join(RAW_DIR, f"{gsm}_*_barcodes.tsv.gz"))
    if not mtx_f or not bar_f:
        log(f"!! {gsm} 文件缺失，跳过"); continue
    with gzip.open(mtx_f[0], "rb") as f:
        mat = scipy.io.mmread(f)
    mat = mat.tocsr()
    # 只取 32 marker + 验证基因行 -> 转置为 cells x (32+n)
    sel_rows = sorted(set(mk_idx.values()) | set(vfy_idx.values()))
    sub = mat[sel_rows, :].T.tocsr().astype(np.float32)
    del mat
    with gzip.open(bar_f[0], "rt") as f:
        barcodes = [l.strip() for l in f]
    X_parts.append(sub)
    bc_parts.append(pd.DataFrame({"gsm": gsm, "barcode": barcodes}))
    log(f"[{gi+1}/{len(gsm_list)}] {gsm}: {sub.shape[0]} 细胞 x {sub.shape[1]} 基因列")
X = sp.vstack(X_parts).tocsr()
bcm = pd.concat(bc_parts, ignore_index=True)
del X_parts, bc_parts
log(f"全矩阵细胞 {X.shape[0]} vs 注释 {len(ann)}")

# ---------- log1p CPM 归一 + 8 类打分 ----------
umi = np.asarray(X.sum(axis=1)).ravel()
Xn = X.multiply(1e6 / np.maximum(umi, 1)[:, None]).tocsr()
Xn.data = np.log1p(Xn.data)
# 矩阵列 -> 全基因索引 -> 基因名（32 panel + 验证基因混合列）
inv_idx = {i: g for g, i in gene_idx.items()}
col_genes = [inv_idx[i] for i in sel_rows]
gpos = {g: j for j, g in enumerate(col_genes)}   # 压缩矩阵第 j 列 = col_genes[j]
vpos = {g: j for j, g in enumerate(col_genes) if g in vfy_idx}
scores = {}
for ct, genes in MARKERS.items():
    cols = [gpos[g] for g in genes if g in gpos]
    scores[ct] = np.asarray(Xn[:, cols].mean(axis=1)).ravel()
sdf = pd.DataFrame(scores)
best = sdf.idxmax(axis=1)
bval = sdf.max(axis=1)
marker_type = np.where(bval > 0, best.values, "Unknown")
log("marker 全图谱定型分布: " + pd.Series(marker_type).value_counts().to_string().replace("\n", " | "))

# ---------- 对齐注释（gsm + barcode 复合键，防跨 GSM barcode 重复）----------
bcm["marker_type"] = marker_type
mrg = ann.merge(bcm, left_on=["gsm", "bc_plain"], right_on=["gsm", "barcode"], how="left")
ann["marker_type"] = mrg["marker_type"].fillna("Unmatched").values
log(f"注释与矩阵匹配率: {(ann['marker_type'] != 'Unmatched').mean()*100:.2f}%")

# ---------- 混淆矩阵 + 每标签指标 ----------
valid = ann[ann["cell_type"].astype(str) != "Unknown"].copy()
conf = pd.crosstab(valid["cell_type"].astype(str), valid["marker_type"])
conf.to_csv(os.path.join(OUT_DIR, "confusion_matrix.csv"))
log("\n混淆矩阵 (官方标签 x marker 类型):\n" + conf.to_string())

SELF = {"Glutamatergic_neuron": "Excitatory", "GABAergic_neuron": "Inhibitory",
        "Astrocyte": "Astrocyte", "Microglia": "Microglia",
        "Oligodendrocyte": "Oligodendrocyte", "OPC": "OPC",
        "Endothelial": "Endothelial", "Pericyte": "Pericyte",
        "VSM": "Pericyte", "Fibroblast": "Pericyte", "Erythrocyte": "Unknown"}
rows = []
for lab in valid["cell_type"].astype(str).unique():
    sub = valid[valid["cell_type"].astype(str) == lab]
    n = len(sub)
    vc = sub["marker_type"].value_counts()
    self_t = SELF.get(lab)
    self_frac = vc.get(self_t, 0) / n if self_t and self_t != "Unknown" else np.nan
    neu_frac = (sub["marker_type"].isin(NEU)).mean()
    rows.append({"official_label": lab, "n": n,
                 "self_marker_frac": round(self_frac * 100, 1) if not np.isnan(self_frac) else None,
                 "neuron_frac": round(neu_frac * 100, 1),
                 "nonneuron_frac_if_neuron_label": round((1 - neu_frac) * 100, 1),
                 "top_marker": vc.index[0], "top_frac": round(vc.iloc[0] / n * 100, 1),
                 "second_marker": vc.index[1] if len(vc) > 1 else "", "second_frac": round(vc.iloc[1] / n * 100, 1) if len(vc) > 1 else None})
res = pd.DataFrame(rows).sort_values("n", ascending=False)
res.to_csv(os.path.join(OUT_DIR, "label_audit.csv"), index=False)
log("\n每官方标签审计:\n" + res.to_string(index=False))

# ---------- 独立验证基因表达（官方标签 × marker 亚型；证伪打分假象）----------
# 若官方 Astrocyte/OPC 标签中"被判 Microglia"的细胞真为微胶质（非打分假象），
# 其 Trem2/Cd68/Lyz2/Apoe 应高表达、且 Gfap/Gja1 应低表达，与官方 Microglia 标签一致。
bcm2 = bcm.reset_index().rename(columns={"index": "row_in_X"})
mrg2 = ann.merge(bcm2, left_on=["gsm", "bc_plain"], right_on=["gsm", "barcode"], how="left")
Xsub = Xn[mrg2["row_in_X"].values]          # 行序 = ann 序, 列 = col_genes
vcols = {g: j for j, g in enumerate(col_genes) if g in vfy_idx}
for grpname, glist in VERIFY.items():
    present = [g for g in glist if g in vcols]
    if not present:
        log(f"!! 验证组 {grpname} 无基因命中"); continue
    col_pos = [vcols[g] for g in present]
    expr = pd.DataFrame(np.asarray(Xsub[:, col_pos].todense()), columns=present)
    ann2 = ann.reset_index(drop=True)
    tab = pd.concat([ann2[["cell_type", "marker_type"]].astype(str), expr], axis=1)
    piv = tab.groupby(["cell_type", "marker_type"], observed=True)[present].mean().round(3)
    piv.to_csv(os.path.join(OUT_DIR, f"verify_{grpname}_expression.csv"))
log("独立验证基因表达表已输出 (verify_*_expression.csv)")

# 汇总 json
import json
summary = {"dataset": "CEREBRI GSE269748", "total_cells": int(len(ann)),
           "matched_pct": round(float((ann["marker_type"] != "Unmatched").mean() * 100), 2),
           "official_labels": [{"label": r["official_label"], "n": int(r["n"]),
                                "self_marker_pct": r["self_marker_frac"],
                                "neuron_marker_pct": r["neuron_frac"],
                                "top_marker": r["top_marker"], "top_pct": r["top_frac"]} for r in rows]}
with open(os.path.join(OUT_DIR, "summary.json"), "w") as f:
    json.dump(summary, f, indent=2, ensure_ascii=False)
log("\nDONE — 输出: " + OUT_DIR)
