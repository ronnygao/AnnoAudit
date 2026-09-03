#!/usr/bin/env python3
# ============================================================
# v16_methodology_enhancement: 方法学增强四方向
#   方向1: C2 判别力量化 — margin 对金标准标签的 AUROC + 混淆矩阵
#          + bootstrap 95% CI + per-cluster/per-cell-type 错分率
#   方向2: 独立 marker 生物学锚定 — 16 基因（与 32-gene panel 零重叠）
#          score 后独立 margin，与 C2 margin 相关 + 二分一致性
#   方向3: 标签置换检验 — AUROC 零分布 + p 值（n_perm=1000）
#   方向4: 跨数据集/跨物种迁移 — 32 基因同源映射 + 两数据集基因覆盖率
#          + 判别力（C2/AUROC）对比
#
# 输入: v15 落盘的 margin CSV、v14 per-cell 标签（CEREBRI c3_clusters.csv）、
#       CEREBRI 原始矩阵（方向2 锚定）、GSE330130 h5ad（重算定型+锚定）
# 输出: v16_methodology/summary.json + 各方向明细 CSV
# 用法: python3 v16_methodology_enhancement.py
# ============================================================
import os, glob, gzip, json, time, warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import scipy.io, scipy.sparse as sp
import anndata as ad
import scanpy as sc
from sklearn.metrics import roc_auc_score, confusion_matrix

T0 = time.time()
def log(m): print(f"[{time.time()-T0:7.1f}s] {m}", flush=True)

BASE = "/Users/apple/Desktop/颅脑损伤专病数据库/CEREBRI_KCNC3生信分析"
OUT = os.path.join(BASE, "v16_methodology")
os.makedirs(OUT, exist_ok=True)

MARKERS_MOUSE = {
    "Excitatory": ["Slc17a7","Rbfox3","Neurod6","Tbr1","Camk2a"],
    "Inhibitory": ["Gad1","Gad2","Sst","Pvalb","Npy"],
    "Astrocyte": ["Gfap","Aqp4","Slc1a3","Aldh1l1"],
    "Microglia": ["C1qb","Tyrobp","Cx3cr1","Hexb","Aif1"],
    "Oligodendrocyte": ["Mbp","Plp1","Mog","Mag"],
    "OPC": ["Pdgfra","Olig1","Cspg4"],
    "Endothelial": ["Pecam1","Cldn5","Flt1"],
    "Pericyte": ["Pdgfrb","Rgs5","Vtn"],
}
MARKERS_HUMAN = {k: [g.upper() for g in v] for k, v in MARKERS_MOUSE.items()}
NEU_TYPES = ["Excitatory", "Inhibitory"]
GLIA_TYPES = [t for t in MARKERS_MOUSE if t not in NEU_TYPES]

# ---- 独立锚定 marker（16 基因，逐一核对与 32-gene panel 零重叠）----
IND_MOUSE = {
    "NEU": ["Nrgn","Snap25","Stmn2","Syt1","Grin1"],        # 泛神经元
    "AST": ["S100b","Glul","Gja1"],                          # 星形胶质
    "MIC": ["Itgam","Ctss","Lgals3"],                        # 小胶质
    "OLI": ["Cnp","Olig2","Sox10"],                          # 少突胶质
    "END": ["Vwf"],                                          # 内皮
    "PERI": ["Kcnj8"],                                       # 周细胞
}
IND_HUMAN = {k: [g.upper() for g in v] for k, v in IND_MOUSE.items()}
IND_NEU = ["NEU"]
IND_GLIA = [t for t in IND_MOUSE if t not in IND_NEU]

# 32-gene panel 全集，用于零重叠断言
PANEL_ALL = {g.lower() for gs in MARKERS_MOUSE.values() for g in gs}

# ---- 同源映射（32 基因，1:1 直系同源；全部为保守大小写差异）----
ORTHO = {}
for m, h in zip([g for gs in MARKERS_MOUSE.values() for g in gs],
                [g for gs in MARKERS_HUMAN.values() for g in gs]):
    ORTHO[m] = h
assert len(ORTHO) == 32, "同源映射应为 32 对"

# 零重叠断言
for ct, gs in IND_MOUSE.items():
    for g in gs:
        assert g.lower() not in PANEL_ALL, f"重叠基因! {g} 已在 panel"

# ------------------------------------------------------------
def margin_from_scores(S, glial_types):
    neu_max = S[[f"mod_{t}" for t in NEU_TYPES]].max(axis=1).values
    glia_max = S[[f"mod_{t}" for t in glial_types]].max(axis=1).values
    return glia_max - neu_max

def ind_margin_from_scores(S):
    neu = S["ind_NEU"].values
    glia = S[[f"ind_{t}" for t in IND_GLIA]].max(axis=1).values
    return glia - neu

def binary_label(types, type_series):
    """Excitatory/Inhibitory -> 0 (神经元); 其余 8 类 -> 1 (胶质); Unknown -> -1 (剔除)"""
    neu = set(NEU_TYPES)
    lab = np.where(type_series.isin(neu), 0,
           np.where(type_series.isin(types), 1, -1))
    return lab.astype(int)

def auroc_ci(margin, lab, n_boot=500, seed=0):
    m = margin[lab >= 0]; l = lab[lab >= 0]
    if len(np.unique(l)) < 2 or len(m) < 10:
        return None
    auc = float(roc_auc_score(l, m))
    rng = np.random.default_rng(seed)
    n = len(l)
    boots = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        boots[i] = roc_auc_score(l[idx], m[idx])
    return auc, float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))

def perm_test(margin, lab, n_perm=1000, seed=1):
    m = margin[lab >= 0]; l = lab[lab >= 0].copy()
    auc = float(roc_auc_score(l, m))
    rng = np.random.default_rng(seed)
    null = np.empty(n_perm)
    for i in range(n_perm):
        null[i] = roc_auc_score(rng.permutation(l), m)
    p = float((np.abs(null - 0.5) >= np.abs(auc - 0.5)).mean())
    return auc, p, float(null.mean()), float(null.std())

def conf_at_tau(margin, lab, tau=1.5):
    m = margin[lab >= 0]; l = lab[lab >= 0]
    pred = (m > tau).astype(int)
    tn, fp, fn, tp = confusion_matrix(l, pred).ravel()
    sens = tp / max(tp + fn, 1); spec = tn / max(tn + fp, 1)
    bal = (sens + spec) / 2
    return {"TP": int(tp), "FP": int(fp), "FN": int(fn), "TN": int(tn),
            "sensitivity": round(sens * 100, 2), "specificity": round(spec * 100, 2),
            "balanced_accuracy": round(bal * 100, 2),
            "C2_pct": round((pred == 1).mean() * 100, 2)}

def per_group_rate(margin, group_series, tau=1.5, drop=("Unknown",)):
    df = pd.DataFrame({"margin": margin, "grp": group_series.values})
    df = df[~df["grp"].isin(drop)]
    g = df.groupby("grp").agg(n=("margin", "size"),
                              C2_glia_pct=("margin", lambda s: (s > tau).mean() * 100),
                              median_margin=("margin", "median"))
    return g.round(2)

# ------------------------------------------------------------
def run_cerebri_disk_only():
    """方向1+3（CEREBRI）：纯磁盘，不加载矩阵"""
    marg = pd.read_csv(os.path.join(BASE, "v15_tau_sensitivity/cerebri_margin_distribution.csv"))
    margin = marg["margin_glia_minus_neu"].values.astype(float)
    cl = pd.read_csv(os.path.join(BASE, "v14_cerebri_final/c3_clusters.csv"))
    assert len(cl) == len(margin), f"CEREBRI 行数不一致: {len(cl)} vs {len(margin)}"
    res = {"dataset": "CEREBRI", "n_cells": len(margin)}

    lab_marker = binary_label(MARKERS_MOUSE, cl["marker_type"])
    lab_cluster = binary_label(MARKERS_MOUSE, cl["cluster_type"])
    res["gold_standard_A_marker_typing_n_usable"] = int((lab_marker >= 0).sum())
    res["gold_standard_B_cluster_typing_n_usable"] = int((lab_cluster >= 0).sum())

    # 方向1: AUROC + CI（两个金标准）
    for gs_name, lab in [("marker_typing", lab_marker), ("cluster_typing", lab_cluster)]:
        ac = auroc_ci(margin, lab)
        if ac:
            res[f"AUROC_vs_{gs_name}"] = {"auc": round(ac[0], 4),
                                          "CI95": [round(ac[1], 4), round(ac[2], 4)]}
        else:
            res[f"AUROC_vs_{gs_name}"] = None
        log(f"CEREBRI AUROC vs {gs_name}: {ac}")

    # 方向1: 混淆矩阵 @ tau=1.5（金标准 A）
    res["confusion_tau1.5_vs_marker_typing"] = conf_at_tau(margin, lab_marker, 1.5)
    log(f"CEREBRI 混淆 @1.5: {res['confusion_tau1.5_vs_marker_typing']}")

    # 方向1: per-cluster 错分率（C2 判定 glia 的细胞在各 C3 cluster 中占比）
    res["per_cluster_C2_rate"] = per_group_rate(margin, cl["cluster"], 1.5).to_dict("index")
    log("CEREBRI per-cluster C2 rate:\n" + str(res["per_cluster_C2_rate"]))

    # 方向3: 置换检验（金标准 A）
    pt = perm_test(margin, lab_marker)
    res["permutation_test_vs_marker_typing"] = {"observed_AUROC": round(pt[0], 4),
                                                "p": pt[1], "null_mean": round(pt[2], 4),
                                                "null_sd": round(pt[3], 4)}
    log(f"CEREBRI 置换检验: AUROC={pt[0]:.4f} p={pt[1]:.2e}")
    return res

# ------------------------------------------------------------
def load_and_normalize_cerebri():
    RAW_DIR = "/Users/apple/Workbuddy/2026-07-22-16-57-18/cerebri_data/raw"
    ANNOT_CSV = os.path.join(BASE, "v3分析结果/cell_annotations.csv")
    ann = pd.read_csv(ANNOT_CSV)
    ann["bc_plain"] = ann["barcode"].str.replace(r"_GSM\d+$", "", regex=True)
    sel = ann[ann["cell_type"].astype(str) == "Glutamatergic_neuron"].copy()
    qc = (sel["n_genes"] > 200) & (sel["mito_pct"] < 20)
    sel = sel[qc].reset_index(drop=True)
    feat_files = sorted(os.listdir(RAW_DIR))
    feat_path = next(f for f in feat_files if f.endswith("_features.tsv.gz"))
    with gzip.open(os.path.join(RAW_DIR, feat_path), "rt") as f:
        features = [l.strip().split("\t") for l in f if not l.startswith("#")]
    gene_names = [x[1] if len(x) > 1 else x[0] for x in features]
    X_parts = []
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
    X = sp.vstack(X_parts).tocsr()
    del X_parts
    sym = pd.Index(gene_names)
    dup = sym.duplicated(keep=False)
    if dup.any():
        dup_pos = np.where(dup)[0]
        sums = np.asarray(X[:, dup].sum(axis=0)).ravel()
        order = pd.Series(sums, index=dup_pos).groupby(sym[dup].values).idxmax()
        keep_cols = np.zeros(len(sym), dtype=bool)
        keep_cols[np.asarray(order.values, dtype=int)] = True
        keep_cols[~dup] = True
        X = X[:, keep_cols].tocsr()
        genes = sym[keep_cols].values
    else:
        genes = sym.values
    log(f"CEREBRI 全基因矩阵 {X.shape}")
    umi = np.asarray(X.sum(axis=1)).ravel()
    X_cpm = X.multiply(1e6 / np.maximum(umi, 1)[:, None]).tocsr()
    X_cpm.data = np.log1p(X_cpm.data)
    A = ad.AnnData(X=X_cpm)
    A.obs_names = [f"c{i}" for i in range(X.shape[0])]
    A.var_names = list(genes)
    return A

# ------------------------------------------------------------
def run_anchor(A, ind_markers, tag):
    """方向2：独立 marker score + 与 C2 margin 的一致性"""
    cov = {}
    for ct, gs in ind_markers.items():
        use = [g for g in gs if g in A.var_names]
        cov[ct] = f"{len(use)}/{len(gs)}"
        log(f"{tag} 独立锚定 {ct}: {len(use)}/{len(gs)} 基因在数据中")
        sc.tl.score_genes(A, gene_list=use, score_name=f"ind_{ct}")
    S = A.obs[[c for c in A.obs.columns if c.startswith("ind_")]]
    ind_marg = ind_margin_from_scores(S)
    return ind_marg, cov

def anchor_vs_c2(margin, ind_marg, lab, tag):
    valid = lab >= 0
    mm = margin[valid]; ii = ind_marg[valid]; ll = lab[valid]
    r = np.corrcoef(mm, ii)[0, 1]
    auc = roc_auc_score(ll, ii)
    agree = np.sign(mm) == np.sign(ii)
    res = {"pearson_margin_corr": round(float(r), 4),
           "AUROC_independent_markers": round(float(auc), 4),
           "sign_agreement_pct": round(float(agree.mean() * 100), 2),
           "n_usable": int(valid.sum())}
    log(f"{tag} 锚定: corr={r:.4f} AUROC={auc:.4f} 符号一致率={agree.mean()*100:.1f}%")
    return res

# ------------------------------------------------------------
def run_gse330130():
    H5 = "/Users/apple/Desktop/TBI分子时钟/数据/GSE330130_LSC_SNRNASEQ_COUNTS.h5ad"
    NEURON_LABELS = ["Inhibitory Neuron", "Excitatory  Neuron", "Motor Neurons"]
    a = ad.read_h5ad(H5)
    log(f"GSE330130 loaded {a.shape}")
    a.var["symbol"] = a.var["common_name"].astype(str)
    a.var_names = a.var["symbol"].values
    a = a[:, ~pd.Index(a.var_names).duplicated(keep="first")].copy()
    X = a.X.tocsr()
    umi = np.asarray(X.sum(axis=1)).ravel()
    Xn = X.multiply(1e6 / np.maximum(umi, 1)[:, None]).tocsr()
    Xn.data = np.log1p(Xn.data)
    a.X = Xn
    del X, Xn
    off = a.obs["Cell_Type"].isin(NEURON_LABELS).values
    sub = a[off].copy()
    del a
    log(f"GSE 神经元子集 {sub.shape}")

    # 重算 8 模块 → margin（与 v15 CSV 交叉验证）+ per-cell marker 定型
    np.random.seed(0)
    for ct, genes in MARKERS_HUMAN.items():
        use = [g for g in genes if g in sub.var_names]
        sc.tl.score_genes(sub, gene_list=use, score_name=f"mod_{ct}")
    S = sub.obs[[c for c in sub.obs.columns if c.startswith("mod_")]]
    margin = margin_from_scores(S, GLIA_TYPES)
    v15_marg = pd.read_csv(os.path.join(BASE, "v15_tau_sensitivity/gse330130_margin_distribution.csv"))["margin_glia_minus_neu"].values.astype(float)
    assert len(margin) == len(v15_marg)
    maxdiff = float(np.abs(margin - v15_marg).max())
    assert maxdiff < 1e-6, f"GSE margin 重算与 v15 不一致: maxdiff={maxdiff}"
    log(f"GSE margin 重算与 v15 完全一致 ✓ (maxdiff={maxdiff:.2e})")

    # per-cell marker 定型（v14 C1 逻辑: argmax 且 max>0, 否则 Unknown）
    sc_ = S.values
    names = list(MARKERS_HUMAN.keys())
    assert sc_.shape[1] == len(names), f"模块列数 {sc_.shape[1]} != {len(names)}"
    idx = np.argmax(sc_, axis=1)
    mx = sc_[np.arange(len(sc_)), idx]
    lab_s = [names[idx[i]] if mx[i] > 0 else "Unknown" for i in range(len(sc_))]
    mtype = pd.Series(lab_s)
    lab_marker = binary_label(MARKERS_MOUSE, mtype)
    res = {"dataset": "GSE330130", "n_cells": int(sub.n_obs),
           "gold_standard_A_marker_typing_n_usable": int((lab_marker >= 0).sum())}

    # 方向1: AUROC + CI
    ac = auroc_ci(margin, lab_marker)
    res["AUROC_vs_marker_typing"] = {"auc": round(ac[0], 4), "CI95": [round(ac[1], 4), round(ac[2], 4)]}
    res["confusion_tau1.5_vs_marker_typing"] = conf_at_tau(margin, lab_marker, 1.5)
    log(f"GSE AUROC={ac[0]:.4f} 混淆@1.5={res['confusion_tau1.5_vs_marker_typing']}")

    # 方向1: per-cell-type 错分率（官方 3 类神经元）
    res["per_celltype_C2_rate"] = per_group_rate(margin, sub.obs["Cell_Type"], 1.5).to_dict("index")
    log("GSE per-celltype C2 rate:\n" + str(res["per_celltype_C2_rate"]))

    # 方向3: 置换检验
    pt = perm_test(margin, lab_marker)
    res["permutation_test_vs_marker_typing"] = {"observed_AUROC": round(pt[0], 4),
                                                "p": pt[1], "null_mean": round(pt[2], 4),
                                                "null_sd": round(pt[3], 4)}
    log(f"GSE 置换检验: AUROC={pt[0]:.4f} p={pt[1]:.2e}")

    # 方向2: 独立锚定
    ind_marg, cov = run_anchor(sub, IND_HUMAN, "GSE")
    res["anchor_gene_coverage"] = cov
    res["anchor_vs_C2"] = anchor_vs_c2(margin, ind_marg, lab_marker, "GSE")
    return res, sub

# ------------------------------------------------------------
def run_transfer(cerebri_genes, gse_genes, cerebri_res, gse_res):
    """方向4：跨数据集/跨物种迁移"""
    # 32 基因同源对在各自数据集的覆盖率
    m_in_cerebri = {m: (m in cerebri_genes) for m in ORTHO}
    h_in_gse = {h: (h in gse_genes) for h in ORTHO.values()}
    cov_c = sum(m_in_cerebri.values())
    cov_g = sum(h_in_gse.values())
    common = sum(1 for m, h in ORTHO.items() if m_in_cerebri[m] and h_in_gse[h])
    res = {
        "orthologous_pairs": 32,
        "mouse_genes_present_in_CEREBRI": f"{cov_c}/32",
        "human_orthologs_present_in_GSE330130": f"{cov_g}/32",
        "pairs_present_in_both": f"{common}/32",
        "missing_in_CEREBRI": [m for m, ok in m_in_cerebri.items() if not ok],
        "missing_in_GSE": [h for h, ok in h_in_gse.items() if not ok],
        "CEREBRI": {"C2_tau1.5": cerebri_res["confusion_tau1.5_vs_marker_typing"]["C2_pct"],
                    "AUROC_vs_marker_typing": cerebri_res["AUROC_vs_marker_typing"]["auc"]},
        "GSE330130": {"C2_tau1.5": gse_res["confusion_tau1.5_vs_marker_typing"]["C2_pct"],
                      "AUROC_vs_marker_typing": gse_res["AUROC_vs_marker_typing"]["auc"]},
    }
    log(f"跨数据集迁移: 覆盖率 CEREBRI {cov_c}/32, GSE {cov_g}/32, 两者 {common}/32")
    log(f"判别力对比: CEREBRI C2={res['CEREBRI']['C2_tau1.5']}% AUROC={res['CEREBRI']['AUROC_vs_marker_typing']} | "
        f"GSE C2={res['GSE330130']['C2_tau1.5']}% AUROC={res['GSE330130']['AUROC_vs_marker_typing']}")
    return res

# ------------------------------------------------------------
def main():
    summary = {}
    # ---- CEREBRI: 方向1+3（纯磁盘）----
    cerebri_res = run_cerebri_disk_only()
    summary["CEREBRI"] = cerebri_res

    # ---- CEREBRI: 方向2 独立锚定（需加载矩阵）----
    A = load_and_normalize_cerebri()
    ind_marg, cov = run_anchor(A, IND_MOUSE, "CEREBRI")
    cerebri_res["anchor_gene_coverage"] = cov
    cl = pd.read_csv(os.path.join(BASE, "v14_cerebri_final/c3_clusters.csv"))
    margin = pd.read_csv(os.path.join(BASE, "v15_tau_sensitivity/cerebri_margin_distribution.csv"))["margin_glia_minus_neu"].values.astype(float)
    lab_marker = binary_label(MARKERS_MOUSE, cl["marker_type"])
    cerebri_res["anchor_vs_C2"] = anchor_vs_c2(margin, ind_marg, lab_marker, "CEREBRI")
    cerebri_genes = set(A.var_names)
    del A

    # ---- GSE330130: 方向1+2+3 ----
    gse_res, sub = run_gse330130()
    summary["GSE330130"] = gse_res
    gse_genes = set(sub.var_names)
    del sub

    # ---- 方向4: 跨数据集迁移 ----
    summary["cross_dataset_transfer"] = run_transfer(cerebri_genes, gse_genes, cerebri_res, gse_res)

    with open(os.path.join(OUT, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
    log("V16 DONE -> " + os.path.join(OUT, "summary.json"))

if __name__ == "__main__":
    main()
