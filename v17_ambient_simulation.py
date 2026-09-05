#!/usr/bin/env python3
# ============================================================
# v17_ambient_simulation.py — 判别协议校准：ambient 梯度仿真
# 目的：证明 C1 单通道在 ambient RNA 下高估污染，而"独立基因判别协议"
#       (argmax 异型细胞 → 比较本型 vs 异型独立基因) 可校正 ambient 高估。
# 合成官方标签 = Excitatory(真) + Microglia(真污染 c)，叠加 ambient 剂量 d：
#   X_obs = (1-d)·X_true + d·A          (A = ambient 谱: 微胶质主导 / 少突主导)
# 测量：C1 原始污染估计 vs 判别校正估计 vs 真 c
# ============================================================
import os, glob, gzip, time, json, warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import scipy.io, scipy.sparse as sp

T0 = time.time()
def log(m): print(f"[{time.time()-T0:7.1f}s] {m}", flush=True)

SIM_DIR = "/Users/apple/Desktop/颅脑损伤专病数据库/CEREBRI_KCNC3生信分析/v12_simulation"
OUT_DIR = "/Users/apple/Desktop/颅脑损伤专病数据库/P2_离子通道组景观论文/AnnoAudit_survey/v17_ambient_sim"
os.makedirs(OUT_DIR, exist_ok=True)

RAW_DIR = "/Users/apple/Workbuddy/2026-07-22-16-57-18/cerebri_data/raw"
ANNOT_CSV = "/Users/apple/Desktop/颅脑损伤专病数据库/CEREBRI_KCNC3生信分析/v3分析结果/cell_annotations.csv"
CLEVELS = [0.0, 0.05, 0.10, 0.20, 0.50, 0.80]
DLEVELS = [0.0, 0.10, 0.25, 0.50]      # ambient UMI 占比
REPS = 3
N_SYN = 1500

MK = {"Excitatory": ["Slc17a7","Rbfox3","Neurod6","Tbr1","Camk2a"],
      "Inhibitory": ["Gad1","Gad2","Sst","Pvalb","Npy"],
      "Astrocyte": ["Gfap","Aqp4","Slc1a3","Aldh1l1"],
      "Microglia": ["C1qb","Tyrobp","Cx3cr1","Hexb","Aif1"],
      "Oligodendrocyte": ["Mbp","Plp1","Mog","Mag"],
      "OPC": ["Pdgfra","Olig1","Cspg4"],
      "Endothelial": ["Pecam1","Cldn5","Flt1"],
      "Pericyte": ["Pdgfrb","Rgs5","Vtn"]}
IND_NEU = ["Snap25","Stmn2","Syt1"]
IND_MG  = ["Trem2","Cd68","Lyz2"]

# ---------- 池 ----------
X_pool = sp.load_npz(os.path.join(SIM_DIR, "X_pool.npz"))
info = np.load(os.path.join(SIM_DIR, "pool_info.npz"), allow_pickle=True)
pool_genes = list(map(str, info["pool_genes"]))
pool_mt = info["pool_marker_type"]
gpos = {g: i for i, g in enumerate(pool_genes)}
EXC = np.where(pool_mt == "Excitatory")[0]
MIC = np.where(pool_mt == "Microglia")[0]
log(f"池: Excitatory {len(EXC)}, Microglia {len(MIC)}")

def idx_of(genes): return np.array([gpos[g] for g in genes if g in gpos])

# ---------- Part 1: ambient 谱 ----------
def build_ambient(marker_type_list, label, n_genes):
    """从官方标签细胞构建平均 counts 谱（31,017 维 counts 向量）"""
    ann = pd.read_csv(ANNOT_CSV)
    ann["bc_plain"] = ann["barcode"].str.replace(r"_GSM\d+$", "", regex=True)
    sel = ann[ann["cell_type"].astype(str).isin(marker_type_list)]
    log(f"构建 {label} ambient 谱: {len(sel)} 细胞")
    feat_path = sorted(glob.glob(os.path.join(RAW_DIR, "*_features.tsv.gz")))[0]
    with gzip.open(feat_path, "rt") as f:
        feats = [l.strip().split("\t") for l in f if not l.startswith("#")]
    genes_full = [x[1] if len(x) > 1 else x[0] for x in feats]
    # 只加载这些细胞所在 GSM 的矩阵太大——改：全 340K 无法一次加载
    # 方案：逐 GSM 提取目标细胞行（与 survey 脚本一致），全基因列太宽 → 仅保留池内表达的 31,017 基因? 
    # 简化：直接逐 GSM 加载全矩阵取目标细胞（每 GSM 一次 mmread），列保持全 31,017
    gene_sel = np.array([genes_full.index(g) if g in genes_full else -1 for g in pool_genes])
    assert (gene_sel >= 0).all(), "池基因应全在 features 中"
    cnt = np.zeros(n_genes)
    tot = 0
    for gsm in sorted(sel["gsm"].unique()):
        grp = sel[sel["gsm"] == gsm]
        mtx_f = glob.glob(os.path.join(RAW_DIR, f"{gsm}_*_matrix.mtx.gz"))
        bar_f = glob.glob(os.path.join(RAW_DIR, f"{gsm}_*_barcodes.tsv.gz"))
        with gzip.open(mtx_f[0], "rb") as f: mat = scipy.io.mmread(f)
        mat = mat.tocsr()                      # genes × cells
        with gzip.open(bar_f[0], "rt") as f: barcodes = [l.strip() for l in f]
        bc2row = {b: i for i, b in enumerate(barcodes)}
        rows = [bc2row[b] for b in grp["bc_plain"] if b in bc2row]
        if rows:
            sub = mat[:, rows].tocsr()         # genes × n_sel
            sub = sub[gene_sel, :]             # 31017 × n_sel
            cnt += np.asarray(sub.sum(axis=1)).ravel()
            tot += len(rows)
        del mat
    cnt = cnt / max(tot, 1)
    np.save(os.path.join(OUT_DIR, f"ambient_{label}.npy"), cnt)
    log(f"ambient_{label}: 平均 UMI {cnt.sum():.0f}, top基因: {[pool_genes[i] for i in np.argsort(-cnt)[:5]]}")
    return cnt

# 少突主导 ambient（模拟 GSE330130 脊髓 snRNA 场景）—— 需加载官方 Oligo 细胞
A_oligo = build_ambient(["Oligodendrocyte"], "oligo", X_pool.shape[1])
# 微胶质主导 ambient（池内即可）
Xm = X_pool[MIC].tocsr()
A_mg = np.asarray(Xm.mean(axis=0)).ravel()
np.save(os.path.join(OUT_DIR, "ambient_mg.npy"), A_mg)
log(f"ambient_mg: 平均 UMI {A_mg.sum():.0f}, top: {[pool_genes[i] for i in np.argsort(-A_mg)[:5]]}")

# ---------- Part 2: 网格仿真 ----------
Xp_cpm_raw = X_pool.tocsr()          # counts
mk_idx = {ct: idx_of(gs) for ct, gs in MK.items()}
ind_neu_i = idx_of(IND_NEU); ind_mg_i = idx_of(IND_MG)
all_i = np.concatenate(list(mk_idx.values()))
is_neu_panel = np.zeros(X_pool.shape[1], dtype=bool)
is_neu_panel[np.concatenate([mk_idx["Excitatory"], mk_idx["Inhibitory"]])] = True

def norm_log1p_cpm(X):
    u = np.asarray(X.sum(axis=1)).ravel()
    Xn = X.multiply(1e6 / np.maximum(u, 1)[:, None]).tocsr()
    Xn.data = np.log1p(Xn.data)
    return Xn

def score_cols(Xn, cols):
    return np.asarray(Xn[:, cols].mean(axis=1)).ravel() if len(cols) else np.zeros(Xn.shape[0])

def c1_nonneuron(Xn):
    """32-panel argmax 8类, 阈值>0; 返回非 Excitatory/Inhibitory 比例"""
    scs = {ct: score_cols(Xn, mk_idx[ct]) for ct in MK}
    sdf = pd.DataFrame(scs); best = sdf.idxmax(axis=1).values; bv = sdf.max(axis=1).values
    t = np.where(bv > 0, best, "Unknown")
    return 1 - np.isin(t, ["Excitatory", "Inhibitory"]).mean()

def discriminate(Xn):
    """判别协议：对每个细胞算本型(Exc) vs 异型(Mg)独立基因分；返回 (污染指示, 判别分)"""
    s_n = score_cols(Xn, ind_neu_i)
    s_m = score_cols(Xn, ind_mg_i)
    return s_n, s_m

rows = []
for amb_name, A in [("none", None), ("microglia", A_mg), ("oligo", A_oligo)]:
    for c in CLEVELS:
        for d in DLEVELS:
            if amb_name == "none" and d > 0:
                continue
            est_c1, est_disc, est_mix = [], [], []
            for rep in range(REPS):
                rng = np.random.RandomState(int(c * 1000 + d * 100 + rep))
                n_exc = int(N_SYN * (1 - c)); n_mic = N_SYN - n_exc
                ex_i = rng.choice(EXC, n_exc, replace=True)
                mi_i = rng.choice(MIC, n_mic, replace=True)
                sel = np.concatenate([ex_i, mi_i])
                Xt = X_pool[sel].tocsr()          # N_SYN × 31017 counts
                if amb_name != "none" and d > 0:
                    # ambient 谱取 top-500 高丰度基因稀疏化（ambient 由高丰度基因主导，低丰度贡献可忽略）
                    order = np.argsort(-A)
                    keep = order[:500]
                    A_sp = sp.csr_matrix((A[keep], (np.zeros(500, dtype=int), keep)), shape=(1, X_pool.shape[1]))
                    A_rep = sp.vstack([A_sp] * N_SYN)
                    Xt = Xt.multiply(1 - d) + A_rep.multiply(d)
                    Xt = Xt.tocsr()
                Xn = norm_log1p_cpm(Xt)
                # C1 原始
                e1 = c1_nonneuron(Xn)
                # 判别: 对"argmax 非神经"的细胞，用独立基因判本型(Exc)与否
                scs = {ct: score_cols(Xn, mk_idx[ct]) for ct in MK}
                sdf = pd.DataFrame(scs); best = sdf.idxmax(axis=1).values; bv = sdf.max(axis=1).values
                t = np.where(bv > 0, best, "Unknown")
                nonneu = ~np.isin(t, ["Excitatory", "Inhibitory"])
                s_n = score_cols(Xn, ind_neu_i); s_m = score_cols(Xn, ind_mg_i)
                # 真污染判定: 异型(Mg)独立基因 ≥ τ·本型(Exc)独立基因, τ=1.0（对非神经候选）
                conf = nonneu & (s_m > s_n)
                est_c1.append(e1); est_disc.append(conf.mean())
            rows.append({"ambient": amb_name, "true_c": c, "d": d,
                         "C1_raw": round(float(np.mean(est_c1) * 100), 1),
                         "disc_corrected": round(float(np.mean(est_disc) * 100), 1)})
            log(f"amb={amb_name:8s} c={c:4.0%} d={d:4.0%} → C1_raw={rows[-1]['C1_raw']}% disc={rows[-1]['disc_corrected']}%")
res = pd.DataFrame(rows)
res.to_csv(os.path.join(OUT_DIR, "ambient_simulation.csv"), index=False)
print("\n===== 结果表 =====")
print(res.to_string(index=False))
log("DONE")
