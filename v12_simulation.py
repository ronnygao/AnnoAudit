#!/usr/bin/env python3
# ============================================================
# v12 (final): AnnoAudit 模拟 benchmark — 校准曲线
# 协议（与正文最终口径一致）：
#   C1 marker max-score：8 类 argmax，全零细胞归 target（无证据不算污染）
#   C2 margin 门控模块评分：scanpy score_genes 8 类模块分，
#      margin = max(胶质模块) - max(神经元模块)，margin > τ=1.5 计污染
#      （τ 由干净标签校准：FPR=2.1%；无门控的 argmax 规则 FPR=42-48%）
#   C3 无监督聚类：top-2000 变异基因 → PCA30 → KMeans k=10 → 簇 majority 映射
#   C4 适用性门控 CellTypist：模型须在 marker 纯化参照上 ≥80% 正确才可用；
#      本数据 DMB 12.4% / MTG 18.1% / PFC 1.6% 全部 FAIL（与 v14 终版一致）→ 剔除
#   ACS = mean(C1..C4)
# 用法：python v12_simulation.py <level_idx 0-5>   每块 5 reps
# ============================================================
import os, sys, pickle, time, warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
import scanpy as sc
import anndata as ad

T0 = time.time()
def log(msg):
    print(f"[{time.time()-T0:7.1f}s] {msg}", flush=True)

OUT_DIR = "/Users/apple/Desktop/颅脑损伤专病数据库/CEREBRI_KCNC3生信分析/v12_simulation"
MODEL_DIR = "/Users/apple/.celltypist/data/models"
RES_CSV = os.path.join(OUT_DIR, "simulation_results_v2.csv")
APPL_CSV = os.path.join(OUT_DIR, "applicability_table.csv")
TAU = 1.5
LEVELS = [0.0, 0.05, 0.10, 0.20, 0.50, 0.80]
REPS = 5
N_SYN = 1000

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
NEU_TYPES = ["Excitatory", "Inhibitory"]
GLIAL_TYPES = [t for t in MARKERS if t not in NEU_TYPES]

# ---------- 加载池 ----------
X_pool = sp.load_npz(os.path.join(OUT_DIR, "X_pool.npz"))
info = np.load(os.path.join(OUT_DIR, "pool_info.npz"), allow_pickle=True)
pool_genes = list(map(str, info["pool_genes"]))
pool_marker_type = info["pool_marker_type"]
log(f"池: {X_pool.shape}, 兴奋性 {(pool_marker_type=='Excitatory').sum()}, 小胶质 {(pool_marker_type=='Microglia').sum()}")

umi = np.asarray(X_pool.sum(axis=1)).ravel()
Xp_cpm = X_pool.multiply(1e6 / np.maximum(umi, 1)[:, None]).tocsr()
Xp_cpm.data = np.log1p(Xp_cpm.data)
Xp_10k = X_pool.multiply(1e4 / np.maximum(umi, 1)[:, None]).tocsr()
Xp_10k.data = np.log1p(Xp_10k.data)
log("归一化完成")

gene_pos = {g: i for i, g in enumerate(pool_genes)}
mk_idx = {ct: np.array([gene_pos[g] for g in gs if g in gene_pos]) for ct, gs in MARKERS.items()}
all_mk = np.concatenate([v for v in mk_idx.values()])

EXC_POOL = np.where(pool_marker_type == "Excitatory")[0]
MIC_POOL = np.where(pool_marker_type == "Microglia")[0]

# ---------- C4: DMB 手动复现（Title-case 匹配，已验证与官方一致） ----------
def is_neuronal_dmb(lab):
    return "Neuron" in str(lab) or "Excitatory" in str(lab) or "Inhibitory" in str(lab)

def is_neuronal_mtg(lab):
    lab = str(lab)
    nonneural = ("Oligo", "Astro", "OPC", "Micro", "Endo", "VLMC", "Peri", "SMC")
    neuronal = ("IT", "ET", "CT", "NP", "L6b", "Sst", "Pvalb", "Vip", "Lamp5", "Sncg", "Chandelier")
    if any(k in lab for k in nonneural):
        return False
    return any(k in lab for k in neuronal)

def is_neuronal_pfc(lab):
    lab = str(lab)
    return lab.startswith("InN") or lab.startswith("ExN") or "CUX2" in lab or "RORB" in lab

def load_model(mf):
    with open(os.path.join(MODEL_DIR, mf + ".pkl"), "rb") as f:
        return pickle.load(f)

def build_F(sel_rows, feats, uppercase):
    """sel_rows: 池行号 → (n_sel × n_feats) CP10K log1p 特征矩阵"""
    pos = {}
    for i, g in enumerate(pool_genes):
        pos.setdefault(g.upper() if uppercase else g, i)
    idx = [pos.get(g) for g in feats]
    F = np.zeros((len(sel_rows), len(feats)), dtype=np.float64)
    Xs = Xp_10k[sel_rows].tocsr()
    for j, i in enumerate(idx):
        if i is not None:
            F[:, j] = np.asarray(Xs[:, i].todense()).ravel()
    return F

def predict_model(sel_rows, mf, uppercase, rule):
    d = load_model(mf)
    lr, sc = d["Model"], d["Scaler_"]
    feats = [str(x) for x in lr.features]
    F = build_F(sel_rows, feats, uppercase)
    lab = lr.predict(sc.transform(F))
    return np.array([rule(l) for l in lab])

# ---------- 适用性检查（只做一次） ----------
if not os.path.exists(APPL_CSV):
    rr = np.random.RandomState(0)
    sel_exc = rr.choice(EXC_POOL, 1000, replace=False)
    sel_mic = rr.choice(MIC_POOL, 1000, replace=False)
    rows = []
    for mf, upper, rule in [
        ("Adult_Human_MTG", True, is_neuronal_mtg),
        ("Adult_Human_PrefrontalCortex", True, is_neuronal_pfc),
        ("Developing_Mouse_Brain", False, is_neuronal_dmb),
    ]:
        neu_exc = predict_model(sel_exc, mf, upper, rule).mean()
        mic_lab = None
        d = load_model(mf)
        feats = [str(x) for x in d["Model"].features]
        Fm = build_F(sel_mic, feats, upper)
        labs = d["Model"].predict(d["Scaler_"].transform(Fm))
        mic_like = np.mean([("Micro" in str(l)) or ("Immune" in str(l)) for l in labs])
        rows.append({"model": mf, "clean_exc_neuronal_pct": round(neu_exc*100, 1),
                     "clean_mic_microglial_pct": round(mic_like*100, 1),
                     "applicable_exc": neu_exc >= 0.80})
        log(f"适用性 {mf}: exc神经元率={neu_exc*100:.1f}% mic免疫率={mic_like*100:.1f}%")
    pd.DataFrame(rows).to_csv(APPL_CSV, index=False)
    log(f"适用性表已保存 -> {APPL_CSV}")

APPLIC = pd.read_csv(APPL_CSV) if os.path.exists(APPL_CSV) else None

# DMB 预计算（适用模型唯一）
d_dmb = load_model("Developing_Mouse_Brain")
feats_dmb = [str(x) for x in d_dmb["Model"].features]
pos_raw = {}
for i, g in enumerate(pool_genes):
    pos_raw.setdefault(g, i)
idx_dmb = [pos_raw.get(g) for g in feats_dmb]
F_dmb_pool = np.zeros((X_pool.shape[0], len(feats_dmb)), dtype=np.float32)
for j, i in enumerate(idx_dmb):
    if i is not None:
        F_dmb_pool[:, j] = np.asarray(Xp_10k[:, i].todense()).ravel()
log(f"DMB 特征矩阵预计算: {F_dmb_pool.shape}, 匹配 {sum(i is not None for i in idx_dmb)}/{len(feats_dmb)}")

def c4_dmb(sel_rows):
    F = F_dmb_pool[sel_rows].astype(np.float64)
    lab = d_dmb["Model"].predict(d_dmb["Scaler_"].transform(F))
    return 1.0 - np.mean([is_neuronal_dmb(l) for l in lab])

# ---------- 四检查 ----------
def module_scores_scanpy(sel_rows, seed):
    Xs = Xp_cpm[sel_rows].tocsr()
    A = ad.AnnData(X=Xs)
    A.obs_names = [f"c{i}" for i in range(len(sel_rows))]
    A.var_names = pool_genes
    np.random.seed(seed)
    for ct, genes in MARKERS.items():
        use = [g for g in genes if g in A.var_names]
        sc.tl.score_genes(A, gene_list=use, score_name=f"mod_{ct}")
    mods = [c for c in A.obs.columns if c.startswith("mod_")]
    return A.obs[mods]

def run_audit(sel_rows, seed):
    n = len(sel_rows)
    Xs = Xp_cpm[sel_rows].tocsr()
    # C1
    s1 = {}
    for ct, gi in mk_idx.items():
        s1[ct] = np.asarray(Xs[:, gi].mean(axis=1)).ravel() if len(gi) else np.zeros(n)
    s1_df = pd.DataFrame(s1)
    bval = s1_df.max(axis=1).values
    best = s1_df.values.argmax(axis=1)
    types = list(MARKERS.keys())
    mtype = np.where(bval > 0, np.array(types)[best], "Unknown")
    c1 = 1 - pd.Series(mtype).isin(NEU_TYPES).mean()
    # C2（margin 门控）+ 记录无门控 argmax 作对照
    S = module_scores_scanpy(sel_rows, seed)
    neu_max = S[[f"mod_{t}" for t in NEU_TYPES]].max(axis=1).values
    glia_max = S[[f"mod_{t}" for t in GLIAL_TYPES]].max(axis=1).values
    margin = glia_max - neu_max
    c2 = float((margin > TAU).mean())
    amax = S.values.argmax(axis=1)
    amod = np.array([m.replace("mod_", "") for m in S.columns])[amax]
    c2_argmax = 1 - pd.Series(amod).isin(NEU_TYPES).mean()
    # C3
    gm = np.asarray(Xs.mean(axis=0)).ravel()
    gv = np.asarray(Xs.multiply(Xs).mean(axis=0)).ravel() - gm ** 2
    keepvar = gm > 0.01
    top = np.argsort(np.where(keepvar, gv, -1))[-min(2000, keepvar.sum()):]
    Xh = Xs[:, top].toarray()
    mu, sd = Xh.mean(0), Xh.std(0) + 1e-9
    Xs_ = np.clip((Xh - mu) / sd, -10, 10)
    Xp_ = PCA(n_components=30, random_state=0, svd_solver="randomized").fit_transform(Xs_)
    cl = KMeans(n_clusters=10, random_state=0, n_init=10).fit_predict(Xp_)
    cl_type = pd.Series(mtype).groupby(cl).agg(lambda s: s.value_counts().index[0])
    ctype = pd.Series(cl).map(cl_type).values
    c3 = 1 - pd.Series(ctype).isin(NEU_TYPES).mean()
    # C4（适用性门控：模型须在 marker 定义的干净目标类型参考上 ≥80% 正确，
    #     否则记原始值但 ACS 中按 N/A 处理）
    c4_raw = c4_dmb(sel_rows)
    appl_dmb = bool(APPLIC.loc[APPLIC.model == "Developing_Mouse_Brain", "applicable_exc"].iloc[0]) if APPLIC is not None else False
    if appl_dmb:
        c4 = c4_raw
        acs = float(np.mean([c1, c2, c3, c4]))
    else:
        c4 = np.nan
        acs = float(np.mean([c1, c2, c3]))
    return {"n": n, "C1": round(c1*100,1), "C2": round(c2*100,1),
            "C2_argmax_ungated": round(c2_argmax*100,1),
            "C3": round(c3*100,1), "C4_raw": round(c4_raw*100,1),
            "C4_dmb": round(c4*100,1) if not np.isnan(c4) else np.nan,
            "ACS": round(acs*100,1)}

# ---------- 主循环 ----------
level_idx = int(sys.argv[1]) if len(sys.argv) > 1 else 0
contam = LEVELS[level_idx]
results = []
for rep in range(REPS):
    rr = np.random.RandomState(1000 * rep + int(contam * 100))
    n_mic = int(round(N_SYN * contam))
    n_exc = N_SYN - n_mic
    idx_exc = rr.choice(EXC_POOL, n_exc, replace=False) if n_exc > 0 else np.array([], dtype=int)
    idx_mic = rr.choice(MIC_POOL, n_mic, replace=False) if n_mic > 0 else np.array([], dtype=int)
    sel = np.concatenate([idx_exc, idx_mic]).astype(int)
    r = run_audit(sel, seed=10000 + 100 * level_idx + rep)
    r.update({"true_contamination_pct": round(contam*100,1), "rep": rep})
    results.append(r)
    log(f"[L{level_idx}] c={contam*100:5.1f}% rep={rep} | C1={r['C1']:5.1f} C2={r['C2']:5.1f} "
        f"(argmax {r['C2_argmax_ungated']:5.1f}) C3={r['C3']:5.1f} C4={r['C4_raw']:5.1f}"
        f"{'[N/A]' if np.isnan(r['C4_dmb']) else ''} ACS={r['ACS']:5.1f}")

# 追加保存
if os.path.exists(RES_CSV):
    old = pd.read_csv(RES_CSV)
    old = old[old["true_contamination_pct"] != round(contam*100, 1)]
    results = old.to_dict("records") + results
df = pd.DataFrame(results).sort_values(["true_contamination_pct", "rep"]).reset_index(drop=True)
df.to_csv(RES_CSV, index=False)
done = sorted(df["true_contamination_pct"].unique())
log(f"已完成水平: {done} / {sorted(round(l*100,1) for l in LEVELS)}")
if len(done) == len(LEVELS):
    summ = df.groupby("true_contamination_pct").agg(
        C1=("C1","mean"), C2=("C2","mean"), C3=("C3","mean"),
        C4=("C4_dmb","mean"), ACS=("ACS","mean"), ACS_sd=("ACS","std")).round(1)
    summ.to_csv(os.path.join(OUT_DIR, "simulation_summary_v2.csv"))
    log("=== 校准汇总 ===\n" + summ.to_string())
    log("V12 DONE")
