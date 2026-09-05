#!/usr/bin/env python3
"""v17b：正式 D1 规则的 ambient 仿真（替换 v17 初版 disc）
D1：对 argmax=异型 M 的细胞，比较 M 的独立基因 vs 本型 T 的独立基因（τ=1.5）
少突 ambient → 兴奋性细胞 argmax=Oligo → 比较 MBP组 vs SNAP25组（SNAP25 高 → 不判污染）
"""
import os, glob, gzip, time, warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import scipy.io, scipy.sparse as sp

T0 = time.time()
def log(m): print(f"[{time.time()-T0:7.1f}s] {m}", flush=True)
SIM_DIR = "/Users/apple/Desktop/颅脑损伤专病数据库/CEREBRI_KCNC3生信分析/v12_simulation"
OUT_DIR = "/Users/apple/Desktop/颅脑损伤专病数据库/P2_离子通道组景观论文/AnnoAudit_survey/v17b_ambient_sim"
os.makedirs(OUT_DIR, exist_ok=True)
RAW_DIR = "/Users/apple/Workbuddy/2026-07-22-16-57-18/cerebri_data/raw"
ANNOT_CSV = "/Users/apple/Desktop/颅脑损伤专病数据库/CEREBRI_KCNC3生信分析/v3分析结果/cell_annotations.csv"
CLEVELS = [0.0, 0.05, 0.10, 0.20, 0.50, 0.80]
DLEVELS = [0.0, 0.10, 0.25, 0.50]
REPS = 3
N_SYN = 1500
TAU = 1.5

MK = {"Excitatory": ["Slc17a7","Rbfox3","Neurod6","Tbr1","Camk2a"],
      "Inhibitory": ["Gad1","Gad2","Sst","Pvalb","Npy"],
      "Astrocyte": ["Gfap","Aqp4","Slc1a3","Aldh1l1"],
      "Microglia": ["C1qb","Tyrobp","Cx3cr1","Hexb","Aif1"],
      "Oligodendrocyte": ["Mbp","Plp1","Mog","Mag"],
      "OPC": ["Pdgfra","Olig1","Cspg4"],
      "Endothelial": ["Pecam1","Cldn5","Flt1"],
      "Pericyte": ["Pdgfrb","Rgs5","Vtn"]}
MK8 = list(MK.keys())
IND = {"Excitatory": ["Snap25","Stmn2","Syt1"], "Inhibitory": ["Snap25","Stmn2","Syt1"],
       "Microglia": ["Trem2","Cd68","Lyz2"], "Astrocyte": ["Gja1","Slc1a2","S100b"],
       "Oligodendrocyte": ["Mbp","Plp1","Mog","Mag"], "OPC": ["Pdgfra","Olig1","Cspg4"],
       "Endothelial": ["Pecam1","Cldn5","Flt1"], "Pericyte": ["Pdgfrb","Rgs5","Vtn"]}

X_pool = sp.load_npz(os.path.join(SIM_DIR, "X_pool.npz"))
info = np.load(os.path.join(SIM_DIR, "pool_info.npz"), allow_pickle=True)
pool_genes = list(map(str, info["pool_genes"]))
pool_mt = info["pool_marker_type"]
gpos = {g: i for i, g in enumerate(pool_genes)}
EXC = np.where(pool_mt == "Excitatory")[0]
MIC = np.where(pool_mt == "Microglia")[0]
def idx(genes): return np.array([gpos[g] for g in genes if g in gpos])
mk_idx = {ct: idx(gs) for ct, gs in MK.items()}
ind_idx = {ct: idx(gs) for ct, gs in IND.items()}

# ambient 谱（复用 v17 输出）
A_mg = np.load(os.path.join(os.path.dirname(OUT_DIR), "v17_ambient_sim/ambient_mg.npy"))
A_oligo = np.load(os.path.join(os.path.dirname(OUT_DIR), "v17_ambient_sim/ambient_oligo.npy"))

def norm(X):
    u = np.asarray(X.sum(axis=1)).ravel()
    Xn = X.multiply(1e6 / np.maximum(u, 1)[:, None]).tocsr()
    Xn.data = np.log1p(Xn.data)
    return Xn

def score_cols(Xn, cols):
    return np.asarray(Xn[:, cols].mean(axis=1)).ravel() if len(cols) else np.zeros(Xn.shape[0])

def audit(Xn, known_exc):
    """返回 (marker_only_contam, D1_contam)——本型 Excitatory"""
    scs = {ct: score_cols(Xn, mk_idx[ct]) for ct in MK8}
    sdf = pd.DataFrame(scs); amt = sdf.idxmax(axis=1).values; amv = sdf.max(axis=1).values
    amt = np.where(amv > 0, amt, "Unknown")
    marker_only = 1 - np.isin(amt, ["Excitatory", "Inhibitory"]).mean()
    # D1: argmax 异型 M ≠ Exc/Inh → 比较 ind[M] vs ind[Exc]（本型）
    s_self = score_cols(Xn, ind_idx["Excitatory"])
    conf = np.zeros(Xn.shape[0], dtype=bool)
    neu_ok = np.isin(amt, ["Excitatory", "Inhibitory"])
    for o in MK8:
        if o in ("Excitatory", "Inhibitory"): continue
        m_o = (amt == o)
        if m_o.any():
            s_o = score_cols(Xn, ind_idx[o])
            conf |= m_o & (s_o > TAU * s_self)
    return marker_only, conf.mean()

rows = []
for amb_name, A in [("none", None), ("microglia", A_mg), ("oligo", A_oligo)]:
    for c in CLEVELS:
        for d in DLEVELS:
            if amb_name == "none" and d > 0: continue
            m1s, d1s = [], []
            for rep in range(REPS):
                rng = np.random.RandomState(int(c * 1000 + d * 100 + rep))
                n_exc = int(N_SYN * (1 - c)); n_mic = N_SYN - n_exc
                sel = np.concatenate([rng.choice(EXC, n_exc, replace=True), rng.choice(MIC, n_mic, replace=True)])
                Xt = X_pool[sel].tocsr()
                if amb_name != "none" and d > 0:
                    order = np.argsort(-A); keep = order[:500]
                    A_sp = sp.csr_matrix((A[keep], (np.zeros(500, dtype=int), keep)), shape=(1, X_pool.shape[1]))
                    Xt = Xt.multiply(1 - d) + sp.vstack([A_sp] * N_SYN).multiply(d)
                    Xt = Xt.tocsr()
                mo, d1 = audit(norm(Xt), None)
                m1s.append(mo); d1s.append(d1)
            rows.append({"ambient": amb_name, "true_c": c, "d": d,
                         "C1_raw": round(float(np.mean(m1s) * 100), 1),
                         "D1_corrected": round(float(np.mean(d1s) * 100), 1)})
            log(f"amb={amb_name:8s} c={c:4.0%} d={d:4.0%} → C1={rows[-1]['C1_raw']:5.1f}% D1={rows[-1]['D1_corrected']:5.1f}%")
res = pd.DataFrame(rows)
res.to_csv(os.path.join(OUT_DIR, "ambient_simulation_D1.csv"), index=False)
print("\n===== D1 正式规则 ambient 仿真 =====")
print(res.to_string(index=False))
log("DONE")
