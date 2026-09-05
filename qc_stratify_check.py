#!/usr/bin/env python3
"""QC 分层验证：GSE330130 官方神经元的 ambient 比例是否随检测基因数下降"""
import anndata as ad
import numpy as np, pandas as pd

H5 = "/Users/apple/Desktop/颅脑损伤专病数据库/P2_离子通道组景观论文/AnnoAudit_survey/data/GSE330130_LSC_SNRNASEQ_COUNTS.h5ad"
a = ad.read_h5ad(H5)
a.var["symbol"] = a.var["common_name"].astype(str)
a.var_names = a.var["symbol"].values
a = a[:, ~pd.Index(a.var_names).duplicated(keep="first")].copy()
X = a.X.tocsr().astype(np.float32)
umi = np.asarray(X.sum(axis=1)).ravel()
Xn = X.multiply(1e6 / np.maximum(umi, 1)[:, None]).tocsr()
Xn.data = np.log1p(Xn.data)
VN = a.var_names

def sm(genes):
    cols = [VN.get_loc(g) for g in genes if g in VN]
    return np.asarray(Xn[:, cols].mean(axis=1)).ravel()

s_neur = sm(["SNAP25", "STMN2", "SYT1"])
s_olig = sm(["MBP", "PLP1", "MOG", "MAG"])
nz = np.diff(X.indptr)  # csr: 行 = 细胞
neu = a.obs["Cell_Type"].astype(str).isin(["Inhibitory Neuron", "Excitatory  Neuron", "Motor Neurons"]).values
is_oligo = (a.obs["Cell_Type"].astype(str) == "Oligodendrocytes").values

df = pd.DataFrame({"nz": nz, "s_neur": s_neur, "s_olig": s_olig, "neu": neu, "oligo": is_oligo})
sub = df[df["neu"]].copy()
sub["q"] = pd.qcut(sub["nz"], 4, labels=["Q1_low", "Q2", "Q3", "Q4_high"])
print("官方神经元按检测基因数分层 (总 n=%d):" % len(sub))
print(f"{'layer':<10}{'n':>6}{'nz_med':>8}{'Neuron_like%':>13}{'Oligo_like%':>13}")
for q in ["Q1_low", "Q2", "Q3", "Q4_high"]:
    g = sub[sub["q"] == q]
    nl = (g["s_neur"] > g["s_olig"]).mean() * 100
    ol = (g["s_olig"] > g["s_neur"]).mean() * 100
    print(f"{q:<10}{len(g):>6}{g['nz'].median():>8.0f}{nl:>12.1f}{ol:>13.1f}")

print("\n参考: 真少突标签 s_neur 中位 = %.3f, s_olig 中位 = %.3f" % (
    s_neur[df["oligo"].values].median(), s_olig[df["oligo"].values].median()))
print("官方神经元 s_neur 中位 = %.3f (全体)" % s_neur[neu].median())
# 按层输出 C1-like 估计: 少突样占比 = 判别 "Oligo_like"
