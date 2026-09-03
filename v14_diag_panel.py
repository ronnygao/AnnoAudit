#!/usr/bin/env python3
# GSE330130 C1 panel 敏感性: 32 基因直系同源大 panel (±SLC17A6) vs 20 基因小 panel
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, anndata as ad

H5 = "/Users/apple/Desktop/TBI分子时钟/数据/GSE330130_LSC_SNRNASEQ_COUNTS.h5ad"
NEU = ["Excitatory", "Inhibitory"]
NEURON_LABELS = ["Inhibitory Neuron", "Excitatory  Neuron", "Motor Neurons"]

PANEL_LARGE = {  # v7 人类大 panel（脊髓含 SLC17A6）= 33 基因
    "Excitatory": ["SLC17A7", "SLC17A6", "RBFOX3", "NEUROD6", "TBR1", "CAMK2A"],
    "Inhibitory": ["GAD1", "GAD2", "SST", "PVALB", "NPY"],
    "Astrocyte": ["GFAP", "AQP4", "SLC1A3", "ALDH1L1"],
    "Microglia": ["C1QB", "TYROBP", "CX3CR1", "HEXB", "AIF1"],
    "Oligodendrocyte": ["MBP", "PLP1", "MOG", "MAG"],
    "OPC": ["PDGFRA", "OLIG1", "CSPG4"],
    "Endothelial": ["PECAM1", "CLDN5", "FLT1"],
    "Pericyte": ["PDGFRB", "RGS5", "VTN"],
}
PANEL_ORTHO = {k: ([g for g in v if g != "SLC17A6"] if k == "Excitatory" else v)
               for k, v in PANEL_LARGE.items()}  # 32 基因 = 小鼠大 panel 严格直系同源

a = ad.read_h5ad(H5)
a.var["symbol"] = a.var["common_name"].astype(str)
a = a[:, ~pd.Index(a.var["symbol"]).duplicated(keep="first")].copy()
a.var_names = a.var["symbol"].values
X = a.X.tocsr()
umi = np.asarray(X.sum(axis=1)).ravel()
Xn = X.multiply(1e6 / np.maximum(umi, 1)[:, None]).tocsr()
Xn.data = np.log1p(Xn.data)
off = a.obs["Cell_Type"].isin(NEURON_LABELS).values
Xo = Xn[off].tocsr()

def c1(panel, thresh):
    s = {}
    for ct, gs in panel.items():
        idx = [a.var_names.get_loc(g) for g in gs if g in a.var_names]
        miss = [g for g in gs if g not in a.var_names]
        if miss: print(f"  缺失 {ct}: {miss}")
        s[ct] = np.asarray(Xo[:, idx].mean(axis=1)).ravel()
    sdf = pd.DataFrame(s)
    bv = sdf.max(axis=1).values
    best = sdf.idxmax(axis=1)
    mt = np.where(bv > thresh, best.values, "Unknown")
    vc = pd.Series(mt).value_counts()
    c1v = 1 - pd.Series(mt).isin(NEU).mean()
    return c1v, vc

for name, panel in [("33基因大panel(含SLC17A6)", PANEL_LARGE), ("32基因直系同源", PANEL_ORTHO)]:
    for th in [0.0, 0.1]:
        c1v, vc = c1(panel, th)
        print(f"{name} th={th}: C1={c1v*100:5.1f}%  Exc={vc.get('Excitatory',0):4d} Inh={vc.get('Inhibitory',0):4d} "
              f"Oligo={vc.get('Oligodendrocyte',0):4d} Astro={vc.get('Astrocyte',0):4d} Unk={vc.get('Unknown',0):3d}")
print("\nv7 原始(X原样,th=0):     C1= 56.9%  Exc= 411 Inh= 645 Oligo=1123 Astro= 197 Unk=  8")
print("v14 小panel(th=0.1):     C1= 79.9%  Exc=  91 Inh= 401 Oligo=1602 Astro= 193 Unk= 24")
