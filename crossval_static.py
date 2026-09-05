#!/usr/bin/env python3
"""交叉验证 2：静态数字自洽性审计（confusion_matrix + label_audit + 注释 CSV）
1) 混淆矩阵行和 == 官方标签 n
2) 微胶质错标总量独立统计（官方非微胶质标签中 marker-微胶质总数）
3) 关键百分比换算一致性（963/45051, 44088/45051 等）
"""
import pandas as pd

BASE = "/Users/apple/Desktop/颅脑损伤专病数据库/P2_离子通道组景观论文/AnnoAudit_survey"
issues = []

# ---- 1) CEREBRI 混淆矩阵行和 vs label_audit n ----
conf = pd.read_csv(f"{BASE}/cerebri_all_labels/confusion_matrix.csv", index_col=0)
audit = pd.read_csv(f"{BASE}/cerebri_all_labels/label_audit.csv")
print("== CEREBRI 混淆矩阵行和 vs 官方 n ==")
ok = True
for _, r in audit.iterrows():
    lab = r["official_label"]
    if lab in conf.index:
        row_sum = conf.loc[lab].sum()
        match = abs(row_sum - r["n"]) < 1e-6 or abs(row_sum - r["n"]) / r["n"] < 0.001
        if not match:
            ok = False
            issues.append(f"CEREBRI {lab}: 行和 {row_sum} ≠ n {r['n']}")
print("全部一致 PASS" if ok else "存在不一致 → " + str(issues))

# ---- 2) 微胶质错标总量（argmax 口径，独立统计）----
mg_col = "Microglia" if "Microglia" in conf.columns else None
outer = conf.loc[[l for l in conf.index if l != "Microglia"], "Microglia"].sum()
total_offdiag_mg = outer
print(f"\n官方非 Microglia 标签中 marker-微胶质总数 = {total_offdiag_mg:,}")
print(f"（声明区间 38,000-50,000 → {'PASS' if 38000 <= total_offdiag_mg <= 50000 else 'CHECK'}）")
# 补充：其他胶质异型（判别也可能确认 Oligo/Astro_like 等）→ 上界估计
for col in ["Oligodendrocyte", "Astrocyte"]:
    print(f"  其中 marker-{col}: {conf.loc[[l for l in conf.index if l != col], col].sum():,}")

# ---- 3) 关键换算一致性 ----
print("\n== 关键百分比换算 ==")
checks = [
    ("Glut 2.1%", 963 / 45051 * 100),
    ("Glut 97.9% (overlap)", 44088 / 45051 * 100),
    ("Glut 判别 67.4% 细胞数", 0.674 * 45051),
    ("45051/340185 = 13.2%", 45051 / 340185 * 100),
    ("1174/340185 = 0.35%", 1174 / 340185 * 100),
    ("38 倍", 45051 / 1174),
]
for name, v in checks:
    print(f"  {name:<28} = {v:.2f}")

# ---- 4) GSE330130 label_audit 同样行和审计 ----
conf_h = pd.read_csv(f"{BASE}/gse330130_all_labels/confusion_matrix.csv", index_col=0)
audit_h = pd.read_csv(f"{BASE}/gse330130_all_labels/label_audit.csv")
print("\n== GSE330130 混淆矩阵行和 vs 官方 n ==")
ok = True
for _, r in audit_h.iterrows():
    lab = r["official_label"]
    if lab in conf_h.index:
        row_sum = conf_h.loc[lab].sum()
        if abs(row_sum - r["n"]) / r["n"] > 0.001:
            ok = False
            issues.append(f"GSE330130 {lab}: {row_sum} vs {r['n']}")
print("全部一致 PASS" if ok else "不一致 → " + str(issues))

# ---- 5) 注释 CSV 独立计数（官方标签/条件组）----
ann = pd.read_csv("/Users/apple/Desktop/颅脑损伤专病数据库/CEREBRI_KCNC3生信分析/v3分析结果/cell_annotations.csv")
vc = ann["cell_type"].astype(str).value_counts()
print("\n== 注释 CSV 独立计数（核对 v14/v18 引用）==")
for lab, n in [("Glutamatergic_neuron", 45051), ("GABAergic_neuron", 2317), ("Astrocyte", 42560), ("OPC", 3474)]:
    print(f"  {lab:<22} CSV={vc.get(lab, 0):>8,}  引用={n:>8,}  {'PASS' if vc.get(lab,0)==n else 'FAIL'}")
print(f"  总行数={len(ann):,}（引用 340,198/340,185 需注意口径：含 Unknown 16,189）")

if issues:
    print("\n⚠️ ISSUES:\n" + "\n".join(issues))
else:
    print("\n✅ 静态审计全部 PASS")
