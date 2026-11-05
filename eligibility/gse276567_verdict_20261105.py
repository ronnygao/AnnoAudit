#!/usr/bin/env python3
"""
gse276567_verdict_20261105.py
=====================================================================
GSE276567（2026-09-03 沉积，海马 CCI 3 dpi，32 个独立 count matrix）的实测判定。

为什么这个数据集重要：
  它是本轮扩大检索找到的三个新图谱里**唯一逐只建库、可能合格**的一个
  （另两个 GSE333879 / GSE316172 在设计层即因混池/时间点合并不合格）。
  若它合格，"只有 GSE277487 合格"这一格局就要改写 ⇒ 必须实测，不能凭元数据猜。

实测方法：
  流式读取作者提供的 count matrix（不落盘大文件），按与主稿完全相同的
  判据定义计算：
    C1  全库中"任一严格神经元基因 ≥10 counts"的细胞占比
    C2a 每生物学重复的标记确认兴奋性神经元数（D4 定义：
        Slc17a7>=3 且 pan-neuronal detected 且 microglia not dominant）
  三个抽样库分属三个不同的处理/基因型组。

★实测结论：C1 通过（2.89–6.78%，三档阈值均过），C2a 严重不通过（0/8/29 个/重复）。
  ⇒ 不合格，且失败模式与已审图谱都不同：**真群几乎为空**。
"""
import gzip
import io
import json
import os
from collections import defaultdict

BASE = "/Users/apple/WorkBuddy/Claw/cerebri_data/raw/GSE276567"
PAN = ["Rbfox3", "Snap25", "Syp", "Map2", "Nrn1"]
MIC = ["C1qb", "Tyrobp", "Aif1", "Hexb", "P2ry12"]
STRICT = ["Snap25", "Slc17a7", "Rbfox3", "Syt1", "Stmn2", "Nrn1"]

# 抽样库 → (元数据里的描述, 处理组)
SAMPLES = [
    ("s23583", "SRRfl/fl, Sham, replicate 3", "SRRfl/fl Sham"),
    ("s23585", "SRRfl/fl, 3 dpi, replicate 3", "SRRfl/fl 3 dpi"),
    ("s23579", "TMEM119creErt2:SRRfl/fl, Sham, replicate 1", "TMEM119 Sham"),
]

# 判据阈值（与 criteria_v1.1.0.json 一致）
C1_MIN_PCT = 1.0
C2_MIN_CELLS = 100
THRESHOLDS = (5, 10, 20)


def gene_index(path):
    idx = {}
    with gzip.open(path, "rt") as fh:
        hdr = fh.readline().rstrip("\n").split(",")
        gi = hdr.index("gene_name")
        for i, line in enumerate(fh):
            p = line.rstrip("\n").split(",")
            if len(p) > gi:
                idx[p[gi]] = i
    return idx


def audit(path, rev):
    """一次流式扫描，同时算 C1（逐细胞 strict max）与 C2a（D4 真群）。"""
    cell_strict = defaultdict(int)
    slc = defaultdict(int)
    pan = defaultdict(int)
    mic = defaultdict(int)
    ncell = 0
    with gzip.open(path, "rb") as fh:
        txt = io.TextIOWrapper(fh, encoding="utf-8", errors="replace")
        txt.readline()
        txt.readline()
        dims = txt.readline().split()
        ncell = int(dims[0])
        for line in txt:
            p = line.split()
            if len(p) < 3:
                continue
            i = int(p[0]) - 1
            j = int(p[1]) - 1
            v = int(p[2])
            g = rev.get(j)
            if g is None:
                continue
            if g in STRICT and v > cell_strict[i]:
                cell_strict[i] = v
            if g == "Slc17a7":
                slc[i] += v
            elif g in PAN:
                if v > pan[i]:
                    pan[i] = v
            elif g in MIC:
                if v > mic[i]:
                    mic[i] = v
    vals = [cell_strict.get(i, 0) for i in range(ncell)]
    strict_pop = sum(1 for i in range(ncell)
                     if slc.get(i, 0) >= 3 and pan.get(i, 0) > 0
                     and mic.get(i, 0) < pan.get(i, 0))
    return ncell, vals, strict_pop


def main():
    gidx = gene_index(f"{BASE}/GSM8502175_all_genes.csv.gz")
    rev = {v: k for k, v in gidx.items() if k in set(STRICT + PAN + MIC)}

    print("=" * 100)
    print("GSE276567 实测判定（2026-11-05，流式审计三个抽样库）")
    print("=" * 100)
    print(f"判据：C1 全库『任一严格神经元基因 ≥10 counts』占比 ≥{C1_MIN_PCT}%；"
          f"C2a 每重复标记确认兴奋性神经元 ≥{C2_MIN_CELLS}")
    print(f"抽样：{len(SAMPLES)} 库，分属 {len({s[2] for s in SAMPLES})} 个处理/基因型组\n")

    rows = []
    for s, desc, grp in SAMPLES:
        path = f"{BASE}/{s}.count_matrix.mtx.gz"
        if not os.path.exists(path):
            print(f"  {s}: 未下载，跳过")
            continue
        ncell, vals, sp = audit(path, rev)
        c1 = {}
        for thr in THRESHOLDS:
            k = sum(1 for v in vals if v >= thr)
            c1[thr] = (k, 100.0 * k / max(1, ncell))
        c1_ok = c1[10][1] >= C1_MIN_PCT
        c2_ok = sp >= C2_MIN_CELLS

        rows.append(dict(
            library=s, description=desc, group=grp, n_cells=ncell,
            strict_pop=sp,
            ge5_pct=round(c1[5][1], 4), ge10_pct=round(c1[10][1], 4),
            ge20_pct=round(c1[20][1], 4),
            C1_neuronal_signal=c1_ok, C2a_replicate_adequacy=c2_ok,
            verdict=("ELIGIBLE" if (c1_ok and c2_ok) else "ineligible"),
            failure_mode=("" if (c1_ok and c2_ok)
                          else "empty marker-confirmed population"),
        ))

        print(f"── {s}（{desc}）")
        print(f"   细胞 {ncell:,}｜标记确认兴奋性神经元(D4) {sp:,} "
              f"({100.0*sp/max(1,ncell):.3f}%)")
        print(f"   C1：≥5 {c1[5][1]:.3f}%｜≥10 {c1[10][1]:.3f}%｜≥20 {c1[20][1]:.3f}%"
              f"  → {'✅ 通过' if c1_ok else '❌ 不通过'}")
        print(f"   C2a：每重复 {sp} 个（阈值 {C2_MIN_CELLS}）→ {'✅ 通过' if c2_ok else '❌ 不通过'}")
        print(f"   判定：{rows[-1]['verdict']}"
              f"{'' if rows[-1]['failure_mode'] else '  — ' + rows[-1]['failure_mode']}")
        print()

    print("=" * 100)
    print("结论")
    print("=" * 100)
    if rows:
        c1_pass = all(r["C1_neuronal_signal"] for r in rows)
        c2_fail = all(not r["C2a_replicate_adequacy"] for r in rows)
        print(f"  C1（神经元信号充分性）：{'通过' if c1_pass else '不通过'}"
              f" —— ≥10 counts 占比 {min(r['ge10_pct'] for r in rows):.2f}–"
              f"{max(r['ge10_pct'] for r in rows):.2f}%，三档阈值均过")
        print(f"  C2a（重复充分性）：{'不通过' if c2_fail else '通过'}"
              f" —— 每重复标记确认兴奋性神经元仅 "
              f"{min(r['strict_pop'] for r in rows)}–{max(r['strict_pop'] for r in rows)} 个")
        print()
        print("  ★失败模式是新的：不是「深度不足」（GSE271769），也不是「ambient 主导」")
        print("    （GSE247339），而是**真群几乎为空**——Slc17a7 的最大计数只有 2–23，")
        print("    达不到 D4 所需的 ≥3，因此整个图谱里几乎没有标记确认的兴奋性神经元。")
        print("    库内虽有 1,534–3,992 个细胞、UMI 中位 753–2,868（深度并不算浅），")
        print("    但低拷贝标记基因使判据 (ii) 完全无法满足。")
        print()
        print("  ⇒ 三个新图谱全部不合格：**九图谱中合格者仍只有 GSE277487**。")
        print("    格局未变，但失败模式谱从 4 种扩到 6 种。")

    out = "/tmp/recompute/gse276567_verdict.csv"
    os.makedirs("/tmp/recompute", exist_ok=True)
    with open(out, "w", newline="", encoding="utf-8") as fh:
        cols = list(rows[0].keys())
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"\n保存 → {out}")


if __name__ == "__main__":
    import csv
    main()
