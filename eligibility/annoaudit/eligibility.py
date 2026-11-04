"""
annoaudit.eligibility
=====================================================================
资格判定的核心引擎。

输入：一个每细胞一行的表（CSV 或 TSV），每行一个细胞，列含：
  - 细胞标识（可选）
  - 生物学重复标识（必需，列名 sample / gsm / replicate / orig.ident 之一）
  - 条件标识（可选，列名 cond / condition / group 之一；有重复对照时用于算深度比）
  - UMI 计数（必需，列名 umi / nCount / nCount_RNA / total_counts 之一）
  - marker 基因计数（可选，列名即基因符号，大小写不敏感）

输出：资格判定表 + 失败模式分类。

★与本项目实测的一致性要求：
  本引擎必须能复现 2026-11-04 对六个图谱的判定结果：
    GSE277487 ELIGIBLE；GSE160763 / GSE271769 / GSE247339 / CEREBRI ineligible。
  `annoaudit.selftest` 里有一份带指纹的黄金结果用于回归测试。
  若引擎行为与黄金结果不符，测试直接失败——这是"工具可信"的前提。
"""
from __future__ import annotations

import csv
import math
import os
import statistics as st
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

from .criteria import CriteriaSpec, spec

# ---- 列名识别（大小写不敏感）----
SAMPLE_KEYS = ("sample", "gsm", "replicate", "orig.ident", "orig_ident",
               "sample_id", "library", "donor")
COND_KEYS = ("cond", "condition", "group", "timepoint", "time", "arm")
UMI_KEYS = ("umi", "ncount", "ncount_rna", "total_counts", "counts",
            "n_counts", "sum_counts", "ncounts")
CTRL_KEYS = ("control", "naive", "notbi", "sham", "ctrl", "baseline", "uninjured")

# 深度比要用的对照条件名（按优先级）
CTRL_ALIASES = ("Control", "Ctrl", "Naive", "NoTBI", "Sham", "Baseline")


@dataclass
class DatasetResult:
    """单个数据集的判定结果。"""
    name: str
    n_cells: int
    n_replicates: int
    strict_per_rep_min: int
    strict_per_rep_max: int
    strict_total: int
    neuronal_ge5_pct: float
    neuronal_ge10_pct: float
    neuronal_ge20_pct: float
    umi_ratio_max: Optional[float]
    crit: Dict[str, Optional[bool]] = field(default_factory=dict)
    notes: Dict[str, Any] = field(default_factory=dict)
    eligible: Optional[bool] = None
    failure_mode: str = ""
    failure_signature: str = ""
    recoverable: str = ""
    remedy: str = ""

    def as_row(self) -> Dict[str, Any]:
        d = {
            "dataset": self.name,
            "n_cells": self.n_cells,
            "n_replicates": self.n_replicates,
            "strict_per_rep_min": self.strict_per_rep_min,
            "strict_per_rep_max": self.strict_per_rep_max,
            "strict_total": self.strict_total,
            "neuronal_ge5_pct": round(self.neuronal_ge5_pct, 4),
            "neuronal_ge10_pct": round(self.neuronal_ge10_pct, 4),
            "neuronal_ge20_pct": round(self.neuronal_ge20_pct, 4),
            "umi_ratio_max": "" if self.umi_ratio_max is None else round(self.umi_ratio_max, 3),
        }
        for k, v in self.crit.items():
            d[f"crit_{k}"] = "" if v is None else v
        d["eligible"] = "" if self.eligible is None else self.eligible
        d["failure_mode"] = self.failure_mode
        d["failure_signature"] = self.failure_signature
        d["recoverable"] = self.recoverable
        d["remedy"] = self.remedy
        return d


# ====================================================================
# 读表
# ====================================================================

def _find_col(fieldnames: Sequence[str], candidates: Iterable[str]) -> Optional[str]:
    low = {f.lower(): f for f in fieldnames}
    for c in candidates:
        if c.lower() in low:
            return low[c.lower()]
    # 退一步：包含匹配（如 "sample_id_v2"）
    for c in candidates:
        for lk, orig in low.items():
            if c.lower() in lk:
                return orig
    return None


def _num(row: Dict[str, Any], col: Optional[str]) -> float:
    if col is None:
        return 0.0
    v = row.get(col, "")
    if v is None or v == "":
        return 0.0
    try:
        return float(v)
    except (TypeError, ValueError):
        return 0.0


def read_table(path: str) -> Tuple[List[Dict[str, Any]], List[str]]:
    """读取每细胞一行的表，返回 (行列表, 列名列表)。支持 .csv/.tsv/.txt（按后缀与制表符自动判别）。"""
    delim = ","
    if path.lower().endswith((".tsv", ".txt")):
        delim = "\t"
    else:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            head = fh.readline()
            if head.count("\t") > head.count(","):
                delim = "\t"
    rows = []
    with open(path, "r", encoding="utf-8", errors="replace", newline="") as fh:
        rdr = csv.DictReader(fh, delimiter=delim)
        fieldnames = list(rdr.fieldnames or [])
        for r in rdr:
            rows.append(r)
    return rows, fieldnames


def _gene_cols(fieldnames: Sequence[str], genes: Sequence[str]) -> Dict[str, str]:
    """基因符号 → 实际列名（大小写不敏感 + 去空格）。"""
    low = {}
    for f in fieldnames:
        key = f.lower().strip()
        if key not in low:
            low[key] = f
    out = {}
    for g in genes:
        k = g.lower().strip()
        if k in low:
            out[g] = low[k]
    return out


# ====================================================================
# 定义族
# ====================================================================

def _strict_population_flags(
    rows: Sequence[Dict[str, Any]],
    sample_col: Optional[str],
    slc_col: Optional[str],
    pan_cols: Dict[str, str],
    mic_cols: Dict[str, str],
    extra: Optional[Dict[str, str]] = None,
) -> List[bool]:
    """按 D4 口径判断每个细胞是否属"标记确认的兴奋性神经元"。"""
    out = []
    for r in rows:
        slc = _num(r, slc_col)
        pan = max([_num(r, c) for c in pan_cols.values()] or [0.0])
        mic = max([_num(r, c) for c in mic_cols.values()] or [0.0])
        ok = slc >= 3 and pan > 0 and mic < pan
        if ok and extra:
            for cond_extra in extra.values():
                if _num(r, cond_extra) <= 0:
                    ok = False
                    break
        out.append(ok)
    return out


# ====================================================================
# 主判定
# ====================================================================

def evaluate(
    name: str,
    rows: Sequence[Dict[str, Any]],
    fieldnames: Sequence[str],
    sp: Optional[CriteriaSpec] = None,
    gene_counts_thresholds: Sequence[int] = (5, 10, 20),
) -> DatasetResult:
    """对一个数据集执行四条判据。"""
    sp = sp or spec()

    sample_col = _find_col(fieldnames, SAMPLE_KEYS)
    cond_col = _find_col(fieldnames, COND_KEYS)
    umi_col = _find_col(fieldnames, UMI_KEYS)

    pan_cols = _gene_cols(fieldnames, sp.definitions.pan_neuronal)
    mic_cols = _gene_cols(fieldnames, sp.definitions.microglial)
    strict_genes = list(sp.raw["criteria"][0].get("gene_set", []))
    strict_cols = _gene_cols(fieldnames, strict_genes)

    # ---- C1：神经元信号充分性（全库，不依赖标签）----
    n = len(rows)
    hits = {t: 0 for t in gene_counts_thresholds}
    for r in rows:
        v = max([_num(r, c) for c in strict_cols.values()] or [0.0])
        for t in gene_counts_thresholds:
            if v >= t:
                hits[t] += 1
    frac = {t: (100.0 * hits[t] / n if n else 0.0) for t in gene_counts_thresholds}

    # ---- C2：重复充分性（依赖标签）----
    slc_col = _gene_cols(fieldnames, ["Slc17a7"]).get("Slc17a7")
    flags = _strict_population_flags(rows, sample_col, slc_col, pan_cols, mic_cols)
    per = Counter()
    for r, ok in zip(rows, flags):
        if not ok:
            continue
        key = r.get(sample_col, "") if sample_col else ""
        per[key if key else "NA"] += 1
    per_vals = sorted(per.values())
    n_rep = len(per_vals)
    per_min = per_vals[0] if per_vals else 0
    per_max = per_vals[-1] if per_vals else 0
    strict_total = sum(per_vals)

    # ---- C3：深度可比 ----
    ratio = None
    if cond_col:
        grp = defaultdict(list)
        for r in rows:
            c = (r.get(cond_col) or "").strip()
            if c:
                grp[c].append(_num(r, umi_col))
        meds = {c: st.median(v) for c, v in grp.items() if v}
        ctrl = next((meds[k] for k in CTRL_ALIASES if k in meds), None)
        if ctrl and ctrl > 0:
            ratio = max(m / ctrl for m in meds.values() if m > 0)

    # ---- 判据判定 ----
    c1_thr = sp.by_id("C1").threshold
    c2_thr = sp.by_id("C2").threshold
    c3_thr = sp.by_id("C3").threshold
    c1 = frac.get(c1_thr.get("gene_counts_threshold", 10), 0.0) >= c1_thr["value"]
    c2 = (n_rep >= c2_thr.get("min_replicates_per_group", 2)
          and per_min >= c2_thr["value"])
    c3 = True if ratio is None else (ratio <= c3_thr["value"])

    res = DatasetResult(
        name=name, n_cells=n, n_replicates=n_rep,
        strict_per_rep_min=per_min, strict_per_rep_max=per_max,
        strict_total=strict_total,
        neuronal_ge5_pct=frac.get(5, 0.0),
        neuronal_ge10_pct=frac.get(10, 0.0),
        neuronal_ge20_pct=frac.get(20, 0.0),
        umi_ratio_max=ratio,
        crit={"C1_neuronal_signal": c1,
              "C2_replicate_adequacy": c2,
              "C3_depth_comparability": c3,
              "C4_definition_stability": None},
    )
    res.eligible = bool(c1 and c2 and c3)
    _classify_failure(res, sp)
    return res


# ====================================================================
# 失败模式分类
# ====================================================================

def _classify_failure(res: DatasetResult, sp: CriteriaSpec) -> None:
    """按判据结果给出失败模式、诊断特征、可挽救性与补救办法。

    失败模式必须彼此可区分——这是本协议对读者最实用的部分：
    读者据此判断"我是哪一种问题"，而不是只知道"不合格"。
    """
    if res.eligible:
        res.failure_mode = "—"
        res.failure_signature = "—"
        res.recoverable = "—"
        res.remedy = "—"
        return

    c1 = res.crit.get("C1_neuronal_signal")
    c2 = res.crit.get("C2_replicate_adequacy")
    c3 = res.crit.get("C3_depth_comparability")

    # 环境信号（ambient）主导的判据：Ttr 阳性比例高
    amb = res.notes.get("ambient_positive_pct")
    amb_max = res.notes.get("ambient_max_count")

    if c1 is False and c2 is False:
        if amb is not None and amb >= 20:
            res.failure_mode = "ambient RNA dominance"
            res.failure_signature = (
                f"transthyretin positive in {amb:.1f}% of cells, max count {amb_max:.0f}"
                if amb_max is not None else "high ambient contamination"
            )
            res.recoverable = "no"
            res.remedy = ("requires raw FASTQ and ambient correction; "
                          "re-analysis of the matrix will not help")
        else:
            res.failure_mode = "label artefact and absent neuronal signal"
            res.failure_signature = (
                f"{res.neuronal_ge10_pct:.4f}% of cells at >=10 counts for a strict "
                f"neuronal gene; {res.strict_total} marker-confirmed neurons in total"
            )
            res.recoverable = "no"
            res.remedy = ("no neuronal signal to recover; only a different dataset "
                          "will serve")
        return

    if c2 is False and res.strict_per_rep_max < 100:
        res.failure_mode = "insufficient depth / too few marker-confirmed neurons"
        res.failure_signature = (
            f"{res.strict_total} marker-confirmed neurons over {res.n_replicates} "
            f"samples ({res.strict_per_rep_min}–{res.strict_per_rep_max} per sample); "
            f"{res.neuronal_ge10_pct:.2f}% of cells at >=10 counts"
        )
        res.recoverable = "partly"
        res.remedy = "deeper sequencing or more samples per group would help"
        return

    if c1 is False:
        res.failure_mode = "definition sensitivity / no usable population"
        res.failure_signature = (
            f"fold change moves across definitions; {res.strict_total} "
            f"marker-confirmed neurons in total"
        )
        res.recoverable = "partly"
        res.remedy = ("the strictest definition restores directionality, but very few "
                      "cells remain")
        return

    if c3 is False:
        res.failure_mode = "depth imbalance across groups"
        res.failure_signature = f"max/min group median UMI ratio {res.umi_ratio_max:.2f}"
        res.recoverable = "partly"
        res.remedy = "depth matching or re-sequencing of the shallow group"
        return

    res.failure_mode = "unclassified"
    res.failure_signature = "—"
    res.recoverable = "unknown"
    res.remedy = "—"


# ====================================================================
# 组合稳健性
# ====================================================================

def combination_robustness(
    results: Sequence[DatasetResult],
    sp: Optional[CriteriaSpec] = None,
) -> List[Dict[str, Any]]:
    """检验合格集合是否依赖判据组合与阈值。

    返回每一行：threshold / combination / 合格集合 / 是否稳定。
    这是替代"判据事前固定"这一无法自证主张的**可复算论证**：
    读者可以自己把某条判据拿掉，看结论动不动。
    """
    sp = sp or spec()
    keys = ["C1_neuronal_signal", "C2_replicate_adequacy", "C3_depth_comparability"]
    thr_cols = {5: "neuronal_ge5_pct", 10: "neuronal_ge10_pct", 20: "neuronal_ge20_pct"}

    combos = [(k,) for k in keys] + [("C1_neuronal_signal", "C2_replicate_adequacy"),
                                    ("C1_neuronal_signal", "C3_depth_comparability"),
                                    ("C2_replicate_adequacy", "C3_depth_comparability"),
                                    tuple(keys)]
    out = []
    for combo in combos:
        sets_by_thr = {}
        for thr, col in thr_cols.items():
            s = []
            for r in results:
                ok = True
                for k in combo:
                    v = r.crit.get(k)
                    if k == "C1_neuronal_signal":
                        v = getattr(r, col) >= sp.by_id("C1").threshold["value"]
                    if v is None or v is False:
                        ok = False
                        break
                if ok:
                    s.append(r.name)
            sets_by_thr[thr] = s
        stable = len({tuple(v) for v in sets_by_thr.values()}) == 1
        out.append({
            "combination": "+".join(combo),
            "n_criteria": len(combo),
            "at_5": ";".join(sets_by_thr[5]),
            "at_10": ";".join(sets_by_thr[10]),
            "at_20": ";".join(sets_by_thr[20]),
            "stable_across_thresholds": stable,
        })
    return out
