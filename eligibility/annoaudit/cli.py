"""
annoaudit.cli
=====================================================================
命令行入口。

用法：
  # 判定一个数据集（每细胞一行的 CSV/TSV）
  annoaudit eligibility --table cells.csv --name GSE277487 --out out/

  # 判定多个数据集
  annoaudit eligibility --table a.csv --name A --table b.csv --name B --out out/

  # 查看当前判据版本、哈希与全文
  annoaudit criteria show
  annoaudit criteria show --json

  # 判据组合稳健性（哪些判据在承重）
  annoaudit robustness --table a.csv --name A --out out/

  # 自测：用带指纹的黄金结果回归验证引擎行为
  annoaudit selftest

设计取舍：
  - 不依赖 scanpy / anndata：只读扁平表（每细胞一行）。
    理由：读者手上通常已有 matrix/annotation，但未必有 h5ad；
    扁平表是最大的公约数，且可用任何工具导出。
  - 判据来自 protocol/criteria_vX.Y.Z.json，不写死在代码里（见 criteria.py 注释）。
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from typing import List

from . import __version__
from .criteria import available_versions, criteria_hash, load_criteria, spec
from .eligibility import (combination_robustness, evaluate, read_table)


def _outdir(p: str) -> str:
    os.makedirs(p, exist_ok=True)
    return p


# ====================================================================
# eligibility
# ====================================================================

def cmd_eligibility(args) -> int:
    outdir = _outdir(args.out)
    sp = spec(args.criteria_version)
    results = []
    for path, name in zip(args.table, args.name):
        rows, fields = read_table(path)
        r = evaluate(name, rows, fields, sp=sp)
        results.append(r)

    cols = ["dataset", "n_cells", "n_replicates", "strict_per_rep_min",
            "strict_per_rep_max", "strict_total", "neuronal_ge5_pct",
            "neuronal_ge10_pct", "neuronal_ge20_pct", "umi_ratio_max",
            "crit_C1_neuronal_signal", "crit_C2_replicate_adequacy",
            "crit_C3_depth_comparability", "crit_C4_definition_stability",
            "eligible", "failure_mode", "failure_signature",
            "recoverable", "remedy"]
    vp = os.path.join(outdir, "eligibility_verdict.csv")
    with open(vp, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in results:
            w.writerow(r.as_row())

    print(f"判据版本 {sp.version}（{sp.date}，sha256:{sp.hash[:16]}）")
    print("=" * 92)
    print(f"{'数据集':22}{'重复':>5}{'真群/重复':>13}{'≥10占比':>10}{'UMI比':>8}  判定")
    for r in results:
        per = (f"{r.strict_per_rep_min}~{r.strict_per_rep_max}"
               if r.strict_per_rep_min or r.strict_per_rep_max else "—")
        v = "ELIGIBLE" if r.eligible else ("—" if r.eligible is None else "ineligible")
        print(f"{r.name:22}{r.n_replicates:>5}{per:>13}"
              f"{r.neuronal_ge10_pct:>10.4f}"
              f"{('' if r.umi_ratio_max is None else f'{r.umi_ratio_max:.2f}'):>8}  {v}")
    print(f"\n→ {vp}")

    # 人类可读报告
    rp = os.path.join(outdir, "eligibility_report.md")
    with open(rp, "w", encoding="utf-8") as fh:
        fh.write(f"# AnnoAudit eligibility report\n\n")
        fh.write(f"Criteria version **{sp.version}** ({sp.date}), "
                 f"sha256 `{sp.hash}`.\n\n")
        for r in results:
            fh.write(f"## {r.name}\n\n")
            fh.write(f"- Cells: {r.n_cells:,}\n")
            fh.write(f"- Biological replicates: {r.n_replicates}\n")
            fh.write(f"- Marker-confirmed excitatory neurons per replicate: "
                     f"{r.strict_per_rep_min}–{r.strict_per_rep_max} "
                     f"(total {r.strict_total:,})\n")
            fh.write(f"- Neuronal signal (>=10 counts, strict neuronal gene): "
                     f"{r.neuronal_ge10_pct:.4f}%\n")
            if r.umi_ratio_max is not None:
                fh.write(f"- Max/min group median UMI ratio: {r.umi_ratio_max:.3f}\n")
            fh.write(f"- **Verdict: "
                     f"{'ELIGIBLE (confirmatory)' if r.eligible else 'ineligible'}"
                     f"**\n")
            if not r.eligible:
                fh.write(f"- Failure mode: **{r.failure_mode}**\n")
                fh.write(f"- Diagnostic signature: {r.failure_signature}\n")
                fh.write(f"- Recoverable in principle: {r.recoverable}\n")
                fh.write(f"- Remedy: {r.remedy}\n")
            fh.write("\n")
    print(f"→ {rp}")
    return 0


# ====================================================================
# robustness
# ====================================================================

def cmd_robustness(args) -> int:
    outdir = _outdir(args.out)
    sp = spec(args.criteria_version)
    results = []
    for path, name in zip(args.table, args.name):
        rows, fields = read_table(path)
        results.append(evaluate(name, rows, fields, sp=sp))
    rowsout = combination_robustness(results, sp)
    vp = os.path.join(outdir, "criteria_combination_robustness.csv")
    with open(vp, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rowsout[0].keys()))
        w.writeheader()
        for r in rowsout:
            w.writerow(r)
    print("判据组合稳健性（合格集合是否随判据选择/阈值改变）")
    print("=" * 92)
    print(f"  {'组合':34}{'≥5':>10}{'≥10':>10}{'≥20':>10}   稳定")
    for r in rowsout:
        def brief(s):
            return s.replace(";", ",") if s else "—"
        print(f"  {r['combination']:34}{brief(r['at_5']):>10}{brief(r['at_10']):>10}"
              f"{brief(r['at_20']):>10}   {'是' if r['stable_across_thresholds'] else '否'}")
    print(f"\n→ {vp}")
    return 0


# ====================================================================
# criteria
# ====================================================================

def cmd_criteria(args) -> int:
    if args.action == "list":
        for v in available_versions():
            print(f"{v}  sha256:{criteria_hash(v)}")
        return 0
    crit = load_criteria(args.criteria_version)
    if args.json:
        print(json.dumps(crit, ensure_ascii=False, indent=2))
        return 0
    sp = spec(args.criteria_version)
    print(f"AnnoAudit criteria v{sp.version}  ({sp.date})  sha256:{sp.hash}")
    print("=" * 92)
    for c in sp.criteria:
        print(f"\n{c.id}  {c.name}   [role: {c.role}]")
        print(f"   {c.statement}")
        print(f"   理由：{c.rationale}")
    print("\n组合规则：", sp.raw["combination"]["statement"])
    print("\n定义族：")
    d = sp.definitions
    print("  事前固定的嵌套族：", ", ".join(d.nested.keys()))
    print("  事后加入的独立族：", ", ".join(d.independent.keys()))
    return 0


# ====================================================================
# selftest
# ====================================================================

def cmd_selftest(args) -> int:
    from .selftest import run_selftest
    return run_selftest(args)


# ====================================================================
# parser
# ====================================================================

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="annoaudit",
        description="AnnoAudit: confirmatory-eligibility criteria for neuron-level "
                    "inference in public single-cell atlases.",
    )
    p.add_argument("--version", action="version", version=f"annoaudit {__version__}")
    p.add_argument("--criteria-version", default="1.0.0",
                   help="判据文件版本（默认 1.0.0）")
    sub = p.add_subparsers(dest="cmd", required=True)

    e = sub.add_parser("eligibility", help="对一个或多个数据集执行资格判定")
    e.add_argument("--table", action="append", required=True,
                   help="每细胞一行的 CSV/TSV（可重复）")
    e.add_argument("--name", action="append", required=True,
                   help="与 --table 一一对应的数据集名（可重复）")
    e.add_argument("--out", default="annoaudit_out", help="输出目录")
    e.set_defaults(func=cmd_eligibility)

    r = sub.add_parser("robustness", help="检验合格集合对判据组合与阈值的依赖")
    r.add_argument("--table", action="append", required=True)
    r.add_argument("--name", action="append", required=True)
    r.add_argument("--out", default="annoaudit_out")
    r.set_defaults(func=cmd_robustness)

    c = sub.add_parser("criteria", help="查看判据全文")
    c.add_argument("action", choices=["show", "list"], nargs="?", default="show")
    c.add_argument("--json", action="store_true")
    c.set_defaults(func=cmd_criteria)

    s = sub.add_parser("selftest", help="用带指纹的黄金结果回归验证引擎")
    s.add_argument("--verbose", action="store_true")
    s.set_defaults(func=cmd_selftest)

    return p


def main(argv: List[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
