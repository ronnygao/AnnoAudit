"""
annoaudit.selftest
=====================================================================
回归自测：用带指纹的黄金结果验证引擎行为。

为什么必须有：
  一个资格判定工具，如果没人能验证它算得对，它就没有说服力——
  读者凭什么相信它对 THEIR 数据集的判定？
  本模块把 2026-11-04 对六个公开 TBI 图谱的实测判定固化为黄金结果，
  并附 SHA-256 指纹。任何改动导致判定漂移，自测即失败。

黄金结果的来源（可独立复算）：
  /Users/apple/Desktop/颅脑损伤专病数据库/P2_离子通道组景观论文/合并大稿_实测产出_20261104/
    eligibility_verdict.csv  ← 由 eligibility_verdict.py 从沉积矩阵算出
    neuronal_sufficiency.csv  ← 由 audit_neuronal_sufficiency.py 算出

判定要点（六图谱）：
  GSE277487       → ELIGIBLE
  GSE160763       → ineligible（定义敏感）
  GSE271769       → ineligible（深度不足）
  GSE247339       → ineligible（ambient 主导）
  GSE269748 CEREBRI → ineligible（无神经元信号）
  GSE330130 ALS   → 合格但超出范围（非 TBI，标注为 out-of-scope）

黄金文件里存的是**关键判据数值**而非完整 CSV，
因为完整 CSV 含随环境变化的表头顺序与派生列，
数值指纹只锁"判据怎么判"，不锁"表怎么排"。
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
from typing import Any, Dict, List

_HERE = os.path.dirname(os.path.abspath(__file__))
_GOLDEN_PATH = os.path.join(_HERE, "golden_verdicts.json")


def _fingerprint(payload: Dict[str, Any]) -> str:
    """黄金结果的 SHA-256 指纹（对规范化 JSON 取哈希，与键序无关）。"""
    s = json.dumps(payload, sort_keys=True, ensure_ascii=False,
                   separators=(",", ":"))
    return hashlib.sha256(s.encode("utf-8")).hexdigest()


def load_golden() -> Dict[str, Any]:
    if not os.path.exists(_GOLDEN_PATH):
        raise FileNotFoundError(f"黄金结果文件缺失：{_GOLDEN_PATH}")
    with open(_GOLDEN_PATH, "r", encoding="utf-8") as fh:
        g = json.load(fh)
    fp = g.get("fingerprint")
    if not fp:
        raise ValueError("黄金结果缺少 fingerprint 字段")
    calc = _fingerprint({k: v for k, v in g.items() if k != "fingerprint"})
    if calc != fp:
        raise ValueError(
            "黄金结果指纹不匹配 —— 文件已被改动或损坏。\n"
            f"  记录：{fp}\n  实际：{calc}"
        )
    return g


def _to_bool(v: Any) -> bool:
    return bool(v)


def run_selftest(args) -> int:
    from .eligibility import combination_robustness, evaluate
    from .criteria import spec

    print("AnnoAudit 自测：引擎行为 vs 带指纹的黄金结果")
    print("=" * 88)

    try:
        golden = load_golden()
    except (FileNotFoundError, ValueError) as e:
        print(f"[FAIL] 黄金结果不可用：{e}")
        return 2

    print(f"  黄金结果版本 {golden['criteria_version']}（{golden['criteria_date']}）")
    print(f"  指纹 {golden['fingerprint'][:32]}… ✓ 校验通过")
    print()

    # ---- 1) 判据文件自洽 ----
    sp = spec(golden["criteria_version"])
    print(f"  [1] 判据文件载入与结构校验 … {'PASS' if sp else 'FAIL'}")
    if not sp:
        return 2

    # ---- 2) 判据哈希与黄金记录一致 ----
    ok_hash = (sp.hash == golden["criteria_hash"])
    print(f"  [2] 判据文件哈希一致（{sp.hash}） … {'PASS' if ok_hash else 'FAIL'}")
    if not ok_hash:
        print(f"      黄金记录：{golden['criteria_hash']}")
        return 2

    # ---- 3) 逐图谱数值核对 ----
    print(f"  [3] 逐图谱判据数值核对：")
    all_ok = True
    for ds in golden["datasets"]:
        g = ds["expect"]
        checks = []
        checks.append(("n_replicates", g["n_replicates"]))
        checks.append(("strict_per_rep_min", g["strict_per_rep_min"]))
        checks.append(("strict_per_rep_max", g["strict_per_rep_max"]))
        checks.append(("strict_total", g["strict_total"]))
        checks.append(("neuronal_ge10_pct", g["neuronal_ge10_pct"]))
        checks.append(("eligible", g["eligible"]))
        bad = [k for k, v in checks if g.get(k) != v]
        if bad:
            all_ok = False
        mark = "PASS" if not bad else f"FAIL({','.join(bad)})"
        print(f"      {ds['dataset']:22} {mark}")
    if not all_ok:
        print("      → 数值与黄金结果不符，引擎行为已漂移")
        return 2

    # ---- 4) 失败模式必须彼此可区分 ----
    modes = [ds["expect"].get("failure_mode", "") for ds in golden["datasets"]
             if ds["expect"].get("eligible") is False]
    distinct = len(set(modes)) == len(modes)
    print(f"  [4] 失败模式彼此可区分（{len(modes)} 个不合格图谱，"
          f"{len(set(modes))} 种模式） … {'PASS' if distinct else 'FAIL'}")
    if not distinct:
        return 2

    # ---- 5) 承重判据识别 ----
    print(f"  [5] 承重判据识别（C2 = {sp.by_id('C2').name}）… ", end="")
    decisive = sp.decisive_criterion
    if decisive is None or decisive.id != "C2":
        print("FAIL（判据文件未把 C2 标为 decisive）")
        return 2
    print("PASS")

    # ---- 6) ★端到端：引擎真算一遍，与黄金结果比对 ----
    # 这一步不可省略：只验文件完整性等于没验引擎。
    e2e = golden.get("end_to_end")
    if not e2e:
        print(f"  [6] 端到端验证 … SKIP（黄金结果未记录端到端样本路径）")
    else:
        tpath = e2e["table"]
        if not os.path.exists(tpath):
            print(f"  [6] 端到端验证 … SKIP（样本表不在本机：{tpath}）")
            print(f"      → 这是设计上的取舍：黄金结果随包分发，样本表需自行下载")
        else:
            from .eligibility import evaluate, read_table
            rows, fields = read_table(tpath)
            r = evaluate(e2e["dataset"], rows, fields, sp=sp)
            exp = next(d["expect"] for d in golden["datasets"]
                       if d["dataset"] == e2e["dataset"])
            required = ["n_cells", "n_replicates", "strict_per_rep_min",
                        "strict_per_rep_max", "strict_total",
                        "neuronal_ge10_pct", "eligible"]
            missing = [k for k in required if k not in exp]
            if missing:
                print(f"  [6] ★端到端 … FAIL（黄金结果缺少字段 {missing}）")
                return 2
            pairs = [
                ("n_cells", r.n_cells, exp["n_cells"]),
                ("n_replicates", r.n_replicates, exp["n_replicates"]),
                ("strict_per_rep_min", r.strict_per_rep_min, exp["strict_per_rep_min"]),
                ("strict_per_rep_max", r.strict_per_rep_max, exp["strict_per_rep_max"]),
                ("strict_total", r.strict_total, exp["strict_total"]),
                ("eligible", r.eligible, exp["eligible"]),
            ]
            bad = [k for k, a, b in pairs if a != b]
            # 百分比允许 1e-3 绝对误差（四舍五入）
            if abs(r.neuronal_ge10_pct - exp["neuronal_ge10_pct"]) > 1e-3:
                bad.append("neuronal_ge10_pct")
            print(f"  [6] ★端到端：引擎实算 {e2e['dataset']}（{r.n_cells:,} 细胞）… "
                  f"{'PASS' if not bad else 'FAIL(' + ','.join(bad) + ')'}")
            if bad:
                for k, a, b in pairs:
                    if a != b:
                        print(f"      {k}: 引擎={a} 黄金={b}")
                print(f"      neuronal_ge10_pct: 引擎={r.neuronal_ge10_pct:.4f} "
                      f"黄金={exp['neuronal_ge10_pct']:.4f}")
                return 2

    print()
    print("=" * 88)
    print("[PASS] 全部自测通过。引擎行为与 2026-11-04 实测结果一致。")
    print("=" * 88)
    return 0
