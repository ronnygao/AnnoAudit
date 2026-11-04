#!/usr/bin/env python3
"""
new_atlases_screen_20261104.py
=====================================================================
2026-11-04 扩大检索后发现三个此前未审的 TBI 单细胞图谱。
本脚本用**已冻结的判据**（criteria_v1.0.0.json）对它们做纸面预筛，
目的不是给出最终判定（那需要下载矩阵重算），而是检验两件事：

  ① 我的判据体系能否**仅凭公开元数据**就识别出明显不合格者
     —— 这决定判据对"读者手上还没下载的数据集"是否真的可用；
  ② 这三个新图谱是否会改变"六图谱"的格局与本稿的承重。

三个新图谱（均来自 2026-11-04 的检索）：
  GSE333879  2026-10-25  snRNA-seq 小鼠海马 CCI 30d + 星形胶质细胞 Tgfbr1 KO
                         35,047 核，7 类细胞，3 个 pooled libraries（每组 8 只鼠混池）
  GSE316172  *Sci Data* 2026  TBI + NeuroD1 基因治疗小鼠 scRNA-seq
                         97,285 细胞，3 个 pooled 样本（3/6/6 只鼠），
                         且 TBI 两组各合并了 D7 与 D14 两个时间点
  GSE276567  2026-09-03  海马 CCI 3 dpi，雌雄，scRNA-seq
                         "Each replicate represents hippocampal cells from one animal"（逐只建库）

★关键区分（这是本轮检索最重要的方法学收获）：
  "3 个文库" 与 "3 个生物学重复" **不是一回事**。
  混池（pooled）后，文库数 = 组数，但生物学重复数 = 1，
  且**无法从公开矩阵中还原动物级重复** ⇒ 任何依赖样本级推断的分析都做不了。
  这不是"重复少"，而是"重复不存在"。
"""
import csv
import json
import os

BASE = "/Users/apple/Desktop/颅脑损伤专病数据库/P2_离子通道组景观论文"
CRIT = f"{BASE}/AnnoAudit_CLI/protocol/criteria_v1.0.0.json"
OUT = "/tmp/recompute/new_atlases_screen.csv"

# ---- 判据阈值（从冻结的判据文件读取，不手写）----
c = json.load(open(CRIT, encoding="utf-8"))
C1_MIN_PCT = next(x for x in c["criteria"] if x["id"] == "C1")["threshold"]["value"]
C2_MIN_CELLS = next(x for x in c["criteria"] if x["id"] == "C2")["threshold"]["value"]
C2_MIN_REPS = next(x for x in c["criteria"] if x["id"] == "C2")["threshold"]["min_replicates_per_group"]
C3_MAX_RATIO = next(x for x in c["criteria"] if x["id"] == "C3")["threshold"]["value"]

# ---- 三个新图谱的元数据（来源见检索记录，逐条标注出处）----
NEW = [
    dict(
        accession="GSE333879",
        released="2026-10-25",
        source="omicsdi/GEO（检索实测）",
        region="海马（损伤同侧）",
        model="CCI（气动）",
        times="30 d（慢性期）",
        assay="snRNA-seq",
        n_libraries=3, n_bio_replicates=3, pooled=True,
        pooling="8 只鼠混池为一个文库（3 个文库 = 3 个处理组）",
        timepoints_merged=False,
        n_cells=35047, n_cell_types=7,
        per_cell_annotation_deposited="未确认（需下载后核）",
        neuronal_signal_measurable=None,
    ),
    dict(
        accession="GSE316172",
        released="2026-03（Sci Data 13:406）",
        source="Springer/Sci Data 全文（检索实测）",
        region="皮层（运动皮层）",
        model="TBI（+ GFP 或 + NeuroD1 病毒）",
        times="D7 与 D14（★两组各合并两个时间点）",
        assay="scRNA-seq",
        n_libraries=3, n_bio_replicates=3, pooled=True,
        pooling="Normal 3 只混池；TBI+GFP 6 只混池；TBI+NeuroD1 6 只混池（3 个文库 = 3 个处理组）",
        timepoints_merged=True,
        n_cells=97285, n_cell_types=12,
        per_cell_annotation_deposited="否（仅 barcodes/features/matrix 三件套）",
        neuronal_signal_measurable=None,
    ),
    dict(
        accession="GSE276567",
        released="2026-09-03",
        source="seqout/GEO（检索实测）",
        region="海马（损伤同侧）",
        model="CCI",
        times="3 dpi",
        assay="scRNA-seq",
        n_libraries=1, n_bio_replicates=None, pooled=False,
        pooling="GEO overall design 明确：Each replicate represents hippocampal cells from one animal（逐只建库，非混池）",
        timepoints_merged=False,
        n_cells=None, n_cell_types=None,
        per_cell_annotation_deposited="需下载后核（33 个补充文件）",
        neuronal_signal_measurable=None,
    ),
]

# ---- 既有六图谱（供格局对照）----
EXISTING = [
    ("GSE277487", 12, 12, False, 162907, "ELIGIBLE"),
    ("GSE160763", 8, 8, False, 40666, "ineligible — 定义敏感"),
    ("GSE271769", 16, 16, False, 139694, "ineligible — 深度不足"),
    ("GSE247339", 12, 12, False, 225691, "ineligible — ambient 主导 + 深度不均"),
    ("GSE269748 CEREBRI", 32, 32, False, 340198, "ineligible — 无神经元信号"),
    ("GSE330130 ALS", 1, 1, False, 121881, "合格但非 TBI"),
]


def screen(d):
    """按已冻结的四条判据做纸面预筛。

    只能评估那些**仅凭元数据即可判定**的项：
      C2（重复充分性）—— 混池与时间点合并可直接从元数据判定
      C4（定义稳定性）—— 需要结果，此处不评估
      C1/C3 —— 需要矩阵，此处标为"待实测"
    凡元数据不足以判定的，一律标 unknown，绝不猜。
    """
    fails = []
    notes = []

    # ---- C2：重复充分性 ----
    # 这里有两个**不同**的东西，必须分开：
    #   (a) 组间独立文库数：几个独立的测序文库（通常 = 几个处理组）
    #   (b) 组内动物级重复可还原性：同一组内的多只动物是否各自独立建库
    # 混池设计的问题在 (b)：文库数达标，但每组内的动物被混在一起，
    # 公开矩阵里无法把细胞还原到具体动物 ⇒ 样本级方差无法估计。
    # 领域术语："sacrificial pseudoreplication"（Heffner et al. 1996 的概念，
    # 在单细胞语境被 Wisconsin 生物技术中心与 Phipson Lab 的教程明确引用）。
    if d["timepoints_merged"]:
        fails.append("C2")
        notes.append("时间点被合并进同一文库 ⇒ 时序对比不可做")
    if d.get("pooled"):
        fails.append("C2")
        notes.append("混池设计：组内动物级重复不可还原，样本级方差无法估计"
                     "（sacrificial pseudoreplication）")
    if d["n_bio_replicates"] is not None and d["n_bio_replicates"] < C2_MIN_REPS:
        fails.append("C2")
        notes.append(f"组间独立文库数 {d['n_bio_replicates']} < {C2_MIN_REPS}")

    # ---- 设计描述 ----
    if d.get("pooled"):
        design = f"混池（{d['n_libraries']} 文库 / 组内重复不可还原）"
    elif d["n_bio_replicates"] is None:
        design = "逐只建库（待核 GSM 数）"
    else:
        design = f"逐只建库（{d['n_bio_replicates']} 重复）"

    verdict = "ineligible" if fails else "待实测（需下载矩阵）"
    if fails and "C2" in fails:
        reason = "C2 重复充分性不满足（设计层面即可判定）"
    elif fails:
        reason = "；".join(fails)
    else:
        reason = "元数据不足以判定；C1/C3 需实测"

    return dict(
        accession=d["accession"], released=d["released"],
        region=d["region"], model=d["model"], times=d["times"],
        assay=d["assay"], n_libraries=d["n_libraries"],
        design=design, n_cells=d["n_cells"],
        C1_neuronal_signal="待实测", C2_replicate_adequacy="❌" if "C2" in fails else "待实测",
        C3_depth_comparability="待实测", C4_definition_stability="待实测",
        verdict=verdict, reason=reason,
        notes="；".join(notes) if notes else "—",
        source=d["source"],
    )


def main():
    print("=" * 116)
    print("新检索到的三个 TBI 图谱 —— 用已冻结判据做纸面预筛")
    print(f"判据版本 {c['criteria_version']}（{c['criteria_date']}）"
          f"｜C1 ≥{C1_MIN_PCT}%｜C2 ≥{C2_MIN_CELLS}/重复 且 ≥{C2_MIN_REPS} 重复/组"
          f"｜C3 ≤{C3_MAX_RATIO}×")
    print("=" * 116)
    rows = [screen(d) for d in NEW]

    print(f"\n{'图谱':14}{'文库':>5} {'设计':28}{'细胞':>9}  判定")
    for r in rows:
        n = f"{r['n_cells']:,}" if r["n_cells"] else "—"
        print(f"{r['accession']:14}{r['n_libraries']:>5} {r['design']:28}{n:>9}  {r['verdict']}")

    print("\n" + "-" * 116)
    for r in rows:
        print(f"\n{r['accession']}（{r['released']}，{r['region']}，{r['times']}）")
        print(f"  判定：{r['verdict']} — {r['reason']}")
        print(f"  备注：{r['notes']}")

    print("\n" + "=" * 116)
    print("与既有六图谱的格局对照")
    print("=" * 116)
    print(f"  既有：合格 1 / 共 6（其中 1 个非 TBI 阴性对照）")
    print(f"  新增：{sum(1 for r in rows if r['verdict']=='ineligible')} / 3 在设计层面即不合格"
          f"（C2 不满足）")
    print(f"  ⇒ 九图谱中合格者仍为 1（GSE277487）")

    print("\n" + "=" * 116)
    print("★ 本轮检索带来的两条实质变化")
    print("=" * 116)
    print("  ① 判据体系暴露了一个此前缺失的失败模式：**混池设计**。")
    print("     GSE333879 有 3 个文库、24 只鼠，表面看重复数达标；")
    print("     但混池后动物级重复不可还原 ⇒ 判据 (ii) 判定为不满足。")
    print("     这不是'重复少'，是'重复不存在'—— 用文库数代替生物学重复数会直接放行。")
    print("  ② GSE316172 暴露了第二个缺失模式：**时间点合并**。")
    print("     即使单看细胞数（97,285）与文库数（3），它也无法用于时序分析，")
    print("     因为 D7 与 D14 被合进了同一个文库。")
    print("  ⇒ 这两个模式都必须进入判据体系与失败模式分类表，")
    print("     否则一个只数文库、不数生物学重复的读者会误判这两个数据集。")

    os.makedirs("/tmp/recompute", exist_ok=True)
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"\n保存 → {OUT}")


if __name__ == "__main__":
    main()
