# Cover Letter（合并大稿，B 定位）· 2026-11-05

> **用法**：本文件是母本。投 PLOS Computational Biology 用「版本 A」，
> 投 GigaScience 用「版本 B」，两版仅第 1 段（称谓与定位）与第 5 段（契合理由）不同，
> 其余段落完全一致——**不要各写一遍**，避免两版内容漂移。
>
> **写作纪律（三条，来自本项目的历史教训）**：
> 1. **凡是数字，必须与正文一致**。本母本所有数字均取自
>    `合并大稿_实测产出_20261104/` 的实测 CSV，不得凭印象改写。
> 2. **凡是"我们做过什么"，都要能核**。三次 desk reject、预印本 DOI、判据协议版本与
>    tag，全部是可核事实，主动披露优于被动被查。
> 3. **不得出现任何第三方标注的归因**。本稿的标签全部是我方自产，全文 `official` 只出现在
>    "不得写"的检查项里。CL 同样遵守。
>
> **⚠️ 投稿前必须替换的占位符**：
> `[JOURNAL_NAME]`、`[JOURNAL_ADDRESS]`、`[DATE]`、`10.64898/2026.09.03.749125`（本稿用 v3 的
> 10.64898/2026.09.03.749125）、`https://www.biorxiv.org/content/10.64898/2026.09.03.749125v3`、`10.5281/zenodo.22274752`、
> `[ZENODO_VERSION_DOI]`、`[CODE_URL]`、`[CODE_TAG]`。

---

## 共同的第 2–4 段（两版一致）

### 第 2 段：问题与本文做什么

Two questions are routinely conflated when a published single-cell atlas is re-analysed: whether
a cell-type label is faithful to the markers it was derived from, and whether the deposited
counts contain enough of the relevant cell type, at sufficient depth and replication, for a
population-level comparison to mean anything. The second question has received almost no
attention, and it is the one that determines whether a conclusion can be confirmed at all. Our
manuscript separates the two, and addresses the second with explicit, pre-stated criteria.

We applied a marker-based audit (AnnoAudit) to **nine public traumatic brain injury atlases**,
of which one — GSE277487, 12 biological replicates and 162,907 nuclei — qualified. Across the
nine, the fraction of cells reaching ten counts for a canonical neuronal gene ranged from
**0.0012%** (CEREBRI, GSE269748) to **14.68%** (GSE277487), and it did not track cell number,
sequencing depth, or label quality. In the one qualifying atlas, KCNC3 was acutely up-regulated
in marker-confirmed excitatory neurons: the detection rate rose from **28.90%** in controls to
**47.91%** at 6 h, depth matching gave a stratum-wise difference of **+23.37 percentage points**
with the same sign in **10/10** UMI deciles, the direction was concordant in **15/15**
comparisons across five pre-specified excitatory-neuron definitions and independently
corroborated under three structurally independent definitions, and the result survived all
**twelve** leave-one-replicate-out resamples.

### 第 3 段：为什么读者需要它（可复用性，而非一次性发现）

The deliverable is a decision procedure rather than a verdict about one dataset. For each
excluded atlas we report the failure mode, its diagnostic signature, whether it is recoverable
in principle, and what remedy it would require, because the diagnosis is what tells a reader
whether their own data can be salvaged. The procedure is released as a **versioned, dated,
hashed protocol file** (`criteria_v1.1.0.json`, 2026-11-04) with a command-line
implementation that reads its thresholds from that file, so the tool cannot diverge from the
protocol it claims to implement. The package carries a fingerprinted golden result and a
self-test that recomputes one atlas end-to-end from the per-cell table, so a reader can verify
the implementation before applying it to a dataset the authors have never seen. A further
command reports, for every individual criterion and every two-criterion combination, which
datasets pass at each of three neuronal-signal thresholds — the check a sceptical reader runs
when they suspect a verdict of being an artefact of the criteria chosen.

The criteria themselves are **not presented as fixed before the analysis**, because we cannot
demonstrate that from a published record and a reader is entitled to know. Instead we substitute
a check they can run: the replicate-adequacy criterion alone selects the same atlas at every
threshold examined, and every combination that includes it gives the same verdict, whereas
criterion combinations without it are not stable. We regard this as stronger than an assurance
about our own sequence of decisions, because it can be recomputed independently.

### 第 4 段：六种失败模式（可引用的实质内容）

The audit distinguishes **six ways a public atlas fails to support neuron-level inference**, of
which four are not recoverable by re-analysis: label artefact with absent neuronal signal;
ambient RNA dominance; insufficient sequencing depth; a design in which several animals are
pooled into one library, so that the between-animal variance cannot be estimated at all; time
points merged into a single library, which precludes any temporal comparison; and a nearly empty
marker-confirmed population in which the neuronal signal is present in the library but the
excitatory population cannot be defined from it. Three of the six were identified only after we
screened atlases deposited in 2026, including one that passes every design-level requirement
and nevertheless fails because `Slc17a7` reaches a maximum count of only 2–23 where the
reference definition requires three.

---

## 版本 A：PLOS Computational Biology

> **选它的理由（内部判断，不写进 CL）**：
> PLOS Comp Biol 明确接受预印本（"posting on preprint servers does not preclude
> consideration"），初审 format-free（"waives all formatting requirements until your manuscript
> has received a provisional Editorial Accept decision"），投递成本低；
> GigaScience 走 OUP/ScienceOne 全流程，投稿包需重做。GigaScience 的
> "reproducibility, usability, utility" 判据与本稿更贴，但作为第二顺位。
> **风险（必须披露）**：Bioinformatics（BIOINF-2026-3215）明文禁投本刊任何版本，
> 但该禁投只约束 Bioinformatics 本身，不影响投 PLOS。

**[DATE]**

The Editors
PLOS Computational Biology
Public Library of Science
[JOURNAL_ADDRESS]

Dear Editors,

We submit our manuscript, "AnnoAudit: a marker-based audit with confirmatory-eligibility
criteria identifies which public traumatic brain injury atlases support neuron-level inference,"
for consideration as a Research Article in PLOS Computational Biology.

〔插入共同第 2–4 段〕

**Why this journal.** PLOS Computational Biology publishes work of exceptional significance
that uses computational methods to further our understanding of living systems, and its
preprint policy is explicit that a manuscript posted to a preprint server remains eligible for
consideration. That is the situation here, and the relationship is stated in full below. The
contribution fits the journal in a specific sense: this is a computational method paper whose
output is a decision procedure, released with the software that implements it and with the
verification needed to trust it, and it addresses the reproducibility of neuron-level analysis
in public atlases rather than a single biological question.

**Disclosure of previous submissions and related work.** We disclose the following
voluntarily, so that the editors can judge the manuscript's standing and we can be judged on
the disclosure rather than on the record. An earlier version of this work, restricted to a
single atlas, was posted as a preprint (doi:10.64898/2026.09.03.749125, https://www.biorxiv.org/content/10.64898/2026.09.03.749125v3) and appeared as an
earlier version of that record; the present manuscript is the peer-reviewed extension of that
preprint, and the Methods sections that the preprint established are reproduced verbatim. The
present manuscript adds the cross-atlas eligibility analysis, the threshold and
combination-robustness analysis, the failure-mode taxonomy, three atlases deposited in 2026,
and a released implementation of the criteria. The preprint does not report a pass-or-fail
verdict for any dataset; the present manuscript reports one, and the Discussion states
explicitly which two conditions the preprint was missing that the present work satisfies. The
manuscript is not under consideration elsewhere, and no dataset or claim in it overlaps with a
manuscript currently in review. We would also note, for completeness, that earlier drafts of
this work were rejected without review by three other journals (CNS Neuroscience & Therapeutics,
Briefings in Bioinformatics, and Bioinformatics, manuscript BIOINF-2026-3215); the rejection
letters we received concerned novelty and breadth of interest rather than any defect in the
data or analysis, and the present manuscript is a deliberate change of scope — from one atlas
to a general decision procedure — in response.

**Data, code and materials.** All nine atlases are publicly available from the Gene Expression
Omnibus under the accessions listed in the Methods. The analysis code, the marker panels, the
per-cell audit results, and the versioned criteria protocol are available under an MIT license
at [CODE_URL] (tag [CODE_TAG]); the protocol file used in this manuscript is
`criteria_v1.1.0.json`, dated 2026-11-04, with its SHA-256 recorded in the accompanying
human-readable document, and the archived release carries Zenodo concept DOI
10.5281/zenodo.22274752 (this version: [ZENODO_VERSION_DOI]). The criteria are additionally
distributed as an installable command-line tool with a self-test that recomputes one audited
atlas end-to-end, so that the implementation can be verified before use.

**Author contributions, competing interests and AI disclosure.** All authors have approved
the submission and declare no competing interests. We have no personal or institutional
relationship with the authors of any atlas examined, and we did not suggest any of them as
reviewers. During the preparation of this work the authors used AI-assisted tools for language
editing and for drafting analysis code; after using these tools the authors reviewed and
edited the content as needed and take full responsibility for the content of the publication.

Thank you for your consideration.

Sincerely,

Rong Gao (高嵘), MD, PhD
Director, Department of Neurosurgery; PI, Translational Medicine Center
The First People's Hospital of Zhangjiagang (Suzhou University Affiliated Zhangjiagang Hospital)
Suzhou, Jiangsu, China
Email: gaorong1234@suda.edu.cn
Phone: +86 13812677924

---

## 版本 B：GigaScience

**[DATE]**

The Editors
GigaScience
Oxford University Press
[JOURNAL_ADDRESS]

Dear Editors,

We submit our manuscript, "AnnoAudit: a marker-based audit with confirmatory-eligibility
criteria identifies which public traumatic brain injury atlases support neuron-level inference,"
for consideration as a Research Article in GigaScience.

〔插入共同第 2–4 段〕

**Why this journal.** GigaScience states that its criteria for publication are
reproducibility, usability and utility rather than a subjective assessment of impact, and that
it publishes the research objects — data, software, workflows — that accompany an article
alongside it. This manuscript is built to be judged on exactly those terms. The reusable object
here is not only a set of criteria but a versioned, dated, hashed protocol file with a
command-line implementation that reads its thresholds from that file, a fingerprinted golden
result, and a self-test that recomputes one atlas end-to-end from the per-cell table, so that
usability and reproducibility can be checked rather than asserted. The utility is that the
procedure returns a usable answer for a dataset the authors have never seen: two of the six
failure modes we describe are detected from the sample metadata alone, before a matrix is
downloaded.

**Disclosure of previous submissions and related work.** An earlier version of this work,
restricted to a single atlas, was posted as a preprint (doi:10.64898/2026.09.03.749125, https://www.biorxiv.org/content/10.64898/2026.09.03.749125v3);
the present manuscript is the peer-reviewed extension of that preprint, and the Methods sections
the preprint established are reproduced verbatim. The present manuscript adds the cross-atlas
eligibility analysis, the threshold and combination-robustness analysis, the failure-mode
taxonomy, three atlases deposited in 2026, and a released implementation of the criteria.
The preprint does not report a pass-or-fail verdict for any dataset; the present manuscript
reports one, and the Discussion states explicitly which two conditions the preprint was missing
that the present work satisfies. The manuscript is not under consideration elsewhere. For
completeness: earlier drafts were rejected without review by three other journals (CNS
Neuroscience & Therapeutics, Briefings in Bioinformatics, and Bioinformatics, manuscript
BIOINF-2026-3215); the rejection letters concerned novelty and breadth of interest rather than
any defect in the data or analysis, and the present manuscript is a deliberate change of scope
— from one atlas to a general decision procedure — in response.

**Data, code and materials.** All nine atlases are publicly available from GEO under the
accessions listed in the Methods. Code, marker panels, per-cell audit results and the
versioned criteria protocol are available under an MIT license at [CODE_URL] (tag [CODE_TAG]);
the protocol file used here is `criteria_v1.1.0.json`, dated 2026-11-04, with its SHA-256
recorded in the accompanying human-readable document. The archived release carries Zenodo concept
DOI 10.5281/zenodo.22274752 (this version: [ZENODO_VERSION_DOI]). We would be glad to deposit the
per-cell audit outputs in GigaDB, and the criteria tool can be presented through GigaGalaxy if
the editors judge that useful.

**Author contributions, competing interests and AI disclosure.** All authors have approved the
submission and declare no competing interests. We have no personal or institutional
relationship with the authors of any atlas examined, and we did not suggest any of them as
reviewers. During the preparation of this work the authors used AI-assisted tools for language
editing and for drafting analysis code; after using these tools the authors reviewed and
edited the content as needed and take full responsibility for the content of the publication.

Thank you for your consideration.

Sincerely,

Rong Gao (高嵘), MD, PhD
Director, Department of Neurosurgery; PI, Translational Medicine Center
The First People's Hospital of Zhangjiagang (Suzhou University Affiliated Zhangjiagang Hospital)
Suzhou, Jiangsu, China
Email: gaorong1234@suda.edu.cn
Phone: +86 13812677924

---

## 附：数字核对表（投稿前逐项对）

| CL 中的数 | 来源 | 值 |
|---|---|---|
| 九图谱 | R5 判定表 | 9 |
| 合格图谱数 | `eligibility_verdict.csv` | 1（GSE277487）|
| GSE277487 重复数 / 核数 | `audit_gse277487.csv` | 12 / 162,907 |
| 神经元信号占比范围 | `neuronal_sufficiency.csv` | 0.0012% – 14.68% |
| KCNC3 检出率 | `sensitivity_definitions.csv` | 28.90% → 47.91%（6h）|
| 深度匹配中位差 | 同上 | +23.37 pp，10/10 桶同向 |
| 定义一致性 | 同上 | 15/15 |
| 独立定义 | `probe_independent_defs.csv` | N1/N2/N3，6h +24.80/+19.88/+21.07 pp |
| LOO | `probe_loo_replicates.csv` | 36/36 同向 |
| 失败模式数 | R5.1 | 6（4 种不可救）|
| 新图谱 | `gse276567_verdict.csv` | GSE333879 / GSE316172 / GSE276567 全部不合格 |
| GSE276567 Slc17a7 max | 同上 | 2–23 |
| 判据版本 / 日期 / tag | `criteria_v1.1.0.json` | 1.1.0 / 2026-11-04 / `v1.1.0-criteria` |
| Zenodo 概念 DOI | `10.5281/zenodo.22274752` | 待填 10.5281/zenodo.22274752 |
| 预印本 DOI | v3 | 10.64898/2026.09.03.749125 |
| 三次 desk reject | Bioinformatics 状态卡 | CNSNT(08-25) / BiB BIB-26-2125(09-02) / BIOINF-2026-3215(09-22) |

---

## ⚠️ 投稿前必须处理的一件事（与 CL 无关但同样要紧）

**Zenodo V2 的描述字段写着 "manuscript submitted to GigaScience"，这是错的** ——
GigaScience 路线从未执行（状态卡明确："GigaScience 路线未执行，09-10 改投本刊"）。
这条描述是公开可见的元数据，若编辑或读者比对，会与 CL 中的披露不一致。

⇒ **发 Zenodo V3 时必须改掉这句**，并顺便把 v1.1.0 的判据协议与 CLI 一起打包进去。
（Zenodo V3 需老高本人网页端操作；V2 是手工上传 zip，非 GitHub 自动发版。）
