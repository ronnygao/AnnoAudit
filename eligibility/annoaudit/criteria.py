"""
annoaudit.criteria
=====================================================================
资格判据的单一真相源。

设计原则（为什么判据放在 JSON 里而不是代码里）：
  审稿人/读者拿到的是「一份带日期的判据文件」，不是「一段代码」。
  如果阈值写死在代码里，那么"我们用的判据是什么"只能靠读代码确认；
  如果判据在一个带版本号、带日期的 JSON 里，那么：
    ① 工具与协议永远不会不一致（工具从 JSON 读）
    ② 可以对判据文件本身做哈希，写进正文
    ③ 将来改判据 = 发新版本文件，旧的仍可引用

这个模块只负责：
  - 载入 criteria_vX.Y.Z.json
  - 提供判据的机器可读表示
  - 提供版本与哈希查询

具体判定逻辑在 eligibility.py。
"""
from __future__ import annotations

import hashlib
import json
import os
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

_HERE = os.path.dirname(os.path.abspath(__file__))
_PROTOCOL_DIR = os.path.join(os.path.dirname(_HERE), "protocol")

DEFAULT_CRITERIA_VERSION = "1.0.0"


def protocol_path(version: str = DEFAULT_CRITERIA_VERSION) -> str:
    """返回指定版本判据文件的路径。"""
    return os.path.join(_PROTOCOL_DIR, f"criteria_v{version}.json")


def load_criteria(version: str = DEFAULT_CRITERIA_VERSION) -> Dict[str, Any]:
    """载入判据文件并校验其结构。"""
    path = protocol_path(version)
    if not os.path.exists(path):
        avail = available_versions()
        raise FileNotFoundError(
            f"判据文件不存在：{path}\n可用版本：{avail or '（无）'}"
        )
    with open(path, "r", encoding="utf-8") as fh:
        crit = json.load(fh)
    _validate(crit)
    return crit


def available_versions() -> List[str]:
    """列出本地可用的判据版本。"""
    if not os.path.isdir(_PROTOCOL_DIR):
        return []
    out = []
    for fn in sorted(os.listdir(_PROTOCOL_DIR)):
        if fn.startswith("criteria_v") and fn.endswith(".json"):
            out.append(fn[len("criteria_v"):-len(".json")])
    return out


def criteria_hash(version: str = DEFAULT_CRITERIA_VERSION) -> str:
    """判据文件的 SHA-256（十六进制前 16 位），用于正文引用与版本核验。"""
    path = protocol_path(version)
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()[:16]


def _validate(crit: Dict[str, Any]) -> None:
    """结构校验：判据文件缺字段会让工具静默给出错误结论，必须在载入时就失败。"""
    required_top = ["criteria_version", "criteria_date", "criteria", "combination"]
    for k in required_top:
        if k not in crit:
            raise ValueError(f"判据文件缺少必需字段：{k}")
    ids = [c.get("id") for c in crit["criteria"]]
    if sorted(ids) != ["C1", "C2", "C3", "C4"]:
        raise ValueError(f"判据 ID 应为 C1–C4，实得 {ids}")
    for c in crit["criteria"]:
        for k in ("id", "name", "statement", "threshold", "rationale"):
            if k not in c:
                raise ValueError(f"判据 {c.get('id')} 缺少字段 {k}")
    if crit["combination"].get("rule") != "ALL":
        raise ValueError(
            "本工具实现的是合取判据（rule=ALL）。若要把判据改为可选项，"
            "必须同时改 eligibility.py，否则工具行为与文件声明不一致。"
        )


@dataclass
class Criterion:
    """单条判据的机器可读表示。"""
    id: str
    name: str
    statement: str
    rationale: str
    threshold: Dict[str, Any]
    role: str = "gate"
    depends_on_labels: bool = True

    @property
    def is_decisive(self) -> bool:
        return self.role == "decisive"

    def describe(self) -> str:
        return f"{self.id} [{self.name}] — {self.statement}"


@dataclass
class Definitions:
    """兴奋性神经元定义族。"""
    nested: Dict[str, Any] = field(default_factory=dict)
    independent: Dict[str, Any] = field(default_factory=dict)
    pan_neuronal: List[str] = field(default_factory=list)
    microglial: List[str] = field(default_factory=list)

    @property
    def all_ids(self) -> List[str]:
        return list(self.nested.keys()) + list(self.independent.keys())

    @staticmethod
    def _definition_ids(block: Dict[str, Any]) -> Dict[str, Any]:
        """从定义块里挑出真正的定义条目（排除 pre_specified / added / evidential_weight
        这类元字段）。JSON 里的元字段与定义混在同一层，不做区分会把它们当成定义。"""
        meta = {"pre_specified", "added", "evidential_weight", "note", "gene_sets"}
        return {k: v for k, v in block.items() if k not in meta}


@dataclass
class CriteriaSpec:
    """完整判据规格。"""
    raw: Dict[str, Any]
    version: str
    date: str
    criteria: List[Criterion]
    definitions: Definitions
    combination_rule: str

    def by_id(self, cid: str) -> Optional[Criterion]:
        for c in self.criteria:
            if c.id == cid:
                return c
        return None

    @property
    def decisive_criterion(self) -> Optional[Criterion]:
        for c in self.criteria:
            if c.is_decisive:
                return c
        return None

    @property
    def hash(self) -> str:
        return criteria_hash(self.version)


def spec(version: str = DEFAULT_CRITERIA_VERSION) -> CriteriaSpec:
    """载入并结构化判据规格。"""
    crit = load_criteria(version)
    criteria = [
        Criterion(
            id=c["id"],
            name=c["name"],
            statement=c["statement"],
            rationale=c["rationale"],
            threshold=c["threshold"],
            role=c.get("role", "gate"),
            depends_on_labels=c.get("depends_on_labels", True),
        )
        for c in crit["criteria"]
    ]
    ed = crit.get("excitatory_definitions", {})
    gs = ed.get("gene_sets", {})
    defs = Definitions(
        nested=Definitions._definition_ids(ed.get("nested", {})),
        independent=Definitions._definition_ids(ed.get("independent", {})),
        pan_neuronal=gs.get("pan_neuronal", []),
        microglial=gs.get("microglial", []),
    )
    return CriteriaSpec(
        raw=crit,
        version=crit["criteria_version"],
        date=crit["criteria_date"],
        criteria=criteria,
        definitions=defs,
        combination_rule=crit["combination"]["rule"],
    )
