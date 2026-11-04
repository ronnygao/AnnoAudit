"""
AnnoAudit — confirmatory-eligibility criteria for neuron-level inference
in public single-cell atlases.

AnnoAudit 判定一个公开单细胞图谱能否支撑**确认性的、神经元层面**的
结论（例如某基因的时间动态），判据全部写在带版本号与日期的协议文件里
（protocol/criteria_vX.Y.Z.json），工具从该文件读取，因此工具与协议
永远一致。

快速开始：
    pip install -e .
    annoaudit criteria show
    annoaudit selftest
    annoaudit eligibility --table cells.csv --name MY_DATASET --out out/

设计取舍：
  - 只读扁平表（每细胞一行），不要求 scanpy / anndata。
  - 判据在 JSON 里，不写死在代码里（见 criteria.py 的说明）。
  - 内置带指纹的黄金结果回归自测（selftest.py）。
"""
from .criteria import (  # noqa: F401
    DEFAULT_CRITERIA_VERSION,
    CriteriaSpec,
    available_versions,
    criteria_hash,
    load_criteria,
    protocol_path,
    spec,
)
from .eligibility import (  # noqa: F401
    DatasetResult,
    combination_robustness,
    evaluate,
    read_table,
)

__version__ = "1.0.0"
__all__ = [
    "spec", "load_criteria", "criteria_hash", "available_versions",
    "protocol_path", "CriteriaSpec", "DEFAULT_CRITERIA_VERSION",
    "evaluate", "combination_robustness", "read_table", "DatasetResult",
    "__version__",
]
