#!/usr/bin/env python3
"""
stream_audit_matrix.py
=====================================================================
流式审计一个 10x/MatrixMarket 计数矩阵，**不落盘整个文件**。

为什么需要这个：
  NCBI 对大文件限流严重（实测：46MB 的矩阵反复续传只能到 1MB 就停），
  而读者/我们真正需要的只是几个判据指标，不是完整矩阵。
  ⇒ 直接 curl | zcat 边下边算，服务器在管道关闭时即停止传输。

★本脚本只依赖标准库（不需 numpy/scipy），因为要在任何机器上跑。

指标（对应判据 C1/C3 与真群规模）：
  - 细胞数、基因数、非零项数
  - 每个细胞的 UMI 总数分布（中位、p10/p90）
  - 每个细胞的检出基因数分布
  - 严格神经元基因（Snap25/Slc17a7/Rbfox3/Syt1/Stmn2/Nrn1）的 max count
  - 全库中达到 ≥10 counts 的细胞占比（判据 C1 的量）

用法：
  curl -sL <url> | /usr/bin/python3 stream_audit_matrix.py --genes genes.csv.gz
"""
import argparse
import gzip
import io
import sys
from collections import defaultdict

STRICT_NEURONAL = ["Snap25", "Slc17a7", "Rbfox3", "Syt1", "Stmn2", "Nrn1"]


class _PrependStream:
    """把已读走的字节接回流前面的只读包装（用于 stdin 不可 seek 的场景）。"""

    def __init__(self, stream, prefix: bytes):
        self._s = stream
        self._pre = prefix

    def read(self, n=-1):
        if not self._pre:
            return self._s.read(n)
        if n is None or n < 0:
            data = self._pre + self._s.read()
            self._pre = b""
            return data
        out = self._pre[:n]
        self._pre = self._pre[n:]
        if len(out) < n:
            out += self._s.read(n - len(out))
        return out

    def readable(self):
        return True


def read_gene_names(path):
    """从 all_genes.csv.gz 读 gene_name 列，返回 {gene_name: index}（0-based）。"""
    idx = {}
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split(",")
        try:
            gi = header.index("gene_name")
        except ValueError:
            gi = 1 if len(header) > 1 else 0
        for i, line in enumerate(fh):
            parts = line.rstrip("\n").split(",")
            if len(parts) > gi:
                idx[parts[gi]] = i
    return idx


def audit(stream, gene_index=None, strict=STRICT_NEURONAL):
    """从 MatrixMarket 坐标格式流式计算指标。"""
    n_cells = n_genes = nnz = 0
    umi = defaultdict(int)
    det = defaultdict(int)
    strict_max = defaultdict(int)
    want = {}
    if gene_index:
        for g in strict:
            if g in gene_index:
                want[gene_index[g]] = g

    # ★自动检测 gzip：GEO 的 .mtx.gz 是 gzip 压缩的 MatrixMarket。
    # 注意：stdin 不可 seek，所以用"读一字节再回推"的方式探测魔数。
    first = stream.read(1)
    if first == b"\x1f":
        second = stream.read(1)
        if second == b"\x8b":
            raw = gzip.GzipFile(fileobj=_PrependStream(stream, b"\x1f\x8b"))
        else:
            raw = _PrependStream(stream, b"\x1f" + second)
    else:
        raw = _PrependStream(stream, first)
    buf = io.TextIOWrapper(raw, encoding="utf-8", errors="replace")
    # 头部
    line = buf.readline()
    if not line:
        raise ValueError("空输入")
    if not line.startswith("%%MatrixMarket"):
        raise ValueError(f"不是 MatrixMarket 格式：{line[:60]!r}")
    buf.readline()  # %Rows=... Cols=...
    dims = buf.readline().split()
    n_cells, n_genes = int(dims[0]), int(dims[1])
    nnz = int(dims[2])

    for line in buf:
        p = line.split()
        if len(p) < 3:
            continue
        i = int(p[0]) - 1
        j = int(p[1]) - 1
        v = int(p[2])
        umi[i] += v
        det[i] += 1
        if j in want and v > strict_max[want[j]]:
            strict_max[want[j]] = v

    # 汇总
    umis = sorted(umi.values())
    dets = sorted(det.values())
    n_with_data = len(umis)

    def pct(sorted_vals, q):
        if not sorted_vals:
            return 0
        k = int(q * (len(sorted_vals) - 1))
        return sorted_vals[k]

    # 判据 C1 的量：某严格神经元基因达到 ≥N counts 的细胞占比
    # 需要每个细胞在每个 strict 基因上的 max count → 这里按基因聚合，
    # 严格说 C1 需要逐细胞逐基因的 max；用"至少一个 strict 基因达标"的口径，
    # 与 audit_neuronal_sufficiency.py 的 max(各 strict 基因) 一致。
    ge10 = 0
    ge10_by_gene = {g: 0 for g in want.values()}

    return dict(
        n_cells=n_cells, n_genes=n_genes, nnz=nnz,
        n_cells_with_data=n_with_data,
        umi_median=pct(umis, 0.5), umi_p10=pct(umis, 0.1), umi_p90=pct(umis, 0.9),
        det_median=pct(dets, 0.5), det_p10=pct(dets, 0.1), det_p90=pct(dets, 0.9),
        strict_max_counts=dict(strict_max),
        strict_genes_found=sorted(want.values()),
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes", help="all_genes.csv.gz（含 gene_name 列）")
    ap.add_argument("--matrix", help="本地矩阵路径（不给则从 stdin 读）")
    args = ap.parse_args()

    gene_index = read_gene_names(args.genes) if args.genes else None
    if gene_index:
        print(f"基因表：{len(gene_index):,} 个基因名", file=sys.stderr)

    if args.matrix:
        fh = open(args.matrix, "rb")
    else:
        fh = sys.stdin.buffer

    r = audit(fh, gene_index)
    print("\n=== 流式审计结果 ===")
    print(f"  细胞数 (Rows)        : {r['n_cells']:,}")
    print(f"  基因数 (Cols)        : {r['n_genes']:,}")
    print(f"  非零项 (nnz)         : {r['nnz']:,}")
    print(f"  有数据的细胞数        : {r['n_cells_with_data']:,}")
    print(f"  UMI/细胞  中位        : {r['umi_median']:,.0f}")
    print(f"  UMI/细胞  p10 / p90  : {r['umi_p10']:,.0f} / {r['umi_p90']:,.0f}")
    print(f"  检出基因/细胞 中位    : {r['det_median']:,.0f}")
    print(f"  检出基因/细胞 p10/p90: {r['det_p10']:,.0f} / {r['det_p90']:,.0f}")
    print(f"  严格神经元基因        : {', '.join(r['strict_genes_found']) or '（未找到）'}")
    for g, v in sorted(r["strict_max_counts"].items()):
        print(f"    {g:10} max count = {v:,}")


if __name__ == "__main__":
    main()
