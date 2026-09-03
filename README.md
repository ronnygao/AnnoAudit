# AnnoAudit — analysis scripts

This directory contains the Python scripts implementing the AnnoAudit protocol
described in "AnnoAudit: a marker-based protocol for auditing single-cell atlas
annotations reveals annotation-driven artifacts in a widely used traumatic
brain injury resource" (accompanying manuscript). All code is released under
the MIT license (see LICENSE).

All scripts were executed with Python 3.10+ (conda) using scanpy 1.12.x,
anndata 0.13.x, scikit-learn, scipy, numpy, and pandas. No additional
domain-specific packages are required. Exact versions are pinned in
`requirements.txt` (this directory); install with:

```
pip install -r requirements.txt
```

Runtime is a few minutes to ~40 minutes per script depending on dataset size.

## Script map (manuscript section → script)

| Manuscript section | Script |
|---|---|
| Methods — The AnnoAudit protocol / composite ACS | `v14_cerebri_final.py` (mouse), `v14_gse330130_final.py` (human) |
| Methods — Simulation benchmark | `v12_simulation.py` |
| Methods — Marker-based cell-type annotation (leave-one-marker-out) | `v9b_marker_robustness.py` |
| Methods — Marker-based annotation (module-scoring cross-check) | `v9c_scorer_sensitivity.py` |
| Methods — Signal-migration / provenance fingerprint (Fig 3) | `v10_provenance.py` |
| Methods — ACS (C2 margin-gate tau sensitivity, Supplementary) | `v15_tau_sensitivity.py` |
| Methods — Discriminative performance / biological anchor / transferability (Supplementary Table S4) | `v16_methodology_enhancement.py` |
| Results R4 — panel sensitivity (32- vs 20-gene) | `v14_diag_panel.py` |

## Input data (all public)

- **CEREBRI** (GSE269748): 10x Genomics matrix files
  (`*_matrix.mtx.gz`, `*_barcodes.tsv.gz`, `*_features.tsv.gz`) and the official
  cell annotation table, downloaded from GEO/NCBI and from the CEREBRI
  interactive portal. Local path used in scripts:
  `.../2026-07-22-16-57-18/cerebri_data/raw` + `v3分析结果/cell_annotations.csv`.
- **GSE330130**: human ALS lumbar spinal cord snRNA-seq AnnData object,
  `GSE330130_LSC_SNRNASEQ_COUNTS.h5ad` (GEO). Local path:
  `.../TBI分子时钟/数据/`.
- **CellTypist models**: `Developing_Mouse_Brain.pkl`, `Adult_Human_MTG.pkl`,
  `Adult_Human_PrefrontalCortex.pkl` (CellTypist model zoo;
  `~/.celltypist/data/models/`).

The scripts contain absolute local paths. To reproduce on another machine,
download the datasets above, update the path constants at the top of each
script, and run:

```
python v14_cerebri_final.py      # CEREBRI final audit -> final_audit.json
python v14_gse330130_final.py    # GSE330130 final audit -> final_audit.json
python v12_simulation.py         # contamination simulation calibration
python v15_tau_sensitivity.py    # C2 margin-gate tau sensitivity (both datasets)
python v16_methodology_enhancement.py  # AUROC/permutation/anchor/transferability
```

Per-cell audit results (marker typing, cluster assignment, module scores,
margin distributions) are written as CSV/JSON next to each script output
directory; key aggregates are reported to stdout and in `final_audit.json`.

The `v16_methodology/` subdirectory contains `summary.json`, the key output of
`v16_methodology_enhancement.py` (AUROC, permutation *P* values, and 16-gene
anchor results reported in Supplementary Table S4).

## Version note

`v14_*` scripts implement the final protocol version used for all numbers in
the manuscript (32-gene marker panel; C1 argmax threshold > 0; C2 margin gate
τ = 1.5; C3 KMeans; C4 applicability-gated CellTypist). Earlier v1–v13 scripts
were exploratory and are not part of the submission.
