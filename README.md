# AnnoAudit — analysis scripts

This directory contains the Python scripts implementing the AnnoAudit protocol
described in "AnnoAudit: a marker-based protocol for auditing single-cell atlas
annotations reveals systematic, state-dependent annotation failure in a widely used traumatic
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

## v2 discrimination, ambient and cross-validation scripts (added 2026-09)

These scripts implement the additional analyses of the revised manuscript: the independent-gene discrimination protocol (D1), the ambient-RNA simulation, state-dependence (injury-time) stratification, the QC-stratified re-analysis of GSE330130, and the cross-validation audits behind the reported numbers.

```
python v18_final_audit.py                    # final full-label audit (definitive numbers:
                                             #   CEREBRI glutamatergic 97.8% marker-only /
                                             #   67.4% D1-confirmed; GSE330130 1.8-6.1%)
python v17b_ambient_simulation_D1.py         # formal D1 rule under oligodendrocyte-dominant
                                             #   ambient (Fig 9; replaces the v17 draft)
python v17_ambient_simulation.py             # ambient gradient: marker-only C1 saturates ~100%
                                             #   when ambient >= 10% of the UMI pool
python discriminate_contamination.py         # D1 classifier (independent genes, tau = 1.5)
python classify_all_labels.py                # per-cell identity across all official labels
python survey_cerebri_all_labels.py          # C1 audit, all CEREBRI official labels
python survey_gse330130_all_labels.py        # C1 audit, all GSE330130 official labels
python time_structure_check.py               # contamination by injury time (state-dependent
                                             #   annotation failure, Table 8)
python qc_stratify_check.py                  # GSE330130 library-quality stratification
                                             #   (ambient-dilution signature, R4)
python crossval_static.py                    # static self-consistency audit
python crossval_tau.py                       # tau sensitivity (1.0-3.0), R1
python crossval_time_disc.py                 # D1 time-stratification cross-validation
```

Run order follows the manuscript Methods: `v18_final_audit.py` for the headline audit numbers, `v17b_ambient_simulation_D1.py` + `v17_ambient_simulation.py` for the ambient benchmark, `crossval_*.py` for the robustness statements, `time_structure_check.py` for the state-dependence results.

## Version note

`v14_*` scripts implement the protocol baseline used across the manuscript (32-gene marker panel; C1 argmax threshold > 0; C2 margin gate tau = 1.5; C3 KMeans; C4 applicability-gated CellTypist). The v17/v17b/v18 + discrimination/cross-validation scripts above (2026-09) supersede the earlier `v14-*` "final" runs for the revised manuscript's discrimination, ambient and state-dependence analyses. Earlier v1-v13 scripts were exploratory and are not part of the submission.
