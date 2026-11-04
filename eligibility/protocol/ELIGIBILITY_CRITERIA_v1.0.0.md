# AnnoAudit confirmatory-eligibility criteria

**Version 1.0.0** · **2026-11-04** · SHA-256 `244a46300fdcb1ffd7da18a8fd3c0345d3dbd93dbf78dfed1504c331754728dd`

> AnnoAudit confirmatory-eligibility criteria for neuron-level inference in public single-cell atlases

This file is the human-readable companion to `criteria_v1.0.0.json`, which is the machine-readable single source of truth. The `annoaudit` command-line tool reads the JSON file, so the tool and this document cannot disagree. Cite the version and date above.

## Scope

Decides whether a public single-cell atlas can support a confirmatory, neuron-level (or other cell-type-level) claim about a gene's temporal dynamics, on the basis of the deposited count matrix and the declared biological replicates.

**Out of scope:**
- Judging the correctness of a depositing laboratory's cell-type annotation
- Any claim about resources whose matrices were not examined
- Choosing a marker panel appropriate for a non-neuronal cell type (the criteria below are written for excitatory neurons; the neuronal-signal clause must be re-derived for other lineages)

## The four criteria

A dataset is eligible only if it satisfies **every** criterion (combination rule: `ALL`).

| # | Criterion | Threshold | Label-dependent | Role |
|---|---|---|---|---|
| C1 | neuronal-signal sufficiency | At least 1% of all cells in the library reach >= 10 counts for a strict canonical neuronal gene. | no | gate |
| C2 | replicate adequacy | At least 100 marker-confirmed excitatory neurons per biological replicate, and at least two biological replicates per group. | yes | decisive |
| C3 | depth comparability | The ratio of the maximum to the minimum group median UMI count is no greater than 1.5. | no | precondition |
| C4 | definition-stability | The direction of the effect is concordant across all five pre-specified nested excitatory-neuron definitions (D1-D5). | yes | gate |

### C1 · neuronal-signal sufficiency

**Statement.** At least 1% of all cells in the library reach >= 10 counts for a strict canonical neuronal gene.

**Rationale.** The count matrix must contain neuronal signal for a neuron-level analysis to mean anything. Computed on the whole library and without reference to any label, so that ambient contamination of a labelled subset cannot influence it.

### C2 · replicate adequacy

**Statement.** At least 100 marker-confirmed excitatory neurons per biological replicate, and at least two biological replicates per group.

**Rationale.** Twice the ~50 cells per sample-cell-type combination at which aggregate profiles become noisy (muscat). We take twice that margin because the design has only three replicates per group. This is the criterion that carries the decision: applied alone it selects only the atlas that passes, at every threshold examined, and every criterion combination that includes it gives the same verdict.

### C3 · depth comparability

**Statement.** The ratio of the maximum to the minimum group median UMI count is no greater than 1.5.

**Rationale.** Depth differences generate apparent detection-rate differences for sparsely expressed transcripts, so a comparison across groups of unequal depth is not interpretable without depth matching. In the atlas set examined here this criterion has limited discriminating power (three of five atlases pass), although it does independently reject one of them (GSE247339, max/min group median UMI ratio 1.61). It is retained because it is a precondition for the depth-matched comparison, not because it excludes anything.

### C4 · definition-stability

**Statement.** The direction of the effect is concordant across all five pre-specified nested excitatory-neuron definitions (D1-D5).

**Rationale.** Applied to the result rather than to the dataset. Because the five nested definitions share one marker, agreement among them is partly tautological; three structurally independent definitions (N1-N3) were added afterwards as corroboration and are reported separately from this criterion.

## How the criteria combine

A dataset is eligible only if it satisfies every criterion. Criteria are conjoined, not interchangeable: in the atlas set examined, C1 alone would admit the atlas with the highest ambient contamination, and C3 alone admits three of the five atlases assessed.

**Threshold robustness.** verdict unchanged for the adopted combination at all three thresholds (counts thresholds [5, 10, 20]; Supplementary Table S5).

**Combination robustness.** Every single criterion and every two-criterion combination was tested at three neuronal-signal thresholds. All combinations containing C2 select the same single atlas at every threshold; no combination without C2 is stable. Those alternative combinations were not adopted. (Supplementary Table S5b; script `criteria_combination_robustness.py`).

## Excitatory-neuron definitions

All definitions are applied to the deposited count matrix with the gene symbols as written by the depositor; symbol case is matched case-insensitively. 'pan-neuronal detected' means max(count) > 0 over the pan-neuronal set; 'microglia not dominant' means max(microglial count) <= max(pan-neuronal count).

### Nested family (specified before the analyses reported in the manuscript)

| ID | Rule | Note |
|---|---|---|
| D1 | `Slc17a7 > 0` |  |
| D2 | `Slc17a7 >= 1` |  |
| D3 | `Slc17a7 >= 3` |  |
| D4 | `Slc17a7 >= 3 AND pan-neuronal detected AND microglia not dominant` | reference definition |
| D5 | `Slc17a7 >= 5 AND pan-neuronal detected AND microglia not dominant` |  |

### Structurally independent family (added afterwards)

**When and why.** after the nested family had shown 15/15 concordance, as a check on whether that concordance was an artefact of the shared marker. **Evidential weight.** independent corroboration, not a pre-specified criterion.

| ID | Rule | What makes it independent |
|---|---|---|
| N1 | `Snap25 >= 10 AND pan-neuronal detected AND microglia not dominant` | uses no Slc17a7 information at all |
| N2 | `Slc17a7 >= 3 AND Neurod6 > 0 AND microglia not dominant` | adds a corticofugal marker that is not part of the pan-neuronal set |
| N3 | `Slc17a7 >= 1 AND Camk2a > 0 AND pan-neuronal detected AND microglia not dominant AND Gad1 == 0 AND Gad2 == 0` | introduces an orthogonal axis: explicit exclusion of inhibitory neurons |

### Gene sets

- pan-neuronal: `Rbfox3`, `Snap25`, `Syp`, `Map2`, `Nrn1`
- microglial: `C1qb`, `Tyrobp`, `Aif1`, `Hexb`, `P2ry12`

## Reporting requirements

- State the failure mode for every excluded dataset, not only that it was excluded
- State whether the failure mode is recoverable in principle, and what remedy it would require
- Report the neuronal-signal threshold sensitivity (>=5, >=10, >=20)
- Report which criterion or criteria carry the decision, so that a reader applying the rule to a new dataset knows what to check first
- Report the biological replicate identity and count for each group; the replicate, not the cell, is the unit of analysis

## Provenance note

These criteria are written for datasets whose deposited matrices are accompanied by per-cell gene counts and for which the biological replicate of each cell can be stated. They carry no claim about the annotations of any depositing laboratory, and they should not be read as evaluating the correctness of any external cell-type annotation.

## Citation

- Manuscript DOI: pending
- Preprint (self-audit case): `10.64898/2026.09.03.749125`
- Repository: `github.com/ronnygao/AnnoAudit`
- Zenodo concept DOI: `10.5281/zenodo.22274752`

*Generated from `criteria_v1.0.0.json`; file SHA-256 `244a46300fdcb1ffd7da18a8fd3c0345d3dbd93dbf78dfed1504c331754728dd`.*

---

*This file is a frozen protocol version. To propose a change, open an issue or submit a pull request; do not edit a released version in place — issue a new versioned file so that analyses citing v1.0.0 remain reproducible.*
