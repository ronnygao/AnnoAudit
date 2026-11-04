# AnnoAudit — confirmatory-eligibility criteria for neuron-level inference

**Which public single-cell atlases can support a confirmatory, neuron-level claim?**

`annoaudit` implements the confirmatory-eligibility criteria described in the AnnoAudit
manuscript. It answers a question that is routinely conflated with label quality: given a
public atlas and a deposited count matrix, is that atlas able to support a claim about a
gene's temporal dynamics in neurons — or is it shallow, ambient-dominated, or too
under-replicated to carry such a claim?

The criteria live in a **versioned, dated protocol file**, not in the code:
`protocol/criteria_v1.0.0.json`. The tool reads them from there, so the tool and the
protocol can never drift apart, and the protocol file itself can be hashed and cited.

---

## Install

```bash
pip install -e .
```

Requires Python ≥ 3.9. **No third-party dependencies** — standard library only.

## Quick start

```bash
# 1. check the criteria you are applying
annoaudit criteria show

# 2. confirm the tool reproduces the published verdicts
annoaudit selftest

# 3. apply the criteria to your own dataset
annoaudit eligibility --table cells.csv --name MY_ATLAS --out out/

# 4. check whether the verdict depends on which criterion you use
annoaudit robustness --table cells.csv --name MY_ATLAS --out out/
```

## Input format

A flat table, **one cell per row**, with a header line. CSV or TSV (delimiter is
auto-detected). Required and optional columns are matched case-insensitively:

| Purpose | Accepted column names | Required |
|---|---|---|
| Biological replicate | `sample`, `gsm`, `replicate`, `orig.ident`, `sample_id`, `library`, `donor` | **yes** |
| Total UMI count | `umi`, `nCount`, `nCount_RNA`, `total_counts`, `counts`, `n_counts` | **yes** |
| Condition / group | `cond`, `condition`, `group`, `timepoint`, `time`, `arm` | no (needed for the depth ratio) |
| Marker gene counts | any column whose name is the gene symbol (`Slc17a7`, `Snap25`, `Ttr`, …) | no (but needed for a meaningful verdict) |

Marker columns used by the criteria:

- pan-neuronal: `Rbfox3`, `Snap25`, `Syp`, `Map2`, `Nrn1`
- microglial: `C1qb`, `Tyrobp`, `Aif1`, `Hexb`, `P2ry12`
- strict neuronal signal: `Snap25`, `Slc17a7`, `Rbfox3`, `Syt1`, `Stmn2`, `Nrn1`
- excitatory marker: `Slc17a7`

`Ttr` is read if present and is reported as an ambient diagnostic.

## Output

`annoaudit eligibility` writes two files into `--out`:

- **`eligibility_verdict.csv`** — one row per dataset: every criterion, the measured
  quantities behind it, the verdict, and — for excluded datasets — the failure mode, its
  diagnostic signature, whether it is recoverable in principle, and what remedy it needs.
- **`eligibility_report.md`** — the same information in prose, for reading.

`annoaudit robustness` writes **`criteria_combination_robustness.csv`**, which reports, for
every single criterion and every two-criterion combination, which datasets pass at each of
the three neuronal-signal thresholds. This is the check that answers "is the verdict an
artefact of the criteria you chose?".

## The four criteria

| # | Criterion | Threshold | Depends on labels? | Role |
|---|---|---|---|---|
| C1 | neuronal-signal sufficiency | ≥ 1% of all cells reach ≥ 10 counts for a strict neuronal gene | no | gate |
| C2 | replicate adequacy | ≥ 100 marker-confirmed excitatory neurons per biological replicate, and ≥ 2 replicates per group | yes | **decisive** |
| C3 | depth comparability | max/min group median UMI ratio ≤ 1.5 | no | precondition |
| C4 | definition-stability | effect direction concordant across the five pre-specified nested definitions D1–D5 | yes | gate |

A dataset is eligible only if it satisfies **all** criteria (conjunction, not a vote).

### Why C2 carries the decision

In the six-atlas set this work was developed on:

- **C2 alone** selects only the atlas that passes, at every threshold examined (≥5, ≥10, ≥20).
- **C1 alone** admits two atlases at ≥ 10 counts, three at ≥ 5 counts and one at ≥ 20 counts —
  it is threshold-sensitive on its own.
- **C3 alone** admits three of the five atlases assessed and so has limited discriminating
  power here, although it does independently reject one of them (GSE247339, whose max/min
  group median UMI ratio is 1.61).
- Every combination **containing** C2 selects the same single atlas at all three thresholds;
  no combination without it is stable.

The conjunction is therefore not decoration: under C1 alone, the atlas with the highest
ambient contamination (transthyretin positive in 33.0% of cells, maximum count 44,375) would
be admitted. We report this asymmetry rather than describing the four criteria as
interchangeable, because it tells an applying reader which criterion to check first.

## A note on the order of operations

The three *structurally independent* definitions (N1–N3) were **added after** the five nested
definitions (D1–D5) had already shown concordance in all fifteen comparisons. They were
introduced to test whether that concordance was an artefact of the shared marker. They are
reported as **independent corroboration**, not as a pre-specified criterion, and the protocol
file records this explicitly. We prefer this to asserting a pre-specification we could not
demonstrate: a reader can rerun `annoaudit robustness` and see for themselves whether the
verdict moves.

## Self-test

```bash
annoaudit selftest
```

The package ships a **fingerprinted golden result** (`annoaudit/golden_verdicts.json`)
recording the verdicts computed on six public TBI atlases on 2026-11-04 from the deposited
matrices. The self-test verifies:

1. the criteria file parses and passes structural validation;
2. the criteria file hash matches the value recorded in the golden result;
3. every per-atlas verdict value matches;
4. the four failure modes are mutually distinguishable;
5. C2 is recorded as the decisive criterion;
6. **end-to-end**: the engine recomputes GSE277487 from the per-cell table and reproduces
   every number.

Step 6 is the one that matters. A tool nobody can verify is not evidence.

## Provenance and scope

These criteria bear on whether a deposited matrix can support a neuron-level claim. They
carry **no claim about the correctness of any depositing laboratory's cell-type annotation**,
and they say nothing about resources whose matrices were not examined.

The criteria are written for excitatory neurons. The neuronal-signal clause (C1) must be
re-derived if the target lineage differs.

## Citation

If you use these criteria, please cite the AnnoAudit manuscript and the protocol version
you applied (`criteria_vX.Y.Z.json` carries its own date and the tool reports the file's
SHA-256 prefix).

- Repository: `github.com/ronnygao/AnnoAudit`
- Zenodo concept DOI: `10.5281/zenodo.22274752`
- Preprint (self-audit case): `doi:10.64898/2026.09.03.749125`

## License

MIT.
