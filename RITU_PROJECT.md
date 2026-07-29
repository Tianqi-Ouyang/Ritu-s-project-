# Ritu Project — Local/Early Stage Cohort

Creatinine-to-cystatin C discordance and adverse events in patients receiving
platinum-based chemotherapy, restricted to **locally staged / early-stage disease**.

This project is a direct derivative of the Jiaxuan project. Every analysis,
table, figure, model specification and output format is carried over unchanged.
The single substantive difference is the cohort: all advanced/metastatic
patients are excluded.

---

## 1 Cohort

### 1.1 Restriction

```r
final_master <- final_master %>%
  filter(tumor_stage == "Local/Early") %>%
  # Clear the now-empty Advanced/Metastatic level (and any covariate level left
  # with zero patients) so model matrices are not rank-deficient. cmp_status_*
  # factors are excluded: their 0/1/2 coding is relied on downstream.
  mutate(across(where(is.factor) & !starts_with("cmp_status_"), droplevels))
```

`tumor_stage` is a REDCap variable with two levels:

| REDCap code | Label | Stage | Ritu project |
|---|---|---|---|
| 1 | `Local/Early` | I–II (local / neoadjuvant / early) | **kept** |
| 2 | `Advanced/Metastatic` | III–IV | **excluded** |

The restriction is applied at the very end of the data-management pipeline,
after every derived variable has been created. No derived variable in the
pipeline is defined relative to the cohort (there are no median splits,
quantile categories, z-scores, or imputation models), so restricting after
derivation is numerically identical to restricting before it.

### 1.2 Cohort sizes

| Cohort | Jiaxuan (all stages) | Ritu (Local/Early) |
|---|---|---|
| Whole cohort | 943 | **323** |
| Carboplatin | 596 | **190** |
| Carboplatin AUC ≥ 3 | 463 | **143** |
| Cisplatin (not analysed separately) | 347 | 133 |

### 1.3 Exclusion flow

Inherited from the Jiaxuan pipeline, then one additional step:

1. REDCap records restricted to EMPIs present in `N=596 for tianqi.xlsx`
2. Patients with no matching RPDR lab records excluded
3. Missing baseline cystatin C (`cys_c_baseline` = NA) excluded
4. Missing platinum group (`platin_group` = NA) excluded
5. Record 911 dropped (all RPDR pre-baseline labs outside the 45-day window)
6. **New — `tumor_stage != "Local/Early"` excluded (620 patients)**

---

## 2 What changed relative to the Jiaxuan project

Four kinds of change were made. Nothing else in the analysis was touched.

### 2.1 Cohort filter

`ritu_data.qmd` applies the `tumor_stage == "Local/Early"` filter and writes
`ritu_final_master.rds`. The three analysis files read that file instead of
`jiaxuan_final_master.rds`.

### 2.2 `tumor_stage` removed as a covariate

Because `tumor_stage` is constant within this cohort it carries no information
and cannot be estimated. It is therefore removed from:

- the Table 1 row list (it would print a single 100 % row)
- the adjuster set of every Fine-Gray model (`.adj_set()`)
- the adjuster set of every cause-specific hazard model (`csh_model()`)
- the adjuster set of every Cox PH model for 90-day death (`cox_model()`)
- the adjuster set of every linear subHR / HR curve helper
- the missing-data variable list in the data-management QC section

Leaving it in place would produce a rank-deficient design matrix and either an
error or silently dropped `NA` coefficients.

Two further covariate changes follow from the same principle, and are confined to
the **AUC ≥ 3 page**:

- **`cancer_type_cat_7` levels are re-dropped after the cohort filter.**
  `ritu_data.qmd` calls `droplevels()` against the 323 Local/Early patients, but
  the AUC ≥ 3 filter empties the `GI` level again (0 of 143 patients), leaving an
  all-zero dummy column. `ritu_carbo_auc3.qmd` therefore re-applies
  `droplevels()` immediately after its own filter; `ritu_carbo.qmd` does the same
  for symmetry, though the carboplatin cohort happens to be full rank already.
- **`steroids` is dropped from the AUC ≥ 3 adjuster sets.** All 143 patients in
  that subgroup received steroids, so the variable is constant and collinear with
  the intercept. It is kept as a descriptive row in Table 1.

This was not cosmetic. Before the fix the AUC ≥ 3 design matrix had rank 19 of 21
columns, so `cmprsk::crr`'s `solve()` failed for **every** outcome — including
composite grade 2 AEs with 64 events — and the `tidycmprsk` formula interface
silently returned blank coefficient rows for `GI` and `steroids` instead of
erroring. The whole (323) and carboplatin (190) cohorts are full rank across
every model subset and keep the original adjuster set unchanged.

### 2.3 Sparse-outcome guards

Some outcomes have too few events in this cohort for a ~20-covariate model to be
estimated (see §3). An unfittable model raises an R error, which aborts the
entire Quarto render — so every table, summary-row and figure helper is wrapped:

```r
fg_model <- .wrap_tbl(fg_model)   # and .wrap_row / .wrap_plot for the others
```

On success the wrapper returns exactly what the Jiaxuan version returns. On
failure it emits a short "model not estimable" note in place of the table or
figure, records the event count, and lets the render continue. **No model
specification, adjuster set, or estimation setting is changed by the wrapper** —
it only decides what happens to an error that would otherwise be fatal.

The wrap is applied immediately after each function definition rather than at
the end of the chunk, because the Section 7 / 10 summary tables define their row
helpers and call them in the same chunk.

In the finished site the guard fires exactly once per page pair: only
thrombocytopenia grade 3 errors, in the carboplatin (3 events) and AUC ≥ 3 (0
events) summary tables. Everything else estimates.

That is not the same as everything else being *interpretable*. Where a model is
estimable but degenerate — 90-day death fits with 2 events in the whole cohort and
returns a confidence interval hundreds of digits wide — the guard does not fire
and the raw output is shown as-is. See §3.

### 2.4 Cohort-size prose

Hardcoded denominators in narrative text ("0/943 missing") were updated to the
Ritu denominators, and each page carries a header note stating the restriction,
its N, and which of its outcomes are underpowered.

### 2.5 Explicitly unchanged

- Primary exposure: `egfr_diff_pct_per10` — eGFR difference (%) per 10-percentage-point decline
- Sensitivity exposure: `egfr_diff_pct_sens_per10` (CKD-EPI 2012 CysC vs CKD-EPI 2021 Cr)
- Alternative predictor (carbo files): `dose_discrep_per25increase`
- Reference categories: `sex` = Female, `ecog_score_cat2` = "0"
- Competing-risks framework: death = competing event, `tidycmprsk::crr()`
- 90-day death modelled with Cox PH, not competing risks
- eGFR capping rules (only `ckd_epi_gfr_cre_cys_unindex` capped at 125)
- Plot trimming at the 2.5 / 97.5 percentile
- Word (`.docx`) table downloads and editable PowerPoint (`.pptx`) figure downloads
- Section numbering and ordering

---

## 3 Statistical caveats introduced by the restriction

These follow arithmetically from the smaller cohort and are worth reading before
interpreting the output. They are reported here rather than silently changing
the analysis, because the analysis was specified to match the Jiaxuan project.

Event counts below are the modelled counts (`cmp_status_* == 1`), i.e. after
death is recoded as the competing event — these are the numbers the Fine-Gray
models actually see, and they match the `Events` column of the rendered summary
tables.

| Outcome | Whole (323) | Carbo (190) | AUC ≥ 3 (143) |
|---|---|---|---|
| Composite grade 2 AEs | 149 | 81 | 64 |
| Anemia grade 2 | 138 | 79 | 62 |
| Overall hospitalisation | 68 | 36 | 25 |
| Platinum-related hospitalisation | 44 | 20 | 15 |
| Composite grade 3 AEs | 43 | 20 | 14 |
| Anemia grade 3 | 40 | 18 | 14 |
| Thrombocytopenia grade 2 | 37 | 14 | 9 |
| AKI (KDIGO) | 33 | 12 | 6 |
| AKI grade 2 | 26 | 9 | 6 |
| Hosp — infection/sepsis | 25 | 12 | 9 |
| Hosp — planned surgery/port | 22 | 13 | 9 |
| Hosp — immune-related & other | 15 | 10 | 7 |
| Hosp — failure to thrive | 14 | 7 | 4 |
| Hosp — nausea/vomiting/GI | 13 | 5 | 4 |
| Hosp — cancer progression | 13 | 6 | 4 |
| Thrombocytopenia grade 3 | 9 | 3 | **0** |
| Hosp — AKI | 7 | 4 | 1 |
| Hosp — electrolyte abnormalities | 5 | 4 | 3 |
| Hosp — cytopenia | 5 | 2 | 1 |
| Hosp — febrile neutropenia | 5 | 2 | 2 |
| **90-day death** | **2** | **1** | **0** |

AKI grade 3 and hospitalisation-for-immune-related-AE also have zero events, but
neither is modelled on its own in any section, so nothing changes there.

The multivariable models carry roughly 20 covariate degrees of freedom. The
conventional target is 10 events per covariate, so:

- **Grade 2 composite and anemia models** are adequately powered.
- **Grade 3, AKI, thrombocytopenia and hospitalisation models** are
  over-parameterised; multivariable estimates will be unstable with wide
  confidence intervals. Univariable columns remain interpretable.
- **The 90-day death models are not interpretable.** `coxph` does not error on
  these, so the guard never fires — it just returns degenerate numbers:

  | Cohort | Deaths | Reported subHR/HR (95 % CI) |
  |---|---|---|
  | Whole (323) | 2 | 18.04 (0.00, an upper limit hundreds of digits long) |
  | Carboplatin (190) | 1 | 25.19 (0.00, Inf) |
  | AUC ≥ 3 (143) | 0 | NA (NA, NA) |

  None of these rows carries information. The underlying fact — **2/323 (0.6 %)
  90-day mortality** — is expected for early-stage disease and is itself worth
  stating in the write-up.

- **Thrombocytopenia grade 3 is the one outcome the guard actually catches.**
  With 9 / 3 / 0 events it errors in the carboplatin and AUC ≥ 3 summary tables,
  which show a "not estimable" row instead. It still fits in the whole cohort.
- Death being near-absent also means the competing-risks correction has almost
  no effect in this cohort; Fine-Gray subdistribution hazards and
  cause-specific hazards will be nearly identical.

If a reduced adjuster set is wanted for the sparse outcomes, that is a
follow-up decision, not part of the "same as Jiaxuan" specification.

---

## 4 Repository layout

```
Ritu project/
├── RITU_PROJECT.md          this file
├── README.md                GitHub landing page
├── _quarto.yml              website config (output-dir: docs)
├── custom.scss              theme overrides (copied from Jiaxuan project)
├── index.qmd                project summary page
├── qmd/
│   ├── ritu_data.qmd        data management → ritu_final_master.rds
│   ├── ritu_whole.qmd       whole Local/Early cohort (N = 323)
│   ├── ritu_carbo.qmd       carboplatin (N = 190)
│   └── ritu_carbo_auc3.qmd  carboplatin AUC ≥ 3 (N = 143)
├── ritu_final_master.rds    analytic dataset (gitignored — contains PHI)
└── docs/                    rendered site, committed and served by GitHub Pages
    ├── index.html
    ├── qmd/*.html
    ├── tables/{whole,carbo,auc3}/*.docx
    └── pptx/{whole,carbo,auc3}/*.pptx
```

### Path conventions

`_quarto.yml` sets the knitr root directory to the `Ritu project/` folder.
Consequently:

- REDCap inputs are read as `../Data/<file>.csv` (they live in the parent
  `Cystain C redcap/` folder, shared with the Jiaxuan project)
- RPDR `.txt` pulls are read from their absolute Dropbox paths, unchanged
- `ritu_final_master.rds` is written to the `Ritu project/` root
- `docs/tables/…` and `docs/pptx/…` land inside this project, not the parent

---

## 5 Rendering

Render the whole site (data management runs first and produces the `.rds` the
analysis files need):

```bash
cd "/Users/to909/Desktop/Meg/Cystain C redcap/Ritu project" && quarto render
```

Render a single page:

```bash
cd "/Users/to909/Desktop/Meg/Cystain C redcap/Ritu project" && quarto render qmd/ritu_whole.qmd
```

`ritu_data.qmd` re-reads the full RPDR pull (≈340 MB of `.txt` files) and takes
several minutes. The analysis files render in well under a minute each once
`ritu_final_master.rds` exists.

---

## 6 Publishing

| | |
|---|---|
| Repository | `https://github.com/Tianqi-Ouyang/Ritu-s-project-` |
| Pages source | `main` branch, `/docs` folder |
| Live site | `https://tianqi-ouyang.github.io/Ritu-s-project-/` |

This mirrors the Jiaxuan project setup (`Cys-C-project-` →
`https://tianqi-ouyang.github.io/Cys-C-project-/`).

```bash
cd "/Users/to909/Desktop/Meg/Cystain C redcap/Ritu project" && quarto render && git add -A && git commit -m "Update rendered site" && git push
```

---

## 7 Privacy

`.gitignore` excludes `*.rds`, `*.csv`, `*.xlsx`, `*.pdf` and the `Data/`
folder, so no patient-level file is ever committed. The rendered site contains
only aggregate tables, model coefficients and figures — no row-level records
and no identifier columns (`EMPI`, `MRN`, names, dates of birth). `record_id`
appears only in the data-management diagnostics that list which records were
excluded; it is a study-internal sequence number, not a hospital identifier.
