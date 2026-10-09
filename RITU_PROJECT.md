# Ritu Project

Creatinine-to-cystatin C discordance and adverse events in patients receiving
platinum-based chemotherapy.

A direct derivative of the Jiaxuan project. Every analysis, table, figure, model
specification and output format is carried over unchanged. The project runs at
**two cohort levels off one pipeline**: the **entire cohort** (N = 943, all
tumour stages) and the **Local/Early** cohort (N = 322, advanced/metastatic
excluded), which was the original point of the project.

It has also gained, on request and documented in §2: KDIGO CKD stage in Table 1
(§2.5); creatinine sourced from the newer 2022P001091 RPDR pull (§2.6); the
long-horizon CKD / ESKD / AKI outcome family ported from the ICI-CKD project
(§2.7); and the entire-cohort page (§2.8). None of these alters the
Jiaxuan-derived models — the ported kidney outcomes are descriptive additions
that sit outside the 90-day competing-risks framework.

---

## 1 Cohort

The pipeline produces two analytic datasets in a single run. `ritu_full_master.rds`
(N = 943) is saved immediately before the stage restriction below;
`ritu_final_master.rds` (N = 322) immediately after it. Every variable definition
is identical between them.

### 1.1 The Local/Early restriction

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

| Page | Cohort | N |
|---|---|---|
| `ritu_entire.qmd` | all platinum, **all stages** | **943** |
| `ritu_whole.qmd` | all platinum, Local/Early | **322** |
| `ritu_carbo.qmd` | carboplatin, Local/Early | **189** |
| `ritu_carbo_auc3.qmd` | carboplatin AUC ≥ 3, Local/Early | **142** |

Entire cohort composition: 322 Local/Early + 621 Advanced/Metastatic;
596 carboplatin + 347 cisplatin. Within the Local/Early 322: 189 carboplatin,
133 cisplatin. See §2.8 for how the two levels differ analytically.

### 1.3 Exclusion flow

Inherited from the Jiaxuan pipeline, then one additional step:

1. REDCap records restricted to EMPIs present in `N=596 for tianqi.xlsx` → 946
2. Patients with no matching RPDR lab records excluded (2) → 944
3. Missing baseline cystatin C (`cys_c_baseline` = NA) excluded
4. Missing platinum group (`platin_group` = NA) excluded
5. Record 911 dropped (all RPDR pre-baseline labs outside the 45-day window)
   → **943, saved as `ritu_full_master.rds` (entire cohort)**
6. **`tumor_stage != "Local/Early"` excluded (621 patients)
   → 322, saved as `ritu_final_master.rds` (Local/Early cohort)**

---

## 2 What changed relative to the Jiaxuan project

Ten kinds of change were made. Nothing else in the analysis was touched.

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
  `ritu_data.qmd` calls `droplevels()` against the 322 Local/Early patients, but
  the AUC ≥ 3 filter empties the `GI` level again (0 of 142 patients), leaving an
  all-zero dummy column. `ritu_carbo_auc3.qmd` therefore re-applies
  `droplevels()` immediately after its own filter; `ritu_carbo.qmd` does the same
  for symmetry, though the carboplatin cohort happens to be full rank already.
- **`steroids` is dropped from the AUC ≥ 3 adjuster sets.** All 142 patients in
  that subgroup received steroids, so the variable is constant and collinear with
  the intercept. It is kept as a descriptive row in Table 1.

This was not cosmetic. Before the fix the AUC ≥ 3 design matrix had rank 19 of 21
columns, so `cmprsk::crr`'s `solve()` failed for **every** outcome — including
composite grade 2 AEs with 64 events — and the `tidycmprsk` formula interface
silently returned blank coefficient rows for `GI` and `steroids` instead of
erroring. The whole (322) and carboplatin (189) cohorts are full rank across
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

### 2.5 CKD stage added to Table 1

Two derived variables were added on request, in Section 5.6 of `ritu_data.qmd`:

| Variable | Source eGFR |
|---|---|
| `ckd_stage_cre` | `baseline_cre_egfr` — CKD-EPI 2021 creatinine, BSA-indexed |
| `ckd_stage_cre_cys` | `ckd_epi_gfr_cre_cys` — CKD-EPI 2012 combined, BSA-indexed |

Standard KDIGO G categories: G1 ≥ 90, G2 60–89, G3a 45–59, G3b 30–44, G4 15–29,
G5 < 15 mL/min/1.73m². The **indexed** eGFRs are used because that is the scale
KDIGO staging is defined on — the `_unindex` variables are never staged. No
patient in this cohort reaches G5, so the level is cleared by the existing
`droplevels()` step rather than printed as an empty row.

These are descriptive Table 1 rows on all three analysis pages. **No model
adjusts for them** — the continuous eGFRs are already in the adjuster sets, and
adding a categorical recoding of the same measurement would be collinear.

Cell sizes in the whole cohort (322) are thin at the bottom of the range:

| Stage | `ckd_stage_cre` | `ckd_stage_cre_cys` |
|---|---|---|
| G1 (≥ 90) | 142 | 130 |
| G2 (60–89) | 156 | 153 |
| G3a (45–59) | 18 | 27 |
| G3b (30–44) | 4 | 9 |
| G4 (15–29) | 2 | 3 |

G3b and G4 fall below the project's default small-cell threshold of 11, and the
site is public. Collapsing the bottom into a single "G3b–G5 (< 45)" row is a
one-line change if that is preferred.

### 2.6 Creatinine switched to the 2022P001091 RPDR pull

`ritu_data.qmd` §4.7 now reads labs from **two** extracts:

| Extract | Supplies |
|---|---|
| `2022P001091_mes98_072826125900866852` (Aug 2026) | **CRE only** — 48,843 rows |
| `RPDR 2026-03-19/mes98_03182614492143216` | HGB, PLT, ALB, GHBA1C — 128,156 rows; plus all `Dia` and `Med` |

The newer pull is a **creatinine-only** extract. It contains no `Dia` or `Med`
file and none of HGB / PLT / ALB / GHBA1C, so it cannot be used for "all labs"
literally — doing that would delete `pre_HGB_45days` and `pre_PLT_45days` and
with them the anemia, thrombocytopenia and composite grade 2/3 AE outcomes.
Creatinine therefore comes from the new pull and everything else from the
original one.

Both extracts cover the same 973 EMPIs over the same 1989–2026 span. The newer
one has ~2,800 more creatinine values, which is what the longitudinal eGFR / CKD
/ AKI outcomes in §2.7 are built from.

The cohort flow is unchanged: 946 → 944 after the no-RPDR-labs exclusion (same
two records, 161 and 352) → 943 after dropping record 911 → **322** Local/Early.
Record 911 still has no creatinine within 45 days of baseline in *either* pull,
so its exclusion and the Section 6.1 note remain correct.

A handful of `pre_CRE_45days` values shift slightly because the newer pull has
denser creatinine sampling; that is why `ckd_incidence` / `ckd_progression`
differ by one patient from a run against the old pull.

### 2.7 Long-horizon kidney outcomes ported from the ICI-CKD project

Ported from [`Kavita ICI CKD/data_management.qmd`](https://tianqi-ouyang.github.io/ICI-CKD-/data_management.html)
into `ritu_data.qmd` §§3.8–3.11 (functions) and §§5.16–5.19 (derivations).
Function bodies are copied unchanged; only the anchor differs:

| ICI-CKD | Ritu project |
|---|---|
| `ICI_Dose_Date` | `baseline_sample_date` (first platinum dose) |
| `pre_CRE_180days` | `pre_CRE_45days` (this pipeline's pre-baseline window) |

Variable names are kept identical to the source project — including the
`eskd_defination_1` spelling and the `_after_ici` suffix — so the two projects
stay directly comparable. Here `_after_ici` means "after first platinum".

Results in the 322 Local/Early patients:

| Variable | n (%) |
|---|---|
| `aki` | 63 (19.6 %) |
| `ckd_incidence` | 15 (4.6 %) |
| `ckd_progression` | 20 (6.2 %) |
| `eskd_defination_1` | 1 (0.3 %) |
| `esrd_kt` (at/before baseline) | 1 (0.3 %) |
| `esrd_kt_after_ici` | 1 (0.3 %) |
| `eskd_composite` | 1 (0.3 %) |
| `ckd_composite` | 35 (10.8 %) |

**Where they appear.** All of them are **Table 1 rows** on all three analysis
pages, appended after the eGFR block:

```
eGFR_CRE_baseline, ckd_stage_baseline, last_cre_days,
aki, ckd_incidence, ckd_progression, eskd_defination_1,
esrd_kt_after_ici, eskd_composite, ckd_composite
```

plus `esrd_kt` in the comorbidity row. `ckd_stage_baseline` is shown as a
labelled factor; note that `get_ckd_stage()` collapses G3a/G3b into a single
stage 3, unlike `ckd_stage_cre` / `ckd_stage_cre_cys` from §2.5, which keep them
apart. One patient has no post-baseline creatinine and shows as *Unknown* for
`eGFR_CRE_baseline` and `ckd_stage_baseline`.

They also keep a dedicated section — 2.2 on the whole-cohort page, 3.2 on the two
carboplatin pages — which carries what Table 1 cannot: median time to each event,
the length of creatinine follow-up, and the cross-tabulation of `aki` against
`AKI_kdigo_90days`.

**Time horizon.** Every other outcome in this project is censored at 90 days.
These are not: a qualifying CKD episode must be *sustained for ≥ 90 days*, so
they run over the whole RPDR creatinine record (days from baseline to last
creatinine: median 445, max 1167). Being Table 1 rows they are **descriptive
only** — they are in no Fine-Gray, cause-specific or Cox model, because mixing a
full-follow-up outcome into the 90-day competing-risks framework needs its own
censoring and competing-event decisions.

**`aki` is not the existing AKI variable.** `AKI_g2_90days` / `AKI_g3_90days` /
`AKI_kdigo_90days` grade against a fixed `pre_CRE_45days` reference inside 90
days. `aki` compares each creatinine against the patient's own preceding 90-day
median (or 2-day low) across the entire record. In the whole cohort they agree on
19 patients, `aki` alone flags 44, `AKI_kdigo_90days` alone flags 14, and only 18
of the 63 `aki` events fall inside 90 days (median time to AKI 144 days). Both
are kept; they answer different questions.

**ICD matching caveat.** `dia_list_icd()` matches with `grepl()` on the joined
code vector, so matching is by substring and `.` acts as a regex wildcard —
`T86.1` also captures `T86.10`–`T86.19`. This is the source project's behaviour,
kept unchanged. Diagnoses come from the original 2026-03-19 `Dia` extract, the
only one available.

### 2.8 Entire-cohort page added alongside the Local/Early pages

The project now runs at **two cohort levels off one pipeline run**:

| Page | Cohort | Reads | N |
|---|---|---|---|
| `ritu_entire.qmd` | all platinum, **all stages** | `ritu_full_master.rds` | **943** |
| `ritu_whole.qmd` | all platinum, Local/Early | `ritu_final_master.rds` | 322 |
| `ritu_carbo.qmd` | carboplatin, Local/Early | `ritu_final_master.rds` | 190 |
| `ritu_carbo_auc3.qmd` | carboplatin AUC ≥ 3, Local/Early | `ritu_final_master.rds` | 143 |

`ritu_data.qmd` §5.20 saves `ritu_full_master.rds` immediately **before** the
stage restriction and `ritu_final_master.rds` immediately **after** it, so the
two levels are identical in every variable definition — same creatinine pull,
same CKD / ESKD / AKI family, same Table 1 rows.

`ritu_entire.qmd` is `ritu_whole.qmd` with three changes:

1. reads `ritu_full_master.rds`
2. **`tumor_stage` restored** to Table 1 and to all 10 adjuster sets — it is
   constant in the Local/Early cohort and dropped there (§2.2), but it varies in
   the entire cohort, so the entire-cohort model specification is the same as the
   Jiaxuan project's
3. separate `docs/tables/entire` and `docs/pptx/entire` output directories

Event counts are roughly 3× the Local/Early ones, which materially changes what
is interpretable:

| Outcome | Entire (943) | Local/Early (322) |
|---|---|---|
| `aki` | 240 (25.5 %) | 63 (19.6 %) |
| `ckd_composite` | 70 (7.4 %) | 26 (8.1 %) |
| `ckd_progression` | 61 (6.5 %) | 20 (6.2 %) |
| `ckd_incidence` | 46 (4.9 %) | 15 (4.6 %) |
| `esrd_kt` (baseline) | 5 (0.5 %) | 1 (0.3 %) |
| `eskd_composite` | 2 (0.2 %) | 1 (0.3 %) |
| **90-day death** | **41** | **2** |

The 90-day death models are the clearest difference: degenerate at 322 (§3),
actually estimable at 943 — though 41 events against ~21 covariates is still
about 2 events per covariate, so read them as exploratory.

### 2.9 Outcome definitions revised; competing-risk models added

**AKI is now restricted to within 1 year** of baseline. `aki_any` keeps the
unrestricted count for reference and is not modelled.

**CKD definitions**, all requiring the state to be sustained ≥ 90 days:

| Variable | Definition | Change |
|---|---|---|
| `ckd_incidence` | eGFR ≥ 60 at baseline → < 60 with a **> 25 %** drop | threshold was > 30 % |
| `ckd_incidence_any` | eGFR ≥ 60 at baseline → < 60, **no percent requirement** | new variable |
| `ckd_progression` | eGFR **< 60 at baseline** *and* a further **> 30 %** drop | percent requirement is new |
| `eskd_defination_1` | eGFR < 10 | unchanged |
| `ckd_composite` | `ckd_incidence` + `ckd_progression` + `eskd_composite` | components changed |

`ckd_composite` uses `ckd_incidence` (the > 25 % definition), not
`ckd_incidence_any`, as specified. The three incidence definitions nest strictly:
every `ckd_incidence` patient is in `ckd_incidence_any`.

**The progression change is the consequential one.** It previously read
"eGFR < 60 at baseline *and* after", which counted anyone with stable CKD as
progressing. Requiring a further > 30 % drop takes it from 61 to 8 in the entire
cohort (20 → 5 Local/Early).

| Outcome | Entire (943) | Local/Early (322) |
|---|---|---|
| `aki` (≤ 1 yr) | 197 (20.9 %) | 51 (15.8 %) |
| `aki_any` | 240 (25.5 %) | 63 (19.6 %) |
| `ckd_incidence` (> 25 %) | 61 (6.5 %) | 21 (6.5 %) |
| `ckd_incidence_any` | 89 (9.4 %) | 37 (11.5 %) |
| `ckd_progression` (> 30 %) | 8 (0.8 %) | 5 (1.6 %) |
| `ckd_composite` | 70 (7.4 %) | 26 (8.1 %) |

**Competing-risk models.** Fine-Gray models for `egfr_diff_pct_per10` predicting
AKI (≤ 365 d) and the CKD composite (any time post-baseline), with **death as the
competing event**, are in §4.12–4.13 of the whole and entire pages and §5.12–5.13
of the two carboplatin pages. Adjusters, predictor scaling and the
Cockcroft-Gault adjustment are identical to the 90-day Fine-Gray models.

Death for these models comes from the **RPDR `Date_Of_Death`**, not `death_90`:
the latter is capped at 90 days and cannot place a death in time over a
multi-year horizon. The two sources disagree slightly — RPDR flags 269 deaths in
the entire cohort against REDCap's 41 within 90 days, and RPDR dates 37 of those
within 90 days versus REDCap's 41 — because REDCap follow-up ended earlier. Of
269 RPDR-deceased, 17 have no date and are censored rather than counted as the
competing event, which slightly understates the competing risk.

Censoring is at the last creatinine (day 0 → last follow-up), built by
`cmp_variable_lt()` (§3.4b of `ritu_data.qmd`); the existing `cmp_variable()` is
hard-wired to `time_to_death_90` and cannot be reused.

Competing-risk breakdown (censored / event / death):

| Outcome | Entire (943) | Local/Early (322) |
|---|---|---|
| AKI ≤ 1 yr | 642 / 197 / 104 | 264 / 51 / 7 |
| CKD composite | 632 / 70 / 241 | 271 / 26 / 25 |

Results for `egfr_diff_pct_per10` (subHR per 10 percentage-point decline):

| Model | Entire (943) | Local/Early (322) |
|---|---|---|
| AKI ≤ 1 yr | 1.31 (1.18, 1.47), p < 0.001 | 1.13 (0.87, 1.46), p = 0.38 |
| CKD composite | 1.14 (1.00, 1.30), p = 0.045 | 1.56 (1.13, 2.14), p = 0.006 |

With 70 and 26 events against ~21 covariates, the CKD composite models are
over-parameterised; the AKI model in the entire cohort (197 events) is the only
one of the four that is comfortably powered.

### 2.10 REDCap export refreshed (2026-10-09)

`ritu_data.qmd` §4.1 now reads
`CystatinCInPatientsO_DATA_2026-10-09_1107.csv` from
`.../Thoracic Oncology Spreadsheets/Ritu_project/`, replacing
`../Data/CystatinCInPatientsO_DATA_2026-04-27_0819.csv`. Same 974 records and
the same `record_id` set; 11 columns changed.

**Three direct identifier columns are new in this export — `name`, `mgh_mrn`,
`dob` — and are dropped at read time** with
`select(-any_of(c("name","mgh_mrn","dob")))`, followed by a `stopifnot()` that
fails the render if any survive. They never reach `final_master`, the saved
`.rds`, any table, or the published site. **Do not remove that select**; the site
is public.

What actually changed, and what it touched:

| Column | Records changed | Used by the pipeline? |
|---|---|---|
| `datediff` | 776 | no |
| `fu_date` | 708 | **yes** → `last_follow_up_date` |
| `pt_status` | 84 | **yes** → `death_90` |
| `followup_data_complete` | 64 | no |
| `age` | 3 | yes |
| `platin_dose` | 2 | yes |
| `late_patient_exclusion` | 2 | no |
| `tumor_stage` | 1 | **yes** → the cohort split |
| `form_1_complete`, `hosp_reason_summary`, `yr_fu` | 1 each | no |

**Follow-up is substantially longer.** 695 records gained follow-up (median
+205 days, max +3785), 13 lost some, 266 unchanged. Median
`days_last_follow_up` went 378 → 539 days in the entire cohort and 427 → 603 in
the Local/Early cohort. Recorded deaths rose from 207 to 272.

**One patient moved from Local/Early to Advanced/Metastatic**, so:

| | Before | After |
|---|---|---|
| Entire cohort | 943 | 943 |
| Local/Early | 323 | **322** |
| Local/Early carboplatin | 190 | **189** |
| Local/Early carbo AUC ≥ 3 | 143 | **142** |

The exclusion flow is otherwise identical (946 → 944 → 943; same two records
without RPDR labs, same record 911 drop).

**`death_90` is unchanged at 41.** It requires `days_last_follow_up <= 90`, and
the extended follow-up landed almost entirely on living patients — a patient who
died early already had `fu_date` equal to their death date.

**The kidney outcomes are unaffected.** They are built from RPDR creatinine and
the RPDR `Date_Of_Death`, neither of which changed: `aki` 197, `ckd_composite`
70, `ckd_incidence` 61, `ckd_progression` 8 in the entire cohort, all identical
to the previous run. The Local/Early counts are also unchanged (51 / 26 / 21 / 5),
so the reclassified patient had none of these events.

The AUC ≥ 3 rank fix (§2.2) was re-verified at N = 142: `steroids` is still
constant (142/142) and the `cancer_type_cat_7` `GI` level is still empty, so both
corrections remain necessary and all three cohorts are full rank.

### 2.11 Competing-risk models with death not limited to 90 days

`cmp_variable()` builds every 90-day outcome's competing-risk pair from
`time_to_death_90`, which is capped at day 90. A patient who dies on day 150 is
therefore **censored at day 90** rather than counted as a competing death. That
is correct for a 90-day analysis, but it means almost all of the deaths in the
cohort never enter the models: 41 of 279 in the entire cohort, 2 of 31 in
Local/Early.

A parallel set of models was added that places death at the day it happened.

**Section 5.19b of `ritu_data.qmd`** builds:

| Variable | Definition |
|---|---|
| `death_days_any` / `death_any` | Days to death, RPDR `Date_Of_Death` when present, otherwise the REDCap last-follow-up date when REDCap `pt_status` is "Dead" (for a dead patient that date *is* the death date) |
| `time_to_death_all` / `death_all` | All-cause death over the full follow-up — the uncapped counterpart of `time_to_death_90` / `death_90` |
| `cmp_time_<outcome>_dall` / `cmp_status_<outcome>_dall` | The ten primary 90-day outcomes' competing-risk pairs, rebuilt with `death_days_any` as the competing event and censoring at the last follow-up |

The ten outcomes are `anemia_2`, `anemia_3`, `throm_2`, `throm_3`, `aki_2`,
`aki_kdigo`, `grade_2_aes`, `grade_3_aes`, `any_hosp`, `side_effect`.

**New page sections** — §12 on `ritu_whole.qmd` and `ritu_entire.qmd`, §17 on
`ritu_carbo.qmd` and `ritu_carbo_auc3.qmd`. Each has the ten Fine-Gray models
(subsections .1–.10) and an all-cause-death Cox model over the full follow-up
(.11), with the same predictor (`egfr_diff_pct_per10`) and the same
Cockcroft-Gault-adjusted adjuster set as the existing Fine-Gray section. Word
downloads are `sec12_*.docx` / `sec17_*.docx`.

#### Why `death_days_any` rather than `death_lt`

§2.7 introduced `death_lt` / `death_days_lt` from the RPDR `Date_Of_Death` alone.
The refreshed REDCap export (§2.10) records 272 deaths in the entire cohort while
RPDR dates 252, and the two sets do not nest: 27 REDCap deaths have no RPDR date
and 7 RPDR-dated deaths are not flagged dead in REDCap. Taking the union gives
**279** dated deaths, about 11 % more competing events than RPDR alone. In
Local/Early the union is 31 (28 RPDR-dated, 3 REDCap-only).

**Sections 4.12 and 4.13 were left on `death_days_lt`** so their published
estimates do not change silently. Aligning them on `death_days_any` would add 27
competing deaths to the entire-cohort AKI and CKD-composite models and 3 to
Local/Early.

#### What the new section does and does not estimate

Events are still ascertained over the first 90 days only — that is the window in
which the adverse-event and hospitalization outcomes were abstracted. Risk time
now runs to the last follow-up, so the **cumulative incidence** of each event is
lower in the new section than in the 90-day section and should not be read as an
incidence estimate. The **subdistribution hazard ratio for the exposure** is the
quantity to compare across the two sections.

Event counts change very little; the competing-death column changes a lot:

| Outcome | Entire cohort, event / death (90-day cap) | Entire cohort, event / death (any time) | Local/Early (90-day cap) | Local/Early (any time) |
|---|---|---|---|---|
| Anemia G2 | 433 / 17 | 435 / 126 | 137 / 1 | 137 / 12 |
| Anemia G3 | 133 / 28 | 135 / 224 | 40 / 1 | 41 / 25 |
| Thrombocytopenia G2 | 131 / 27 | 132 / 231 | 37 / 1 | 38 / 28 |
| Thrombocytopenia G3 | 44 / 32 | 45 / 257 | 9 / 1 | 10 / 29 |
| AKI G2 | 89 / 33 | 90 / 249 | 26 / 2 | 26 / 29 |
| AKI (KDIGO) | 109 / 31 | 110 / 240 | 33 / 1 | 33 / 27 |
| Composite G2 AEs | 485 / 12 | 487 / 109 | 148 / 1 | 148 / 12 |
| Composite G3 AEs | 155 / 24 | 157 / 214 | 43 / 1 | 44 / 25 |
| Overall hospitalization | 265 / 13 | 270 / 168 | 69 / 1 | 73 / 16 |
| Platinum-related hospitalization | 163 / 18 | 164 / 202 | 45 / 1 | 46 / 20 |

Two reasons the event column moves at all:

1. **Tie-breaking.** `cmp_variable()` gives death precedence when the event and
   the death fall on the same day; `cmp_variable_lt()` gives the event
   precedence. Overall hospitalization gains the most (265 → 270 entire,
   69 → 73 Local/Early) because dying during an admission is common.
2. Nothing else. Extending follow-up cannot create events, since ascertainment
   still stops at day 90.

The all-cause-death Cox model is considerably better powered than the 90-day one:
279 vs 41 events in the entire cohort, 31 vs 2 in Local/Early. The Local/Early
90-day death models remain numerically degenerate (§3); §12.11 is the usable
mortality model on that cohort, though at 31 events against ~20 covariate degrees
of freedom it is still over-parameterised.

#### A note on the carboplatin pages' heading numbers

On `ritu_carbo.qmd` and `ritu_carbo_auc3.qmd` the top-level and subsection
numbers have been out of step since the Jiaxuan original: `# 5 Fine-Gray Models`
carries subsections `## 4.1`–`## 4.11`, and `# 6 Cause-Specific Hazard Models`
carries `## 5.1`–`## 5.11`. The new section's prose and its "Event (Sec 5)"
column refer to the **top-level** section 5, i.e. Fine-Gray — so the right
comparison for §17.4 is the subsection labelled 4.4, not 5.4 (which is the
cause-specific model). Renumbering would move roughly a hundred headings and
break every in-page anchor and download link on those two pages, so it was left
alone.

#### Which of the new models are usable

Each new section opens with the event-count table and then an automatic
`callout-important` listing every outcome with fewer than 11 events in that
cohort, because a model with ~20 covariate degrees of freedom and a handful of
events can converge and still be meaningless. The AUC ≥ 3 thrombocytopenia
grade 3 model is the clearest example: one event, and it returns subHR 19.2
(7.60, 48.3) with `sex` = Male at 0.00 (0.00, 0.00). It is printed, flagged and
should not be read.

Two models do not fit at all — carboplatin §17.4 (thrombocytopenia grade 3,
4 events) and §17.5 (AKI grade 2, 9 events), both with
`system is computationally singular`. Their 90-day counterparts (§5.4, §5.5) fail
the same way on the same page, so this is the cohort's sparsity rather than
anything introduced here; `.wrap_tbl()` emits the standard
"Model not estimable in this cohort" callout in place of the table.

Outcomes with at least 11 events, by cohort: **ten of ten plus death** in the
entire cohort (943); **nine of ten plus death** in Local/Early (322) —
thrombocytopenia grade 3 has 10; **eight of ten plus death** in carboplatin (189)
— thrombocytopenia grade 3 (4) and AKI grade 2 (9) fall short and are also the
two that fail to fit; and **six of ten** in AUC ≥ 3 (142), where all-cause death
has 10 events and is flagged too.

### 2.12 Explicitly unchanged

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

| Outcome | Whole (322) | Carbo (189) | AUC ≥ 3 (142) |
|---|---|---|---|
| Composite grade 2 AEs | 148 | 81 | 64 |
| Anemia grade 2 | 137 | 79 | 62 |
| Overall hospitalisation | 69 | 36 | 25 |
| Platinum-related hospitalisation | 45 | 20 | 15 |
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
  | Whole (322) | 2 | 18.04 (0.00, an upper limit hundreds of digits long) |
  | Carboplatin (190) | 1 | 25.19 (0.00, Inf) |
  | AUC ≥ 3 (143) | 0 | NA (NA, NA) |

  None of these rows carries information. The underlying fact — **2/322 (0.6 %)
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
│   ├── ritu_data.qmd        data management → both .rds files below
│   ├── ritu_entire.qmd      entire cohort, all stages (N = 943)
│   ├── ritu_whole.qmd       whole Local/Early cohort (N = 322)
│   ├── ritu_carbo.qmd       carboplatin, Local/Early (N = 189)
│   └── ritu_carbo_auc3.qmd  carboplatin AUC ≥ 3, Local/Early (N = 142)
├── ritu_full_master.rds     entire cohort, N = 943 (gitignored — contains PHI)
├── ritu_final_master.rds    Local/Early cohort, N = 322 (gitignored — contains PHI)
└── docs/                    rendered site, committed and served by GitHub Pages
    ├── index.html
    ├── qmd/*.html
    ├── tables/{entire,whole,carbo,auc3}/*.docx
    └── pptx/{entire,whole,carbo,auc3}/*.pptx
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
