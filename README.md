# Ritu Project — Local/Early Stage Cohort

Creatinine-to-cystatin C discordance and adverse events in patients receiving
platinum-based chemotherapy, restricted to **locally staged / early-stage
disease** (`tumor_stage == "Local/Early"`, stage I–II).

**Live site:** https://tianqi-ouyang.github.io/Ritu-s-project-/

This is a derivative of the [Jiaxuan
project](https://tianqi-ouyang.github.io/Cys-C-project-/), with two substantive
differences:

1. **Cohort** — all 621 advanced/metastatic patients are excluded from the
   Local/Early pages, leaving **N = 322**. An entire-cohort page (N = 943) runs
   the same analysis unrestricted.
2. **Modelled outcomes** — only **AKI within 1 year** and the **CKD composite**
   are modelled, with Fine-Gray competing risks and death at any time as the
   competing event. The 90-day adverse-event, hospitalization and 90-day-death
   models are no longer on the pages; those outcomes are still derived by the
   pipeline and kept in the saved master. See `RITU_PROJECT.md` §2.12.

## Cohorts

| Page | Cohort | Jiaxuan (all stages) | Ritu (Local/Early) |
|---|---|---|---|
| Whole Cohort | all platinum | 943 | **322** |
| Carboplatin Cohort | carboplatin | 596 | **189** |
| Carboplatin AUC ≥ 3 | carboplatin, AUC ≥ 3 | 463 | **142** |

## Contents

| File | Purpose |
|---|---|
| `RITU_PROJECT.md` | Full project notes — cohort definition, what changed, power caveats, paths |
| `index.qmd` | Site summary page |
| `qmd/ritu_data.qmd` | Data management pipeline → `ritu_final_master.rds` |
| `qmd/ritu_whole.qmd` | Whole Local/Early cohort analysis |
| `qmd/ritu_carbo.qmd` | Carboplatin cohort analysis |
| `qmd/ritu_carbo_auc3.qmd` | Carboplatin AUC ≥ 3 subgroup analysis |
| `docs/` | Rendered website (served by GitHub Pages) |

## Rendering

```bash
quarto render
```

`ritu_data.qmd` renders first and writes the `ritu_final_master.rds` that the
three analysis pages read. It re-reads the full RPDR pull and takes several
minutes; the analysis pages render in well under a minute each.

The project expects to sit in a `Ritu project/` folder alongside the parent
`Cystain C redcap/` working directory, from which it reads `../Data/`. RPDR
`.txt` inputs are read from their original absolute paths. See the *Path
conventions* section of `RITU_PROJECT.md`.

## Data privacy

No patient-level data is committed. `.gitignore` excludes `*.rds`, `*.csv`,
`*.xlsx`, `*.pdf` and the `Data/` folder. The rendered site contains only
aggregate tables, model coefficients and figures — no row-level records and no
identifier columns.
