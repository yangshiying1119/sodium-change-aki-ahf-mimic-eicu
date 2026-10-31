# Early change in serum sodium and subsequent acute kidney injury in acute heart failure: analysis code

Landmark-cohort analysis of two independent critical-care databases — **MIMIC-IV (v3.1)** and **eICU-CRD (v2.0)** — testing whether the direction of serum sodium change during the first 24 ICU hours is associated with KDIGO-defined AKI over hours 24–72 in patients with acute heart failure. Full methods are reported in the manuscript (see Citation below).

## Repository structure

```text
tests/
  run_tests.R           unit-test entry point (69 assertions; run before use)

code/
  sql/     server-side extraction layer (DuckDB 1.2.0 dialect)
           00_schema_*.sql       database schema setup
           10_hf_phenotype_*.sql heart-failure phenotype flags (eICU; the MIMIC-IV equivalent is embedded in 30_cohort_mimic.sql)
           20_labs_*.sql         chemistry laboratory series (sodium, creatinine, glucose)
           21_urine_*.sql        urine-output series (sensitivity outcome)
           30_cohort_*.sql       landmark cohort bases + phenotype tiers A/B/C
           40–49_*.sql           covariates, RRT, Charlson comorbidity

  r/       analysis layer (R)
           cohort.R / common.R  shared landmark-eligibility logic (identical across databases)
           test_common.R        unit tests (69 assertions)
           extract.R            lab-series assembly (set MIMIC_ROOT to your local data path)
           models.R / bias.R / sensitivity.R / competing.R / discrim.R / advance.R
                                 three-tier logistic models, IPW / MI / negative controls,
                                 Fine–Gray competing risk, c-statistic contrasts
           tables.R / figures.R / figures_cross.R / export.R / i18n.R
                                 output tables and figures
           _run*.R              staged pipeline entry points (_runE* target eICU)
           _diag*.R / _probe*.R / _review*.R / _smoke*.R / _sqlcheck*.R / _eicu_recon*.R
                                 development diagnostics and one-off checks, retained for transparency only
```

The same phenotype and laboratory logic is written once per database pair and executed unchanged on each database; MIMIC-IV and eICU-CRD differ only in the SQL extraction layer, which is deliberately structured so that downstream R code is database-agnostic.

The primary cohorts use database-specific phenotype tiers — MIMIC-IV tier B and eICU-CRD tier C — because eICU charts predominantly unspecified heart-failure codes (I50.9 / 428.0), which leaves the strict tier B underpowered in that database. The tier definitions and this choice are documented in `eTable_phenotype_codes.md` and in the manuscript.

## Requirements

- **R ≥ 4.1** (analysis environment: R 4.1.2). Packages actually used by the code: `DBI`, `duckdb` (1.2.0), `mice`, `cmprsk`, `pROC`, `MASS`, `jsonlite`, `survival`, `ragg`, `forestploter`, `openxlsx` (`MASS` and `survival` ship with base R installations as recommended packages; the rest must be installed). Restricted cubic splines are generated with `splines::ns` (natural cubic splines, i.e. restricted cubic splines with linear tails) from base R's recommended `splines` package; `rms` is not used. Installers: `code/r/00_install_pkgs.R` (from CRAN, with archive fallback) and `code/r/00b_install_archived.R` (pinned source version of `duckdb` 1.2.0 for R 4.1.x).
- **DuckDB 1.2.0** for the SQL extraction layer (accessed from R through the `duckdb` package).

Note: the IPW sensitivity analysis uses retesting-probability weights implemented directly through `glm(weights = ...)` in `code/r/bias.R`; no dedicated weighting package is required.

The exact session information of the server analysis environment is archived as `sessionInfo.txt` in the repository root (R 4.1.2 under Ubuntu 22.04.5, with the exact loaded versions of the analysis packages, e.g. duckdb 1.2.0, mice 3.19.0, cmprsk 2.2-12, pROC 1.19.0.1, forestploter 1.2.0, ragg 1.5.2, openxlsx 4.2.8.1); the SQL engine version (`duckdb` 1.2.0) is additionally pinned by `code/r/00b_install_archived.R`.

## Data access

No patient-level data are included. Both source databases are de-identified and publicly available to credentialed investigators:

- MIMIC-IV v3.1: https://physionet.org/content/mimiciv/3.1/
- eICU-CRD v2.0: https://physionet.org/content/eicu-crd/2.0/

Access requires completion of the PhysioNet credentialing process and a signed data use agreement. After obtaining the data, set the `MIMIC_ROOT` path in `code/r/extract.R` (and the corresponding eICU path) to your local copy.

## Reproduction outline

1. Do not run the `code/sql/` files manually in filename order. The R layer (`code/r/extract.R`) reads and executes the SQL scripts against the database in the correct dependency order itself (notably `21_urine_*.sql` runs after `30_cohort_*.sql`, not before). After setting `MIMIC_ROOT` (and the eICU path) in `code/r/extract.R`, start from the R entry points below and the extraction layer runs with them.

2. Run the MIMIC-IV stages in dependency order. The `_run*` numbers are historical labels, **not** an execution order — each stage reads `.rds` files written by earlier stages, so a fresh run should follow: `_run2.R` (extraction, landmark cohort, covariates) → `_run3.R` (three-tier main models) → `_run4.R` (sensitivity analyses, competing risk) → `_run7.R` (urine-output composite, bias analyses) → `_run10.R` (discrimination and calibration) → `_run11.R` (supplementary analyses) → `_run13.R` (negative controls) → `_run5.R` / `_run6.R` (figures, tables) → `_run12.R` (supplementary and English tables; reads the outputs of `_run11.R` and `_run13.R`) → `_run9.R` (final assembly of `results.json`; reads the outputs of `_run3/_run4/_run7/_run10` and therefore runs last — executed early, it silently omits results that do not exist yet). There is no `_run8.R`. `_run1.R` is a standalone cohort-assembly and scale check; `_runA.R` is a strict-phenotype (tier A) side analysis that reuses the extraction cache written by `_run2.R`. Neither feeds the main sequence.

3. For eICU-CRD, run in order: `_runE1.R`, `_runE2.R`, `_runE3.R`, `_runE7.R`, then `_runEC.R` (table export), `_runEF.R` (figures), `_runEInline.R` (English inline tables). Each eICU stage selects the phenotype tier through a command-line argument (A/B/C): pass `C` to reproduce the primary eICU cohort (e.g. `Rscript _runE1.R C`). With no argument, `_runE1–E3` default to the strict tier B — the underpowered sensitivity tier — whereas `_runE7`, `_runEC` and `_runEF` default to `C`; pass the tier explicitly to avoid mixing cohorts.

4. The `_run*` scripts form a staged pipeline: each stage persists intermediate `.rds` results that later stages read back (several later stages only re-assemble saved results and never refit models), so run them in order within a single working directory.

5. Sensitivity and bias analyses (IPW, multiple imputation, negative controls, mutual adjustment, competing risk) are driven by the `_run*` entry points and the functions in `bias.R` / `sensitivity.R`.

6. `code/sql/49_charlson_mimic.sql` is not invoked by the extraction entry points; it is retained as an optional reference script for comorbidity-burden sensitivity analyses (the Charlson index is not part of the prespecified covariate set). It is adapted from the official PhysioNet mimic-code repository (`mimic-iv/concepts/comorbidity/charlson.sql`, 2022-08-12 snapshot), with attribution retained in the file header.

7. Unit tests: run `Rscript tests/run_tests.R` from the repository root before anything else. This entry point executes the 69 assertions defined in `code/r/test_common.R` (shared landmark-eligibility and KDIGO logic). A non-zero exit code indicates failure.

Comments inside the scripts are written in Chinese; identifiers and structure are in English.

## Reproducibility verification

The R analysis layer was verified end-to-end in a second environment (Windows 11, R 4.5.3 — a different R version from the original server, R 4.1.2). Starting from the server-produced extraction files (`mimic_extract_B.rds`, `eicu_B_extract.rds`; the SQL extraction layer itself requires the raw databases and was not re-run locally), the full staged pipeline — cohorts, three-tier models, sensitivity analyses, competing risk, IPW/multiple imputation, negative controls, mutual adjustment, discrimination/calibration, supplementary analyses, tables, figures, and the final `results.json` assembly — was re-executed for both databases, the strict-phenotype side analysis, and the two-database contrast. Leaf-by-leaf comparison of the regenerated `results.json` against the archived server output showed no substantive differences: the only discrepancies were (i) row order within the "ascertainability grouping" table (category–value pairs identical in all 6 rows, both databases) and (ii) floating-point differences below 1e-12 relative magnitude in a handful of ordered-logistic and standardised-probability values. All reported estimates (ORs, CIs, P values, cohort counts) agree digit-for-digit at reported precision.

## Phenotype code lists

The full ICD code lists defining the acute heart failure phenotype (broad and acute-specific codes, exclusion codes) and the A/B/C phenotype tier definitions per database are provided in `eTable_phenotype_codes.md` and as a supplementary table in the manuscript.

## Citation

Manuscript currently under review:

> Early change in serum sodium and subsequent acute kidney injury in acute heart failure: a landmark cohort analysis of MIMIC-IV and eICU. *Under review.*

The final bibliographic citation and article DOI will be added upon publication.

## Software citation

The research code is archived on Zenodo:

Shiying Yang. (2026). *Early change in serum sodium and subsequent acute kidney injury in acute heart failure: analysis code* (v1.0.1). Zenodo.

Version-specific DOI (v1.0.1):
https://doi.org/10.5281/zenodo.23830227

Concept DOI (all versions):
https://doi.org/10.5281/zenodo.23812108

GitHub repository:

https://github.com/yangshiying1119/sodium-change-aki-ahf-mimic-eicu

## License

MIT for code (see LICENSE).

## Contact

Shiying Yang — yxs18221189048@163.com · ORCID 0009-0002-4630-7553

This work received no specific funding. Both databases are used under PhysioNet data use agreements; MIMIC-IV was approved for research release by the Institutional Review Board of Beth Israel Deaconess Medical Center with a waiver of informed consent, and eICU-CRD was released under equivalent agreements from its contributing institutions. The ethics statement applicable to these analyses is reported in the manuscript.
