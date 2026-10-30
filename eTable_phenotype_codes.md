# eTable (Supplementary). Phenotype code lists and definitions

Code lists verified verbatim against the extraction SQL layer (`code/sql/10_hf_phenotype_*.sql`, `code/sql/30_cohort_*.sql`).

## Panel A. Acute-specific heart failure ICD codes (constituting the acute phenotype)

| ICD-10 code | Description | ICD-9 code | Description |
|---|---|---|---|
| I50.21 | Acute systolic heart failure | 428.21 | Acute systolic heart failure |
| I50.23 | Acute on chronic systolic heart failure | 428.23 | Acute on chronic systolic heart failure |
| I50.31 | Acute diastolic heart failure | 428.31 | Acute diastolic heart failure |
| I50.33 | Acute on chronic diastolic heart failure | 428.33 | Acute on chronic diastolic heart failure |
| I50.41 | Acute combined systolic and diastolic heart failure | 428.41 | Acute combined systolic and diastolic heart failure |
| I50.43 | Acute on chronic combined systolic and diastolic heart failure | 428.43 | Acute on chronic combined systolic and diastolic heart failure |
| I50.811 | Acute right heart failure | | |
| I50.813 | Acute on chronic right heart failure | | |

## Panel B. Broad heart-failure codes and the eICU admission-diagnosis path

| Element | Definition |
|---|---|
| Broad HF codes, both databases | Any ICD-10 code starting with `I50` or any ICD-9 code starting with `428` (any position) |
| eICU admission-diagnosis path | `admissionDx.admitdxname = 'CHF, congestive heart failure'` |
| eICU HF union population | Admission-diagnosis path or broad-code path (either) |

## Panel C. Exclusion codes: end-stage kidney disease or dialysis dependence

| Code system | Code | Description |
|---|---|---|
| ICD-10 | N18.6 | End-stage renal disease |
| ICD-9 | 585.6 | End-stage renal disease |
| ICD-10 | Z99.2 | Dependence on renal dialysis (status) |
| ICD-9 | V45.11 | Renal dialysis status |
| ICD-9 | V45.12 | Post-renal transplant status |

Note: in eICU-CRD these codes match few patients; exclusion there relies mainly on documented renal replacement therapy in the treatment tables.

## Panel D. Phenotype tier definitions, as implemented per database

| Tier | MIMIC-IV (uses `diagnoses_icd.seq_num` position) | eICU-CRD (no diagnosis position recorded) |
|---|---|---|
| A (strictest; sensitivity analysis S6) | Principal diagnosis (`seq_num = 1`) is itself an acute-specific code | Admission diagnosis is CHF and an acute-specific code exists (S6: n = 48, could not be fitted) |
| B (primary, MIMIC-IV) | A heart-failure code (I50.x / 428.x) among the first three diagnosis positions and an acute-specific code in any position | HF union population and an acute-specific code in any position (sensitivity analysis: 210 patients, 44 events) |
| C (broadest) | Any I50.x / 428.x in any position (reference only) | HF union population, acute or chronic (primary, eICU-CRD) |

Storage-format notes, as handled in the extraction code:

- MIMIC-IV stores ICD-9 codes **without** decimal points (e.g. `42821`); matching is by exact string.
- eICU-CRD mixes ICD-9 (dotted, e.g. `428.21`) and ICD-10 (e.g. `I50.21`) in a single comma-separated `diagnosis.icd9code` field; the code splits this field into one row per code and matches both formats, so codes in either position are captured.
- The acute-specific list in Panel A is identical across the two databases; only the position/tier logic differs.
