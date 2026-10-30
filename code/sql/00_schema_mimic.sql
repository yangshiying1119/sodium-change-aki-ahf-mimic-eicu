-- =============================================================================
-- 00_schema_mimic.sql —— MIMIC-IV v3.1 视图层
--
-- 由 R 侧（DBI + duckdb）执行；执行前把 {MIMIC_ROOT} 替换为数据根目录。
-- 直接读 .csv.gz，不落地转换、不建中间文件。
--
-- ⚠️ escape 必须显式写成 '"'。R 侧 duckdb 1.2.0 的默认 escape 判定有误，遇到
--    labevents 里形如 "VERIFIED, ""5 NRBCS""/100 WBCS SEEN." 的合法 RFC4180
--    转义字段会误报 "Value with unterminated quote found."（Python 侧 1.5.5 无此
--    问题）。显式指定后行数与 Python 侧逐条一致（50983 均为 4,111,289 行）。
-- =============================================================================

-- ---- hosp 模块 ----
CREATE OR REPLACE VIEW labevents AS
  SELECT * FROM read_csv_auto('{MIMIC_ROOT}/hosp/labevents.csv.gz', header = true, escape = '"');

CREATE OR REPLACE VIEW patients AS
  SELECT * FROM read_csv_auto('{MIMIC_ROOT}/hosp/patients.csv.gz', header = true, escape = '"');

CREATE OR REPLACE VIEW admissions AS
  SELECT * FROM read_csv_auto('{MIMIC_ROOT}/hosp/admissions.csv.gz', header = true, escape = '"');

CREATE OR REPLACE VIEW diagnoses_icd AS
  SELECT * FROM read_csv_auto('{MIMIC_ROOT}/hosp/diagnoses_icd.csv.gz', header = true, escape = '"');

CREATE OR REPLACE VIEW d_icd_diagnoses AS
  SELECT * FROM read_csv_auto('{MIMIC_ROOT}/hosp/d_icd_diagnoses.csv.gz', header = true, escape = '"');

CREATE OR REPLACE VIEW prescriptions AS
  SELECT * FROM read_csv_auto('{MIMIC_ROOT}/hosp/prescriptions.csv.gz', header = true, escape = '"');

CREATE OR REPLACE VIEW emar AS
  SELECT * FROM read_csv_auto('{MIMIC_ROOT}/hosp/emar.csv.gz', header = true, escape = '"');

-- ---- icu 模块 ----
CREATE OR REPLACE VIEW icustays AS
  SELECT * FROM read_csv_auto('{MIMIC_ROOT}/icu/icustays.csv.gz', header = true, escape = '"');

CREATE OR REPLACE VIEW chartevents AS
  SELECT * FROM read_csv_auto('{MIMIC_ROOT}/icu/chartevents.csv.gz', header = true, escape = '"');

CREATE OR REPLACE VIEW outputevents AS
  SELECT * FROM read_csv_auto('{MIMIC_ROOT}/icu/outputevents.csv.gz', header = true, escape = '"');

CREATE OR REPLACE VIEW inputevents AS
  SELECT * FROM read_csv_auto('{MIMIC_ROOT}/icu/inputevents.csv.gz', header = true, escape = '"');

CREATE OR REPLACE VIEW procedureevents AS
  SELECT * FROM read_csv_auto('{MIMIC_ROOT}/icu/procedureevents.csv.gz', header = true, escape = '"');

CREATE OR REPLACE VIEW d_items AS
  SELECT * FROM read_csv_auto('{MIMIC_ROOT}/icu/d_items.csv.gz', header = true, escape = '"');
