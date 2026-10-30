-- =============================================================================
-- 00_schema_eicu.sql —— eICU-CRD 2.0 视图层
--
-- 由 R 侧（DBI + duckdb）执行；执行前把 {EICU_ROOT} 替换为数据根目录。
-- 直接读 .csv.gz，不落地转换、不建中间文件。
--
-- ⚠️ 与 MIMIC 侧的两点差异，都是本库数据本身的问题：
--    1) pastHistory、respiratoryCharting 两个文件的 CSV 含未转义的引号，
--       必须显式 quote = '"'，否则 DuckDB 报 unterminated quote。
--    2) 本库 CSV 的时间字段一律是"相对 ICU 入科分钟的整数偏移"（*offset），
--       没有绝对时间戳。所有时点运算都在偏移量上做，不再做时区/绝对时间转换。
--
-- ⚠️ 刻意**不建** apacheApsVar 与 apachePredVar 的视图：这两张表含
--    sodium / creatinine / dialysis / urine 等列，就是暴露与结局本身，
--    一旦建了视图就难免有人在建模时顺手 join 进来。需要用它们时单独读。
-- =============================================================================

CREATE OR REPLACE VIEW patient AS
  SELECT * FROM read_csv_auto('{EICU_ROOT}/patient.csv.gz', header = true, escape = '"');

CREATE OR REPLACE VIEW hospital AS
  SELECT * FROM read_csv_auto('{EICU_ROOT}/hospital.csv.gz', header = true, escape = '"');

CREATE OR REPLACE VIEW lab AS
  SELECT * FROM read_csv_auto('{EICU_ROOT}/lab.csv.gz', header = true, escape = '"');

CREATE OR REPLACE VIEW diagnosis AS
  SELECT * FROM read_csv_auto('{EICU_ROOT}/diagnosis.csv.gz', header = true, escape = '"');

CREATE OR REPLACE VIEW admissionDx AS
  SELECT * FROM read_csv_auto('{EICU_ROOT}/admissionDx.csv.gz', header = true, escape = '"');

CREATE OR REPLACE VIEW intakeOutput AS
  SELECT * FROM read_csv_auto('{EICU_ROOT}/intakeOutput.csv.gz', header = true, escape = '"');

CREATE OR REPLACE VIEW nurseCharting AS
  SELECT * FROM read_csv_auto('{EICU_ROOT}/nurseCharting.csv.gz', header = true, escape = '"');

CREATE OR REPLACE VIEW vitalAperiodic AS
  SELECT * FROM read_csv_auto('{EICU_ROOT}/vitalAperiodic.csv.gz', header = true, escape = '"');

CREATE OR REPLACE VIEW vitalPeriodic AS
  SELECT * FROM read_csv_auto('{EICU_ROOT}/vitalPeriodic.csv.gz', header = true, escape = '"');

CREATE OR REPLACE VIEW medication AS
  SELECT * FROM read_csv_auto('{EICU_ROOT}/medication.csv.gz', header = true, escape = '"');

CREATE OR REPLACE VIEW infusionDrug AS
  SELECT * FROM read_csv_auto('{EICU_ROOT}/infusionDrug.csv.gz', header = true, escape = '"');

CREATE OR REPLACE VIEW treatment AS
  SELECT * FROM read_csv_auto('{EICU_ROOT}/treatment.csv.gz', header = true, escape = '"');

CREATE OR REPLACE VIEW apachePatientResult AS
  SELECT * FROM read_csv_auto('{EICU_ROOT}/apachePatientResult.csv.gz', header = true, escape = '"');

-- 有创通气判定用（cov_vent）：airwaytype 与人机连接时长 ventstartoffset / ventendoffset。
-- 该表不含任何暴露或结局字段，建全表视图无风险。注意与 respiratoryCharting 不是同一张表，
-- 后者 CSV 含未转义引号、且无 airwaytype，不用于本课题。
CREATE OR REPLACE VIEW respiratoryCare AS
  SELECT * FROM read_csv_auto('{EICU_ROOT}/respiratoryCare.csv.gz', header = true, escape = '"');

-- ⚠️ apacheApsVar 的唯一一个例外视图，**只放三列 GCS 分项**。
--    整张 apacheApsVar 含 urine / sodium / creatinine / dialysis —— 就是本研究的
--    结局与暴露本身。故此处**显式只投影 GCS 三列**，让"误把暴露当协变量"在
--    类型层面就不可能发生；需要其它列时不要扩这张视图，另开一个看得到列名的查询。
--    缺失以 -1 哨兵编码（不是 NULL），取用时必须显式过滤。
--
--    用途定位：本视图在 cov_gcs 中只是**兜底**。24 小时最低 GCS 的主源是
--    nurseCharting 的 'GCS Total'（3,248,821 行、取值范围 3–15，覆盖约为本表的
--    十倍，且是真正的窗内最低值而非入科单次值）。两路取 MIN 合并，本表补上
--    nurseCharting 缺失的那部分单元。
CREATE OR REPLACE VIEW apache_gcs AS
  SELECT patientunitstayid, eyes, motor, verbal
  FROM read_csv_auto('{EICU_ROOT}/apacheApsVar.csv.gz', header = true, escape = '"');
