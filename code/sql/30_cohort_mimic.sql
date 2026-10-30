-- =============================================================================
-- 30_cohort_mimic.sql —— MIMIC-IV 侧界标队列与表型口径
--
-- 节点编号与方案表 6 一致（N0–N8）：
--   N0 数据库内全部 ICU 住院        N1 成人（≥18 岁）
--   N2 急性心衰表型                 N3 首次符合条件的 ICU 住院   ← 在 R 侧完成
--   N4 界标时仍在 ICU 且存活        N5 0–24h ≥2 次血钠            ← 在 R 侧完成
--   N6 界标前可判定为无 AKI         N7 界标后有可判定比较值      ← 在 R 侧完成
--   N8 结局窗内检出 AKI
-- 本文件负责 N1、N2 与 N4 的判定依据，以及排除标准①；其余需化验序列，在 R 侧算。
--
-- 表型定义（方案表 5；主分析原定 A，已修订为 B，见方案修订记录）：
--   A  主诊断（seq_num = 1）本身即急性心衰码
--   B  前 3 位诊断含 I50.x / 428.x **且** 存在急性心衰子码（不限位次）  ← 主分析
--   C  任意位次含 I50.x / 428.x                                       ← 仅参照
--
-- 急性心衰子码（ICD-10 8 个 + ICD-9 6 个）：
--   ICD-10  I5021 I5023 I5031 I5033 I5041 I5043 I50811 I50813
--   ICD-9   42821 42823 42831 42833 42841 42843
-- 注意：MIMIC-IV 的 ICD-9 码**无小数点**，存为 '42821' 而非 '428.21'，直接等值比较。
-- 慢性心衰（I50.0/428.0 等）不得单独构成急性表型。
-- =============================================================================

-- ---- N1 基础层：成人（方案 3.2①）。停留时长与死亡留到 N4 判，不在此处筛 ----
CREATE OR REPLACE TABLE base AS
SELECT
  ie.subject_id,
  ie.hadm_id,
  ie.stay_id,
  ie.intime,
  ie.outtime,
  ie.first_careunit,
  date_diff('minute', ie.intime, ie.outtime) / 60.0 AS los_hours,
  p.anchor_age + (EXTRACT(YEAR FROM ie.intime) - p.anchor_year) AS age,
  p.gender,
  adm.admission_type,
  adm.insurance,
  adm.race,
  adm.deathtime,
  adm.hospital_expire_flag
FROM icustays ie
JOIN patients p
  ON p.subject_id = ie.subject_id
JOIN admissions adm
  ON adm.hadm_id = ie.hadm_id
WHERE ie.intime IS NOT NULL
  AND ie.outtime IS NOT NULL
  AND (p.anchor_age + (EXTRACT(YEAR FROM ie.intime) - p.anchor_year)) >= 18;

-- ---- N4 判定依据：界标（t0+24h）时仍在 ICU 且存活（方案 3.2④） ----
-- 排除原因拆成"24h 内死亡"与"24h 内转出"两类，供流程图分列。
CREATE OR REPLACE TABLE landmark_flags AS
SELECT
  stay_id,
  CASE WHEN outtime > intime + INTERVAL '24 hours' THEN 1 ELSE 0 END AS still_in_icu,
  CASE WHEN deathtime IS NOT NULL AND deathtime <= intime + INTERVAL '24 hours'
       THEN 1 ELSE 0 END AS died_within_24h
FROM base;

-- ---- 排除标准①：透析依赖或终末期肾病（方案 3.4①） ----
-- 两条路径：(a) 直接诊断为终末期肾病（ICD-10 N186 / ICD-9 5856）；
--           (b) 既往诊断显示透析依赖（ICD-10 Z992；ICD-9 V4511 / V4512）。
-- 理由：此类患者肌酐无法按 KDIGO 解读。按 hadm_id 粒度排除，整个住院都不进分析。
CREATE OR REPLACE TABLE excl_esrd AS
SELECT DISTINCT hadm_id
FROM diagnoses_icd
WHERE icd_code IN ('N186', '5856', 'Z992', 'V4511', 'V4512');

-- ---- N2 心衰表型标志 ----
CREATE OR REPLACE TABLE hf_flags AS
SELECT
  b.stay_id,
  MAX(CASE WHEN d.seq_num = 1
            AND (d.icd_code LIKE 'I50%' OR d.icd_code LIKE '428%')
           THEN 1 ELSE 0 END) AS hf_primary,
  MAX(CASE WHEN d.seq_num <= 3
            AND (d.icd_code LIKE 'I50%' OR d.icd_code LIKE '428%')
           THEN 1 ELSE 0 END) AS hf_top3,
  MAX(CASE WHEN d.icd_code LIKE 'I50%' OR d.icd_code LIKE '428%'
           THEN 1 ELSE 0 END) AS hf_any,
  MAX(CASE WHEN d.icd_code IN (
        'I5021','I5023','I5031','I5033','I5041','I5043','I50811','I50813',
        '42821','42823','42831','42833','42841','42843')
       THEN 1 ELSE 0 END) AS acute_hf,
  -- A 口径要求主诊断**本身就是急性心衰码**，而非"主诊断是宽泛 I50/428 且别处另有急性码"
  MAX(CASE WHEN d.seq_num = 1 AND d.icd_code IN (
        'I5021','I5023','I5031','I5033','I5041','I5043','I50811','I50813',
        '42821','42823','42831','42833','42841','42843')
       THEN 1 ELSE 0 END) AS acute_primary
FROM base b
JOIN diagnoses_icd d
  ON d.hadm_id = b.hadm_id
GROUP BY b.stay_id;

-- ---- N2 三个表型口径的队列（带 N4 与排除①标志，供 R 侧逐层筛） ----
CREATE OR REPLACE TABLE cohort_A AS          -- 严格表型，敏感性分析用
SELECT b.*,
       lf.still_in_icu,
       lf.died_within_24h,
       CASE WHEN e.hadm_id IS NULL THEN 0 ELSE 1 END AS esrd_excluded
FROM base b
JOIN hf_flags f USING (stay_id)
JOIN landmark_flags lf USING (stay_id)
LEFT JOIN excl_esrd e ON e.hadm_id = b.hadm_id
WHERE f.acute_primary = 1;

CREATE OR REPLACE TABLE cohort_B AS          -- 主分析
SELECT b.*,
       lf.still_in_icu,
       lf.died_within_24h,
       CASE WHEN e.hadm_id IS NULL THEN 0 ELSE 1 END AS esrd_excluded
FROM base b
JOIN hf_flags f USING (stay_id)
JOIN landmark_flags lf USING (stay_id)
LEFT JOIN excl_esrd e ON e.hadm_id = b.hadm_id
WHERE f.hf_top3 = 1 AND f.acute_hf = 1;

CREATE OR REPLACE TABLE cohort_C AS          -- 仅参照
SELECT b.*,
       lf.still_in_icu,
       lf.died_within_24h,
       CASE WHEN e.hadm_id IS NULL THEN 0 ELSE 1 END AS esrd_excluded
FROM base b
JOIN hf_flags f USING (stay_id)
JOIN landmark_flags lf USING (stay_id)
LEFT JOIN excl_esrd e ON e.hadm_id = b.hadm_id
WHERE f.hf_any = 1;
