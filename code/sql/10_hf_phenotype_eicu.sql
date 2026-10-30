-- =============================================================================
-- 10_hf_phenotype_eicu.sql —— eICU-CRD 2.0 急性心衰表型（方案 2.6 / 表 5）
--
-- 与 MIMIC 侧的对应关系：MIMIC 用 diagnoses_icd 的 seq_num 区分诊断位次，
-- 本库没有位次概念，`diagnosis` 表也不带 seq。因此：
--   · 方案表 5 的"主诊断（seq_num = 1）"在本库**无法直接对应**，
--     以"该次 ICU 单元存在急性心衰特异子码"作为等价的收紧条件。
--   · 宽口径（存在任意 I50.x / 428.x）仍可照搬，用于参照。
--
-- 两条互补的识别路径（方案 2.6），取并集：
--   ① admissionDx.admitdxname = 'CHF, congestive heart failure'
--      —— 入科问诊记录，覆盖较窄但特异性高；
--   ② diagnosis.icd9code 匹配心衰编码 —— 覆盖更广。
--
-- ⚠️ icd9code 字段**同时混装 ICD-9 与 ICD-10**，形如 '428.0, I50.9'
--    （ICD-10 一律在逗号之后）。因此必须双模式匹配：
--      ICD-9  以 '428' 开头（ICD-9 心衰大类）
--      ICD-10 文本中包含 'I50'（心衰大类）
--    只写其中一个会漏掉另一版本编码的患者。
-- ⚠️ ICD-9 代码在本库存为**带小数点**的 '428.0'（与 MIMIC-IV 存 '4280' 不同），
--    故用 LIKE '428%' 前缀匹配，不做等值比较。
-- =============================================================================

-- ---- 路径①：入科问诊记录 ----
CREATE OR REPLACE TABLE hf_admitdx AS
SELECT DISTINCT patientunitstayid
FROM admissionDx
WHERE admitdxname = 'CHF, congestive heart failure';

-- ---- 诊断码拆行：一行一码 ----
-- ⚠️ 这一步是**必须的**，不是化简。diagnosis.icd9code 是一格多码的字符串
--    （形如 '250.00, E11.9, 428.21, I50.31'），直接在该字段上写前缀通配符
--    （'428.21%'）只能匹配到**排在第一位**的那个码——排在逗号之后的同码会被
--    静默漏掉，既不报错也不产生缺失值，只是把一批真患者算成"没有心衰/没有糖尿病"。
--    拆成一行一码后，LIKE 前缀匹配才真正等价于"该码存在"，且 ICD-9 与 ICD-10
--    不再需要两套位置规则。后续所有诊断码匹配（含协变量层的合并症）统一走本表。
CREATE OR REPLACE TABLE dx_code AS
SELECT DISTINCT
  d.patientunitstayid,
  TRIM(unnest(string_split(d.icd9code, ','))) AS code
FROM diagnosis d
WHERE d.icd9code IS NOT NULL
  AND TRIM(d.icd9code) <> '';

-- ---- 路径②：诊断编码（ICD-9 与 ICD-10 双模式） ----
-- 急性心衰特异子码：与 MIMIC 侧同一份清单。
--   ICD-10  I50.21 I50.23 I50.31 I50.33 I50.41 I50.43 I50.811 I50.813
--   ICD-9   428.21 428.23 428.31 428.33 428.41 428.43
-- 每个码都写"带点"与"不带点"两种形式（本库实测为带点，MIMIC-IV 为不带点），
-- 两库口径一致，也免于将来数据版本变化。
CREATE OR REPLACE TABLE hf_dx_flags AS
SELECT
  patientunitstayid,
  MAX(CASE WHEN code LIKE '428%' OR code LIKE 'I50%'
           THEN 1 ELSE 0 END) AS hf_any,
  MAX(CASE WHEN code LIKE '428.21%' OR code LIKE '42821%'
             OR code LIKE '428.23%' OR code LIKE '42823%'
             OR code LIKE '428.31%' OR code LIKE '42831%'
             OR code LIKE '428.33%' OR code LIKE '42833%'
             OR code LIKE '428.41%' OR code LIKE '42841%'
             OR code LIKE '428.43%' OR code LIKE '42843%'
             OR code LIKE 'I50.21%' OR code LIKE 'I50.23%'
             OR code LIKE 'I50.31%' OR code LIKE 'I50.33%'
             OR code LIKE 'I50.41%' OR code LIKE 'I50.43%'
             OR code LIKE 'I50.811%' OR code LIKE 'I50.813%'
           THEN 1 ELSE 0 END) AS acute_hf
FROM dx_code
GROUP BY patientunitstayid;

-- ---- 心衰并集人群（方案 2.6）----
-- hf = 路径① 或 路径②的宽口径命中
CREATE OR REPLACE TABLE hf_union AS
SELECT
  p.patientunitstayid,
  COALESCE(a.patientunitstayid IS NOT NULL, FALSE)              AS by_admitdx,
  COALESCE(d.hf_any, 0)                                         AS by_icd_hf,
  COALESCE(d.acute_hf, 0)                                       AS acute_hf
FROM patient p
LEFT JOIN hf_admitdx  a ON a.patientunitstayid = p.patientunitstayid
LEFT JOIN hf_dx_flags d ON d.patientunitstayid = p.patientunitstayid
WHERE a.patientunitstayid IS NOT NULL OR d.hf_any = 1;

-- ---- 排除标准①：终末期肾病 / 透析依赖（方案 3.4①）----
-- eICU 的 diagnosis 同样双版本混装：ICD-10 N18.6、ICD-9 585.6；
-- 透析依赖状态码 ICD-10 Z99.2、ICD-9 V45.11 / V45.12。
-- ⚠️ 记忆中的勘察结论指出 V45.1 / V56 / Z99.2 / Z49 在本库**全库 0 行**，
--    故本表大概率匹配不到人；保留查询是为了留下"确实查过"的记录，
--    并在论文中说明 eICU 侧 ESRD 的排除主要依赖治疗记录（treatment 表）。
CREATE OR REPLACE TABLE excl_esrd AS
SELECT DISTINCT patientunitstayid
FROM dx_code
WHERE code LIKE '585.6%'  OR code LIKE '5856%'
   OR code LIKE 'N18.6%'
   OR code LIKE 'Z99.2%'
   OR code LIKE 'V45.11%' OR code LIKE 'V45.12%'
   OR code LIKE 'V4511%'  OR code LIKE 'V4512%';
