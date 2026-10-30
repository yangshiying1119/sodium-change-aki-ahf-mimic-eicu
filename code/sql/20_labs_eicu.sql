-- =============================================================================
-- 20_labs_eicu.sql —— 血钠 / 肌酐长表（eICU-CRD 2.0 侧提取适配层）
--
-- 输出与 MIMIC 侧完全同构的长表 (id, t_hours, value)：
--   id      = patientunitstayid（字符）
--   t_hours = 相对 ICU 入科时刻的小时数，**可为负**（入科前、住院期间的实测值）
--   value   = 已按方案 4.1.2 边界清洗后的数值
--
-- 与 MIMIC 侧的两点差异：
--   (1) 本库没有 itemid 体系，检验项目由文本名称 labname 标识。**精确匹配**
--       hospital 侧的标准名，不用 LIKE 模糊匹配——否则会把 urinary sodium /
--       urinary creatinine 一并卷进来（方案表 4 已列明须排除）。
--   (2) 单位判据取 labmeasurenamesystem（本库取值唯一：mmol/L、mg/dL）。
--       labmeasurenameinterface 有 mEq/L、MEQ/L、meq/l、MMOLL 等多种写法，
--       不作判据（方案表 4 的说明）。
--
-- 时点一律用 labresultoffset（相对入科分钟）。窗口 [-7 天, +7 天]，
-- 与 MIMIC 侧一致，左端容纳界标前比较基线（方案 4.2.2）。
-- =============================================================================

-- ---- 血钠 ----
CREATE OR REPLACE TABLE sodium AS
SELECT
  CAST(patientunitstayid AS VARCHAR)   AS id,
  labresultoffset / 60.0               AS t_hours,
  labresult                            AS value
FROM lab
WHERE labname = 'sodium'
  AND labmeasurenamesystem = 'mmol/L'
  AND labresult IS NOT NULL
  AND labresultoffset BETWEEN -10080 AND 10080
  -- 方案 4.1.2：血钠生理范围 90–200 mmol/L
  AND labresult BETWEEN 90 AND 200;

-- ---- 肌酐 ----
CREATE OR REPLACE TABLE creatinine AS
SELECT
  CAST(patientunitstayid AS VARCHAR)   AS id,
  labresultoffset / 60.0               AS t_hours,
  labresult                            AS value
FROM lab
WHERE labname = 'creatinine'
  AND labmeasurenamesystem = 'mg/dL'
  AND labresult IS NOT NULL
  AND labresultoffset BETWEEN -10080 AND 10080
  -- 方案 4.2.2：0.1–50 mg/dL
  AND labresult BETWEEN 0.1 AND 50;

-- ---- 血糖（仅供 S5 的血糖校正血钠） ----
-- 只用于把血钠换算成校正血钠，**不进任何模型**；血糖本身作为协变量由协变量层
-- 从 lab 另取最差值。窗口比 MIMIC 侧窄（[-1 天, +4 天]），
-- 因为本表只需覆盖 0–24h 暴露窗附近的配对。
--
-- ⚠️ 必须两个名字都取：本库血糖分 'glucose'（中心实验室，131.9 万行）与
--    'bedside glucose'（床旁快速，317.6 万行）两套，**床旁那套量更大、且多在
--    床边与血钠同批抽**，只取前者会让 S5 的配对率大幅下降。两者单位同为 mg/dL
--    （已实测），可直接合并。同表另有 'glucose - CSF'（脑脊液，全库仅 4,071 行），
--    用**精确匹配**即可天然排除——写成 LIKE '%glucose%' 就会把它卷进来。
CREATE OR REPLACE TABLE glucose AS
SELECT
  CAST(patientunitstayid AS VARCHAR)   AS id,
  labresultoffset / 60.0               AS t_hours,
  labresult                            AS value
FROM lab
WHERE labname IN ('glucose', 'bedside glucose')
  AND labmeasurenamesystem = 'mg/dL'
  AND labresult IS NOT NULL
  AND labresultoffset BETWEEN -1440 AND 5760
  AND labresult BETWEEN 20 AND 1000;

-- ---- 钾 / 氯 —— 仅作阴性对照暴露，不进任何主分析或主模型 ----
-- 与 MIMIC 侧同名同构，口径一致（0–24h 末次减首次，基础模型）。同理，两库共用
-- common.R 的 delta_na() 构造暴露，只有提取适配层不同。
-- 范围：钾 1.0–10.0 mmol/L、氯 70–150 mmol/L。
CREATE OR REPLACE TABLE potassium AS
SELECT
  CAST(patientunitstayid AS VARCHAR)   AS id,
  labresultoffset / 60.0               AS t_hours,
  labresult                            AS value
FROM lab
WHERE labname = 'potassium'
  AND labmeasurenamesystem = 'mmol/L'
  AND labresult IS NOT NULL
  AND labresultoffset BETWEEN -10080 AND 10080
  AND labresult BETWEEN 1.0 AND 10.0;

CREATE OR REPLACE TABLE chloride AS
SELECT
  CAST(patientunitstayid AS VARCHAR)   AS id,
  labresultoffset / 60.0               AS t_hours,
  labresult                            AS value
FROM lab
WHERE labname = 'chloride'
  AND labmeasurenamesystem = 'mmol/L'
  AND labresult IS NOT NULL
  AND labresultoffset BETWEEN -10080 AND 10080
  AND labresult BETWEEN 70 AND 150;
