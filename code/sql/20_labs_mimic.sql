-- =============================================================================
-- 20_labs_mimic.sql —— 血钠 / 肌酐长表（MIMIC-IV 侧提取适配层）
--
-- 输出统一长表 (id, t_hours, value)：
--   id      = stay_id（字符）
--   t_hours = 相对 ICU 入科时刻的小时数，**可为负**（院前实测值）
--   value   = 已按方案 4.1.2 边界清洗后的数值
--
-- 两个刻意的设计选择：
--   (1) 不做同时间合并——方案 4.1.2 规定"取中位数、偶数条取较低者"，该规则
--       实现在 r/common.R::merge_same_time()，SQL 侧保持原始行，避免两处各写
--       一套。官方 mimic-code 此处用 AVG，本方案不跟。
--   (2) 取值窗为 [intime − 7 天, intime + 7 天]。左端容纳院前实测值用于构造
--       界标前比较基线（方案 4.2.2）；右端覆盖 72 小时结局窗并留余量。
--       实际参与判定的窗口由 common.R 的 48h/7d 滚动窗控制，不在此处截断。
-- =============================================================================

-- ---- 血钠 ----
-- 50983 血清/血浆生化钠（主源）；52623 同义新条目，低频但同源，合并入主分析。
-- 50824 全血（血气）钠与生化钠不完全等价，单列供敏感性分析，不进主分析。
CREATE OR REPLACE TABLE sodium AS
WITH src AS (
  SELECT
    ie.stay_id,
    ie.intime,
    le.charttime,
    le.valuenum,
    le.valueuom,
    le.itemid
  FROM icustays ie
  JOIN labevents le
    ON le.subject_id = ie.subject_id
   AND le.itemid IN (50983, 52623)
   AND le.valuenum IS NOT NULL
   AND le.charttime BETWEEN ie.intime - INTERVAL 7 DAY AND ie.intime + INTERVAL 7 DAY
)
SELECT
  CAST(stay_id AS VARCHAR)                                   AS id,
  date_diff('minute', intime, charttime) / 60.0              AS t_hours,
  valuenum                                                   AS value
FROM src
-- 方案 4.1.2：血钠生理范围 90–200 mmol/L，超出视为录入错误
WHERE valuenum BETWEEN 90 AND 200
  -- 单位断言：MIMIC-IV 侧 50983/52623 的 valueuom 应为 mEq/L（与 mmol/L 数值等价）
  AND valueuom = 'mEq/L';

-- ---- 全血（血气）钠 —— 仅敏感性分析 ----
CREATE OR REPLACE TABLE sodium_blood AS
SELECT
  CAST(ie.stay_id AS VARCHAR)                          AS id,
  date_diff('minute', ie.intime, le.charttime) / 60.0  AS t_hours,
  le.valuenum                                          AS value
FROM icustays ie
JOIN labevents le
  ON le.subject_id = ie.subject_id
 AND le.itemid = 50824
 AND le.valuenum IS NOT NULL
 AND le.charttime BETWEEN ie.intime - INTERVAL 7 DAY AND ie.intime + INTERVAL 7 DAY
WHERE le.valuenum BETWEEN 90 AND 200;

-- ---- 血糖（仅供 S5 的血糖校正血钠） ----
-- 50931 生化血糖；50809 血气血糖。方案 4.1.2 指出 50809 有量级异常记录
-- （最大 1,276,103 mg/dL），故与协变量层同样限幅 20–1000。
-- 本表**不进任何模型**——血糖本身已是扩展模型的协变量（cov_labs.glucose_max_24h），
-- 这里取的是逐时刻序列，只用于把血钠换算成校正血钠。
CREATE OR REPLACE TABLE glucose AS
SELECT
  CAST(ie.stay_id AS VARCHAR)                          AS id,
  date_diff('minute', ie.intime, le.charttime) / 60.0  AS t_hours,
  le.valuenum                                          AS value
FROM icustays ie
JOIN labevents le
  ON le.subject_id = ie.subject_id
 AND le.itemid IN (50931, 50809)
 AND le.valuenum IS NOT NULL
 AND le.charttime BETWEEN ie.intime - INTERVAL 1 DAY AND ie.intime + INTERVAL 4 DAY
WHERE le.valuenum BETWEEN 20 AND 1000;

-- ---- 肌酐 ----
-- 50912 血清/血浆生化肌酐（主源，mg/dL）；52546 同义新条目，方案表 7 明确要求合并
-- 入主分析。注意 52546 的 hadm_id 为空，但本查询只按 subject_id + charttime 关联，
-- 不受影响。
-- 注意：labevents 中 51081/51082/51106/52024 为尿/体液或全血肌酐，**不是血清肌酐，
-- 不得混入主分析**（52024 全血肌酐按方案表 7 仅作敏感性分析）。
CREATE OR REPLACE TABLE creatinine AS
WITH src AS (
  SELECT
    ie.stay_id,
    ie.intime,
    le.charttime,
    le.valuenum,
    le.valueuom
  FROM icustays ie
  JOIN labevents le
    ON le.subject_id = ie.subject_id
   AND le.itemid IN (50912, 52546)
   AND le.valuenum IS NOT NULL
   AND le.charttime BETWEEN ie.intime - INTERVAL 7 DAY AND ie.intime + INTERVAL 7 DAY
)
SELECT
  CAST(stay_id AS VARCHAR)                          AS id,
  date_diff('minute', intime, charttime) / 60.0     AS t_hours,
  valuenum                                          AS value
FROM src
-- 方案 4.2.2：0.1–50 mg/dL，剔除 0.0 与 808 这类录入错误
WHERE valuenum BETWEEN 0.1 AND 50
  AND valueuom = 'mg/dL';

-- ---- 钾 / 氯 —— 仅作阴性对照暴露，**不进任何主分析或主模型** ----
-- 用途：同一套暴露构造（0–24h 内末次减首次）与同一套基础模型，把暴露换成钾或
-- 氯再跑一遍。若与钠无关的电解质也呈现同等"保护"效应，说明观测到的并非钠本身，
-- 而是入院 24h 内化验整体回归正常 / 病情好转的伴随表现；若只有钠阳性，
-- 则支持效应有一定特异性。详见 advance.R 的负对照分析。
--
-- 条目：50971 生化钾、50902 生化氯（血/血浆，mEq/L）；52610 / 52535 是同义新条目，
-- 与血钠（52623）、肌酐（52546）的处理一致，低频但同源，合并入同表。
-- 范围：钾 1.0–10.0 mEq/L、氯 70–150 mEq/L（生理与可报告范围，剔除录入错误）。
CREATE OR REPLACE TABLE potassium AS
WITH src AS (
  SELECT
    ie.stay_id, ie.intime, le.charttime, le.valuenum, le.valueuom
  FROM icustays ie
  JOIN labevents le
    ON le.subject_id = ie.subject_id
   AND le.itemid IN (50971, 52610)
   AND le.valuenum IS NOT NULL
   AND le.charttime BETWEEN ie.intime - INTERVAL 7 DAY AND ie.intime + INTERVAL 7 DAY
)
SELECT
  CAST(stay_id AS VARCHAR)                          AS id,
  date_diff('minute', intime, charttime) / 60.0     AS t_hours,
  valuenum                                          AS value
FROM src
WHERE valuenum BETWEEN 1.0 AND 10.0
  AND valueuom = 'mEq/L';

CREATE OR REPLACE TABLE chloride AS
WITH src AS (
  SELECT
    ie.stay_id, ie.intime, le.charttime, le.valuenum, le.valueuom
  FROM icustays ie
  JOIN labevents le
    ON le.subject_id = ie.subject_id
   AND le.itemid IN (50902, 52535)
   AND le.valuenum IS NOT NULL
   AND le.charttime BETWEEN ie.intime - INTERVAL 7 DAY AND ie.intime + INTERVAL 7 DAY
)
SELECT
  CAST(stay_id AS VARCHAR)                          AS id,
  date_diff('minute', intime, charttime) / 60.0     AS t_hours,
  valuenum                                          AS value
FROM src
WHERE valuenum BETWEEN 70 AND 150
  AND valueuom = 'mEq/L';
