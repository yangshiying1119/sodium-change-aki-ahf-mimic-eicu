-- =============================================================================
-- 42_rrt_mimic.sql —— 肾脏替代治疗事件 + 透析问诊标记
--
-- rrt        长表 (id, t_hours, value = 1)，t_hours 相对 ICU 入科时刻。用途：
--              (a) 排除标准③：界标前（t ≤ 24h）已开始 RRT 者不进分析
--              (b) 次要结局：界标后（24 < t ≤ 72h）启动 RRT 者按 KDIGO 3 期计
-- esrd_chart 排除标准①的第三条路径：入科问诊记录为"透析患者"。
--              单列成表是为了能单独报告它的边际排除人数——它与诊断码路径可能重叠，
--              若过度排除可单独回退，不影响另两条路径。
--
-- ⚠️ 单靠 procedureevents 会**低估**。MIMIC-IV v3.1 中 CVVH / CVVHD / CVVHDF /
--    SCUF 各模式没有独立 itemid，全部记在 225802 "Dialysis - CRRT" 名下
--    （方案原文所列 225803 / 225809 / 225955 在库中 0 行）。其起止时刻另由
--    chartevents 227290 "CRRT mode" 补足。故本表两源并用。
-- =============================================================================

-- ---- RRT 事件时刻 ----
-- 以过程**开始时刻**为准：跨 t0 的持续治疗按"这段 24 小时内已在做"计。
CREATE OR REPLACE TABLE rrt AS
WITH proc AS (
  SELECT
    ie.stay_id,
    ie.intime,
    pe.starttime AS t
  FROM icustays ie
  JOIN procedureevents pe
    ON pe.stay_id = ie.stay_id
   AND pe.itemid IN (225802,   -- Dialysis - CRRT（含 CVVH/CVVHD/CVVHDF/SCUF 各模式）
                     225441,   -- Hemodialysis
                     225805)   -- Peritoneal Dialysis
   AND pe.starttime BETWEEN ie.intime - INTERVAL 7 DAY AND ie.intime + INTERVAL 7 DAY
),
chart AS (
  SELECT
    ie.stay_id,
    ie.intime,
    ce.charttime AS t
  FROM icustays ie
  JOIN ce_slim ce
    ON ce.stay_id = ie.stay_id
   AND ce.itemid = 227290          -- CRRT mode
   AND ce.charttime BETWEEN ie.intime - INTERVAL 7 DAY AND ie.intime + INTERVAL 7 DAY
)
SELECT DISTINCT
  CAST(stay_id AS VARCHAR)                          AS id,
  date_diff('minute', intime, t) / 60.0             AS t_hours,
  1                                                 AS value
FROM (SELECT * FROM proc UNION ALL SELECT * FROM chart);

-- ---- 排除标准① 路径(c)：入科问诊"透析患者" ----
-- 225126 是勾选项，value 记为 'Yes'/'No'（少数为 1/0），两种写法都接。
CREATE OR REPLACE TABLE esrd_chart AS
SELECT DISTINCT stay_id
FROM ce_slim
WHERE itemid = 225126
  AND (LOWER(TRIM(CAST(value AS VARCHAR))) IN ('yes', '1')
       OR valuenum = 1);
