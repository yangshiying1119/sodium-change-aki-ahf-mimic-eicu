-- =============================================================================
-- 21_urine_mimic.sql —— 尿量时序长表（S8：KDIGO 尿量标准）
--
-- 与 20_labs 同一套长表约定：(id, t_hours, value)，t_hours 相对 ICU 入科时刻。
-- 窗为 [intime, intime + 72h]——0–24h 用于"界标前是否已满足 KDIGO 尿量标准"的
-- 判定，24–72h 用于结局；两者必须同源同口径，否则敏感性分析不成立。
--
-- ⚠️ 尿路 itemid 必须显式列举。**曾漏过一次**：43_covariates 里 cov_urine_out
--    当时写成 `oe.itemid <> 227488` 即"除冲洗液入量外的全部出量"，把胸管、
--    腹腔引流、胃管、粪便、伤口负压全算进了"尿量"（0–24h 中位数被抬到 2480 mL）。
--    出量表里非尿液条目占比很高，绝不能靠排除法。
--
-- 227488（GU Irrigant Volume In，膀胱冲洗液入量）以**负值**进入同一条目，
-- 与 227489（GU Irrigant/Urine Volume Out，混合出量）在同一时刻相加即得净尿量；
-- 若两者时刻不同（冲洗与回抽分次记录），按 charttime 分组求和会得到负的中间值，
-- 由 R 侧在滚动窗内累计时自然抵消，不在此处强行配对。
-- =============================================================================

CREATE OR REPLACE TABLE urine_ts AS
SELECT
  CAST(t.stay_id AS VARCHAR)                        AS id,
  date_diff('minute', t.intime, oe.charttime) / 60.0 AS t_hours,
  SUM(CASE WHEN oe.itemid = 227488 THEN -oe.value ELSE oe.value END) AS value
FROM cohort_C t
JOIN outputevents oe
  ON  oe.stay_id = t.stay_id
 AND  oe.itemid IN (
        226559,  -- Foley
        226560,  -- Void
        226561,  -- Condom Cath
        226563,  -- Suprapubic
        226564,  -- R Nephrostomy
        226565,  -- L Nephrostomy
        226567,  -- Straight Cath
        226584,  -- Ileoconduit
        226557,  -- R Ureteral Stent
        226558,  -- L Ureteral Stent
        226627,  -- OR Urine
        226631,  -- PACU Urine
        227489,  -- GU Irrigant/Urine Volume Out
        227488   -- GU Irrigant Volume In（负值）
      )
 AND oe.charttime >= t.intime
 AND oe.charttime <= t.intime + INTERVAL '72 hours'
 AND oe.value IS NOT NULL
GROUP BY t.stay_id, date_diff('minute', t.intime, oe.charttime) / 60.0;
