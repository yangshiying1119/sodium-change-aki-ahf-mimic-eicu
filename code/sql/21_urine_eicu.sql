-- =============================================================================
-- 21_urine_eicu.sql —— 尿量时序长表（S8：KDIGO 尿量标准）
--
-- 与 MIMIC 侧 21_urine_mimic.sql 同构，输出 (id, t_hours, value)：
--   id      = patientunitstayid（字符）
--   t_hours = intakeoutputoffset / 60，相对 ICU 入科小时数
--   value   = 该时点的尿量（mL），同一时点的多条记录先求和
-- 窗为 [0, 72h]——0–24h 用于"界标前是否已满足 KDIGO 尿量标准"，24–72h 用于结局。
--
-- ── 本库的尿量标识方式与 MIMIC 完全不同 ──────────────────────────────────────
-- MIMIC 用 outputevents.itemid 白名单；eICU 只有 flowsheet 的自由文本标签
-- （cellpath + celllabel），且同一含义有大量拼写/大小写变体（实测 'Urine'、
-- 'URINE CATHETER'、'Urine Output-Foley'、'Urine Output (mL)-Urethral Catheter'
-- 及带尾随空格的版本并存）。故不能靠精确等值列举。
--
-- 判据分两层，**两层缺一不可**：
--   (1) cellpath 必须落在 I&O 的出量区（'flowsheet|Flowsheet Cell Labels|I&O|
--       Output (ml)|' 前缀）——把入量区的条目挡在外面；
--   (2) 标签含 'urine'，且**不含**四类非计量条目：
--         count / occurrence —— 记录的是次数不是毫升（实测 13,250 + 5,197 行）
--         incontinen         —— 失禁，未计量
--         stool              —— 'Mixed Urine/Stool Volume'，混入粪便，非纯尿量
--       ⚠️ MIMIC 侧曾因"排除法"把胸管、胃管等非尿引流全算成尿量（0–24h 中位数
--          被抬到 2,480 mL）。此处虽用"含 urine"正向匹配，仍保留四类排除，
--          因为正向词本身也可能落在非计量条目上。
--
-- ⚠️ cellvaluenumeric 才是数值列；cellvaluetext 存的是文本描述。只取前者非空者。
-- =============================================================================

CREATE OR REPLACE TABLE urine_ts AS
SELECT
  CAST(t.patientunitstayid AS VARCHAR)                    AS id,
  io.intakeoutputoffset / 60.0                            AS t_hours,
  SUM(io.cellvaluenumeric)                                AS value
FROM cohort_C t
JOIN intakeOutput io
  ON  io.patientunitstayid = t.patientunitstayid
 AND  io.intakeoutputoffset >= 0
 AND  io.intakeoutputoffset <= 4320                       -- 72 小时
 AND  io.cellvaluenumeric IS NOT NULL
 AND  io.cellpath LIKE 'flowsheet|Flowsheet Cell Labels|I&O|Output (ml)|%'
 AND  LOWER(TRIM(io.celllabel)) LIKE '%urine%'
 AND  LOWER(TRIM(io.celllabel)) NOT LIKE '%count%'
 AND  LOWER(TRIM(io.celllabel)) NOT LIKE '%occurrence%'
 AND  LOWER(TRIM(io.celllabel)) NOT LIKE '%incontinen%'
 AND  LOWER(TRIM(io.celllabel)) NOT LIKE '%stool%'
GROUP BY t.patientunitstayid, io.intakeoutputoffset;
