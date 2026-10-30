-- =============================================================================
-- 41_chartevents_mimic.sql —— chartevents 瘦身提取（全库只扫这一遍）
--
-- chartevents 是 MIMIC-IV 里最大的表（约 3.3 亿行）。RRT 事件、生命体征、GCS、
-- 体重、透析问诊标记都要用它，若每个用途各查一次，同一份 CSV 会被反复解压扫描。
-- 故先在本文件用**一个带常量 itemid 谓词的扫描**产出 ce_slim，后续 42/43 都从
-- 它取数。谓词写成 WHERE 而非 JOIN 的 ON 条件，便于 DuckDB 下推到 CSV 读取层。
--
-- itemid 语义均已逐条比对 d_items 字典：
--   220050  Arterial Blood Pressure systolic   ┐ 收缩压。注意方案表 9 原文写的
--   220179  Non Invasive Blood Pressure systolic ┘ 220052/220181 实为**平均压**
--           （label 为 "... mean"），属方案笔误，此处按正确 itemid 取数并已回改方案。
--   220045  Heart Rate
--   220210  Respiratory Rate      224690 Respiratory Rate (Total)
--   220739/223900/223901  GCS Eye / Verbal / Motor
--   226512  Admission Weight (Kg)   224639 Daily Weight
--   227290  CRRT mode               225126 Dialysis patient
-- =============================================================================

CREATE OR REPLACE TABLE ce_slim AS
SELECT
  stay_id,
  charttime,
  itemid,
  valuenum,
  value
FROM chartevents
WHERE itemid IN (
  220050, 220179,             -- 收缩压（动脉 / 无创）
  220045,                     -- 心率
  220210, 224690,             -- 呼吸频率
  220739, 223900, 223901,     -- GCS 睁眼 / 语言 / 运动
  226512, 224639,             -- 体重（入科 / 每日）
  227290,                     -- CRRT 模式
  225126                      -- 透析患者（入科问诊）
)
AND charttime IS NOT NULL;
