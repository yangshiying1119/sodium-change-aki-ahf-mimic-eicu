-- =============================================================================
-- 30_cohort_eicu.sql —— eICU-CRD 2.0 侧界标队列与表型口径
--
-- 与 MIMIC 侧 30_cohort_mimic.sql 严格同构：输出 cohort_A / cohort_B / cohort_C
-- 三张表，列名一律对齐（**故意把 patientunitstayid 别名为 stay_id**），
-- 使 r/cohort.R 的判定层与 r/extract.R 的装配代码一行都不用改。
--
-- 本文件只负责 N1、N2、N4 的判定依据与排除标准①；N3/N5–N8 依赖检验序列，在 R 侧按
-- 同一套共享判定层（common.R）完成。
--
-- ── 与 MIMIC 侧的口径对应（表型 A/B/C）────────────────────────────────────────
-- MIMIC 用 diagnoses_icd.seq_num 表达诊断位次，**eICU 没有位次概念**，
-- diagnosis 表也不带 seq。故三个口径按下列方式对应，差异在论文中说明：
--   A_eicu = 入科问诊即为 CHF（admissionDx，相当于"主诊断"）
--            **且**存在急性心衰特异子码
--   B_eicu = 心衰并集（问诊 or 编码）**且**存在急性心衰特异子码   ← 敏感性口径（S6：210 例/44 事件）
--   C_eicu = 心衰并集（问诊 or 编码，不限急性/慢性）              ← 主分析口径（B 层在本库检验效能不足，
--                                                                  实际选用由 R 侧 _run9.R 的 EICU_CC="C" 控制）
-- 即：本库用"是否存在急性子码"取代 MIMIC 的"位次"来表达严格程度。
-- 因此 S6（改用严格表型）在两库间并非同一操作，论文中按此说明。
--
-- ── 时点体系 ────────────────────────────────────────────────────────────────
-- 本库所有时间字段都是"相对 ICU 入科分钟的偏移"，t0 = 0。故界标 = 1440 分钟、
-- 结局窗上界 = 4320 分钟，不存在 MIMIC 侧的绝对时间戳运算。
-- =============================================================================

-- ---- N1 基础层：成人（方案 3.2①）。停留时长与死亡留到 N4 判，不在此处筛 ----
-- age 是 VARCHAR：'> 89' 按 HIPAA 去标识规则赋 90；空值/异常为 NULL 后按 <18 剔除。
-- ⚠️ 不能直接 CAST——' > 89' 会静默变成 NULL 并被"<18"当成成人之外的另一种情况，
--    两件事混在一起，流程图上的 N0→N1 差额就说不清是年龄缺失还是未成年。
-- ⚠️ 本库没有绝对时间戳，而 R 侧有两处要求"真时间"：
--    · cohort.R 的 order(subject_id, intime) 选取"首次符合条件的住院"（纳入②）；
--    · competing.R 的 build_competing() 用 difftime(outtime, intime) 定竞争事件时刻。
--    故先造一个**基准时刻** admit_base，再在其上叠加单元内偏移得到三个真 TIMESTAMP。
--    基准取 "出院年份×一年 − hospitaladmitoffset"：年内精度来自 hospitaladmitoffset，
--    其绝对值不超过 525,600 分钟，两个年份的区间不会交叠，跨年意义上基准随真实入科
--    先后单调，(a) 与 (b) 同时满足。同年多次住院者的先后仅由 hospitaladmitoffset
--    （入院→入 ICU 间隔）决定，同年且同偏移的并列极端情形顺序不定——只影响"首次
--    住院"类判定的理论边角，界标/死亡/difftime 判定不受影响。年份缺失归 0（排最前），
--    不让 NA 静默落到
--    order() 的 na.last 位置上。三个时刻同基准，difftime 之差就是偏移量之差，
--    与 MIMIC 侧的绝对时间戳在判定上完全等价。
CREATE OR REPLACE TABLE base AS
WITH admit AS (
  SELECT
    p.*,
    TIMESTAMP '1970-01-01 00:00:00'
      + to_minutes(CAST(COALESCE(p.hospitaldischargeyear, 0) * 525600.0
                        - COALESCE(p.hospitaladmitoffset, 0) AS BIGINT)) AS admit_base
  FROM patient p
)
SELECT
  p.patientunitstayid,
  -- **两个别名是刻意加的**：R 侧判定层按 subject_id 排序取"首次符合条件的住院"、
  -- 按 hadm_id 关联住院级协变量，若此处改名，cohort.R 会静默落到 NULL 上。
  p.uniquepid                 AS subject_id,      -- 跨院患者标识
  p.patienthealthsystemstayid AS hadm_id,         -- 本次住院标识
  p.unitvisitnumber,                              -- 本次住院内的第几次 ICU 单元
  -- 医院标识（2026-10-31 新增，供 cluster_eicu.R 聚类稳健标准误用；
  -- v1.0.1 归档版无此列，复现聚类分析须连带 r/cluster_eicu.R 与 _runEG.R）
  p.hospitalid,
  -- 性别归一成 M/F：本库存 'Male'/'Female'，而 r/tables.R 的基线表按水平 'M'
  -- 统计"男性"一行（模型也照此建因子）。不归一的话 eICU 的基线表会把"男性"
  -- 印成 0（0.0%），既不报错也看不出来——两库并排成表时尤其危险。
  CASE WHEN p.gender = 'Male'   THEN 'M'
       WHEN p.gender = 'Female' THEN 'F'
       ELSE NULL END AS gender,
  p.ethnicity,
  p.unittype,
  p.unitadmitsource,
  -- 两个别名：r/tables.R 的 BASELINE_VARS 按 MIMIC 的列名取"ICU 类型"与"入院类型"，
  -- 此处对齐后，两库的基线表行结构完全一致，可直接并排成表 1。
  p.unittype       AS first_careunit,
  p.unitadmitsource AS admission_type,
  p.admissionweight,
  p.unitdischargeoffset,              -- 分钟；非严格单调，有负值
  p.unitdischargeoffset / 60.0 AS los_hours,
  p.unitdischargestatus,
  p.hospitaldischargeoffset,
  p.hospitaldischargestatus,
  p.hospitaldischargeyear,
  p.hospitaladmitoffset,
  CASE WHEN TRY_CAST(p.age AS INTEGER) IS NOT NULL THEN TRY_CAST(p.age AS INTEGER)
       WHEN TRIM(COALESCE(p.age, '')) LIKE '>%'  THEN 90
       ELSE NULL END AS age,
  p.admit_base                                       AS intime,
  p.admit_base + to_minutes(p.unitdischargeoffset)   AS outtime,
  -- 死亡时刻：本库只有"出院状态 + 偏移"，没有死亡时刻字段。凡院内死亡者，
  -- 其院内出院偏移即死亡时刻（覆盖"在 ICU 内死亡"与"转出后在病房死亡"两类）。
  -- 两个出院状态任一为 Expired 即算死亡；院内出院偏移缺失时（如仅记录单元内死亡）
  -- 回退到单元出院偏移，避免这一类被静默判成"存活"。
  CASE WHEN p.unitdischargestatus = 'Expired'
         OR p.hospitaldischargestatus = 'Expired'
       THEN p.admit_base
              + to_minutes(COALESCE(p.hospitaldischargeoffset, p.unitdischargeoffset))
       ELSE NULL END                                 AS deathtime,
  p.hospitaladmitoffset AS t_admit_offset
FROM admit p
WHERE p.unitdischargeoffset IS NOT NULL
  AND p.unitdischargeoffset > 0                  -- 剔除非正停留（脏数据）
  AND CASE WHEN TRY_CAST(p.age AS INTEGER) IS NOT NULL THEN TRY_CAST(p.age AS INTEGER)
           WHEN TRIM(COALESCE(p.age, '')) LIKE '>%'  THEN 90
           ELSE NULL END >= 18;

-- ---- N4 判定依据：界标（t0+24h）时仍在 ICU 且存活（方案 3.2④） ----
-- ⚠️ eICU **没有死亡时刻字段**，"死亡"只能由出院状态 + 偏移量还原（见 base 的
--    deathtime 一列）。该还原在**方向上是安全的**：它能识别出每一个院内死亡者，
--    但院内出院偏移对"在病房死亡"者是出院时刻而非抢救停止时刻，两者在本库无法区分，
--    故"24 h 内死亡"在极少数个案上可能略有偏早。该局限在论文中说明。
-- ⚠️ 不能用 hospitaldischargestatus 单独判定：它无法区分"在 ICU 内死亡"与
--    "转出后死亡"，而方案表 6 要求把这两种流失分开列。
-- 判定式与 MIMIC 侧逐字相同，只把绝对时间戳换成基准时刻上的同一套运算。
CREATE OR REPLACE TABLE landmark_flags AS
SELECT
  patientunitstayid,
  CASE WHEN outtime > intime + INTERVAL '24 hours' THEN 1 ELSE 0 END AS still_in_icu,
  CASE WHEN deathtime IS NOT NULL
         AND deathtime <= intime + INTERVAL '24 hours'
       THEN 1 ELSE 0 END AS died_within_24h
FROM base;

-- ---- 排除标准①：透析依赖或终末期肾病（方案 3.4①） ----
-- 两条路径，与 MIMIC 侧同构：
--   (a) 诊断码路径（下表 excl_esrd）—— 由 10_hf_phenotype_eicu.sql 建出。
--       ⚠️ 已知本库 Z99.2 / V45.1x 等透析依赖状态码**全库 0 行**，该路径近乎空转，
--       保留是为了留下"确实查过"的记录。
--   (b) 治疗记录路径（下表 esrd_chart）—— 本库 ESRD 的实际主源：
--       treatmentstring 中含 'for chronic renal failure' 的透析记录，
--       即把透析明确记为**慢性肾衰的对症治疗**者。
--       **不取** 'insertion of venous catheter for hemodialysis'（置管操作）与
--       'dialysis access surgery'（造瘘）——那是准备透析，不等于已在长期透析。
CREATE OR REPLACE TABLE esrd_chart AS
SELECT DISTINCT t.patientunitstayid AS stay_id
FROM treatment t
WHERE LOWER(t.treatmentstring) LIKE '%dialysis%'
  AND LOWER(t.treatmentstring) LIKE '%for chronic renal failure%';

-- ---- N2 心衰并集（方案 2.6 的两条识别路径取并集） ----
-- 路径①：入科问诊 admissionDx；路径②：诊断编码（ICD-9 与 ICD-10 双模式）。
-- 急性心衰特异子码清单与 MIMIC 侧同一份，按本库的小数点写法展开。
CREATE OR REPLACE TABLE hf_flags AS
SELECT
  p.patientunitstayid,
  MAX(CASE WHEN a.patientunitstayid IS NOT NULL THEN 1 ELSE 0 END) AS by_admitdx,
  MAX(COALESCE(d.hf_any, 0))                                       AS hf_any,
  MAX(COALESCE(d.acute_hf, 0))                                     AS acute_hf
FROM patient p
LEFT JOIN hf_admitdx  a ON a.patientunitstayid = p.patientunitstayid
LEFT JOIN hf_dx_flags d ON d.patientunitstayid = p.patientunitstayid
GROUP BY p.patientunitstayid;

-- ---- N2 三个表型口径的队列（带 N4 与排除①标志，供 R 侧逐层筛） ----
-- 三张表列结构完全一致，且与 MIMIC 侧同名，R 侧 .all_cohorts() 可原样复用。
-- **stay_id 别名是刻意加的**：R 侧 .all_cohorts() 与 run_mimic_cohort() 一律按
-- stay_id 取人，别名在此对齐即可，判定层与装配层一行都不用改。
CREATE OR REPLACE TABLE cohort_A AS          -- 严格口径，敏感性分析用
SELECT b.*, b.patientunitstayid AS stay_id,
       lf.still_in_icu,
       lf.died_within_24h,
       CASE WHEN e.patientunitstayid IS NULL THEN 0 ELSE 1 END AS esrd_excluded
FROM base b
JOIN hf_flags f ON f.patientunitstayid = b.patientunitstayid
JOIN landmark_flags lf ON lf.patientunitstayid = b.patientunitstayid
LEFT JOIN excl_esrd e ON e.patientunitstayid = b.patientunitstayid
WHERE f.by_admitdx = 1 AND f.acute_hf = 1;

CREATE OR REPLACE TABLE cohort_B AS          -- 敏感性口径（S6）；主分析用 cohort_C
SELECT b.*, b.patientunitstayid AS stay_id,
       lf.still_in_icu,
       lf.died_within_24h,
       CASE WHEN e.patientunitstayid IS NULL THEN 0 ELSE 1 END AS esrd_excluded
FROM base b
JOIN hf_flags f ON f.patientunitstayid = b.patientunitstayid
JOIN landmark_flags lf ON lf.patientunitstayid = b.patientunitstayid
LEFT JOIN excl_esrd e ON e.patientunitstayid = b.patientunitstayid
WHERE (f.by_admitdx = 1 OR f.hf_any = 1) AND f.acute_hf = 1;

CREATE OR REPLACE TABLE cohort_C AS          -- 主分析口径（R 侧 EICU_CC="C"）
SELECT b.*, b.patientunitstayid AS stay_id,
       lf.still_in_icu,
       lf.died_within_24h,
       CASE WHEN e.patientunitstayid IS NULL THEN 0 ELSE 1 END AS esrd_excluded
FROM base b
JOIN hf_flags f ON f.patientunitstayid = b.patientunitstayid
JOIN landmark_flags lf ON lf.patientunitstayid = b.patientunitstayid
LEFT JOIN excl_esrd e ON e.patientunitstayid = b.patientunitstayid
WHERE f.by_admitdx = 1 OR f.hf_any = 1;
