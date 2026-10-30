-- =============================================================================
-- 43_covariates_mimic.sql —— 方案 4.4 规定的协变量（stay_id 粒度）
--
-- 变量集合严格按方案 4.4 定死，**不做单变量筛选、不做逐步回归**：
--   基础模型：年龄、性别、入科首次血钠、界标前最低肌酐、慢性肾病、糖尿病、
--             早期最低收缩压、高血压、冠心病、心房颤动、心脏瓣膜病
--   扩展模型：基础 + 机械通气、血管活性药、利尿剂、液体出入量、
--             血糖最差值、血红蛋白最低值、心率/呼吸频率最差值、GCS 最低值
-- 其中年龄/性别来自 30_cohort 的 base；首次血钠、界标前最低肌酐、血钠检测次数
-- 由长表在 R 侧算；本文件负责其余部分。
--
-- 口径：一律取与暴露窗 [0, 24h]**有时间重叠**的事件。血管活性药与机械通气用
--       区间重叠（starttime < t0+24h AND endtime > t0），因为跨 t0 的持续输注/
--       通气在临床上属于"这 24 小时内在用"；其余按事件时刻落在窗内。
--
-- ⚠️ 每个聚合各自成表，最后才 LEFT JOIN 汇总。**绝不能**把 inputevents 与
--    outputevents 放进同一条 SELECT 的多个 LEFT JOIN 里——会按行做笛卡尔积，
--    SUM 被成倍放大。
--
-- 目标人群取 cohort_C：A ⊆ B ⊆ C，一份协变量表覆盖三个表型口径。
--
-- 有意不纳入（方案 4.4）：任何含肌酐或血钠的整体严重程度评分（SAPS/APACHE/
-- OASIS/SOFA）。MIMIC-IV 侧本就无现成评分表，此处也不自行构造。
-- 另注：Charlson 指数**不在方案 4.4 的变量集合内**，故不进主模型；其脚本
-- （49_charlson_mimic.sql）保留备查，仅作额外的合并症负担敏感性分析。
-- =============================================================================

-- ---- 目标人群 = cohort_C（含 subject_id 与窗口上界） ----
CREATE OR REPLACE TABLE cov_target AS
SELECT
  stay_id,
  hadm_id,
  subject_id,
  intime,
  intime + INTERVAL '24 hours' AS t24
FROM cohort_C;

-- ---- 合并症（diagnoses_icd，hadm 粒度） ----
-- 编码规则按方案表 9。MIMIC-IV 的 ICD 码无小数点，故 ICD-10 用前 3 位、ICD-9 用
-- 前 3 位或全码匹配。
CREATE OR REPLACE TABLE cov_icd AS
SELECT
  hadm_id,
  MAX(CASE WHEN SUBSTR(icd_code, 1, 3) = 'N18'
            OR SUBSTR(icd_code, 1, 3) = '585'      THEN 1 ELSE 0 END) AS ckd,
  MAX(CASE WHEN SUBSTR(icd_code, 1, 3) IN ('E10','E11','E12','E13','E14')
            OR SUBSTR(icd_code, 1, 3) = '250'      THEN 1 ELSE 0 END) AS diabetes,
  MAX(CASE WHEN SUBSTR(icd_code, 1, 3) IN ('I10','I11','I12','I13','I15')
            OR SUBSTR(icd_code, 1, 3) IN ('401','402','403','404','405')
                                                   THEN 1 ELSE 0 END) AS hypertension,
  MAX(CASE WHEN SUBSTR(icd_code, 1, 3) IN ('I20','I21','I22','I23','I24','I25')
            OR SUBSTR(icd_code, 1, 3) IN ('410','411','412','413','414')
                                                   THEN 1 ELSE 0 END) AS chd,
  MAX(CASE WHEN SUBSTR(icd_code, 1, 3) = 'I48'
            OR icd_code = '42731'                  THEN 1 ELSE 0 END) AS atrial_fib,
  MAX(CASE WHEN SUBSTR(icd_code, 1, 3) IN ('I05','I06','I07','I08',
                                           'I34','I35','I36','I37','I39')
            OR SUBSTR(icd_code, 1, 3) IN ('394','395','396','397','424')
                                                   THEN 1 ELSE 0 END) AS valvular
FROM diagnoses_icd
GROUP BY hadm_id;

-- ---- 生命体征：0–24h 最低收缩压 / 心率与呼吸频率最差值 ----
-- 收缩压取 220050（动脉）∪ 220179（无创）。范围限幅剔除明显伪值。
CREATE OR REPLACE TABLE cov_vitals AS
SELECT
  t.stay_id,
  MIN(CASE WHEN ce.itemid IN (220050, 220179)
            AND ce.valuenum BETWEEN 20 AND 300 THEN ce.valuenum END) AS sbp_min_24h,
  MIN(CASE WHEN ce.itemid = 220045
            AND ce.valuenum BETWEEN 20 AND 300 THEN ce.valuenum END) AS hr_min_24h,
  MAX(CASE WHEN ce.itemid = 220045
            AND ce.valuenum BETWEEN 20 AND 300 THEN ce.valuenum END) AS hr_max_24h,
  MAX(CASE WHEN ce.itemid IN (220210, 224690)
            AND ce.valuenum BETWEEN 1 AND 100 THEN ce.valuenum END) AS rr_max_24h
FROM cov_target t
LEFT JOIN ce_slim ce
  ON  ce.stay_id = t.stay_id
  AND ce.charttime >= t.intime
  AND ce.charttime <= t.t24
GROUP BY t.stay_id;

-- ---- GCS 最低值 ----
-- 睁眼/语言/运动三项在同一时刻齐备时才求和，避免用缺项之和充当总分。
-- 注意：气管插管患者的语言项常记为 0 或 1，该局限在论文中说明。
CREATE OR REPLACE TABLE cov_gcs AS
SELECT stay_id, MIN(gcs_total) AS gcs_min_24h
FROM (
  SELECT
    ce.stay_id,
    ce.charttime,
    SUM(ce.valuenum)                        AS gcs_total,
    COUNT(DISTINCT ce.itemid)               AS n_comp
  FROM cov_target t
  JOIN ce_slim ce
    ON  ce.stay_id = t.stay_id
    AND ce.charttime >= t.intime
    AND ce.charttime <= t.t24
   AND ce.itemid IN (220739, 223900, 223901)
   AND ce.valuenum IS NOT NULL
  GROUP BY ce.stay_id, ce.charttime
)
WHERE n_comp = 3
GROUP BY stay_id;

-- ---- 体重（供液体出入量按 mL/(kg·h) 标化） ----
CREATE OR REPLACE TABLE cov_weight AS
SELECT
  t.stay_id,
  COALESCE(
    MIN(CASE WHEN ce.itemid = 226512 AND ce.valuenum BETWEEN 20 AND 400
             THEN ce.valuenum END),
    MIN(CASE WHEN ce.itemid = 224639 AND ce.valuenum BETWEEN 20 AND 400
             THEN ce.valuenum END)
  ) AS weight_kg
FROM cov_target t
LEFT JOIN ce_slim ce
  ON  ce.stay_id = t.stay_id
  AND ce.charttime <= t.t24
GROUP BY t.stay_id;

-- ---- 血管活性药 + 高渗盐水：是否使用 + 最大泵速 ----
-- 泵速只汇总 itemid 属于血管活性药清单、且 rateuom = 'mcg/kg/min' 的记录。
-- **两个条件缺一不可**：只按 rateuom 取 MAX 会把同样以 mcg/kg/min 计量的异丙酚、
-- 胰岛素、右美托咪定等非血管活性药的最大泵速混进来，该变量就不再是"血管活性药
-- 强度"了。
CREATE OR REPLACE TABLE cov_vaso AS
SELECT
  t.stay_id,
  MAX(CASE WHEN ie.itemid IN (221906, 221289, 221662, 221653,
                              221749, 222315, 221986) THEN 1 ELSE 0 END) AS vasoactive_any,
  MAX(CASE WHEN ie.itemid IN (221906, 221289, 221662, 221749, 222315)
           THEN 1 ELSE 0 END) AS vasopressor_any,
  MAX(CASE WHEN ie.itemid IN (221653, 221986) THEN 1 ELSE 0 END) AS inotrope_any,
  MAX(CASE WHEN ie.itemid IN (225161, 228341) THEN 1 ELSE 0 END) AS hypertonic_saline
FROM cov_target t
LEFT JOIN inputevents ie
  ON  ie.stay_id = t.stay_id
  AND ie.starttime < t.t24
  AND ie.endtime   > t.intime
GROUP BY t.stay_id;

-- ---- 血管活性药剂量：去甲肾上腺素当量最大速率（方案表 9 的"最大速率"连续项）----
-- **不能把 7 种药的 rate 放在同一列取最大**：多巴酚丁胺常用量 2.5–20 mcg/kg/min，
-- 去甲肾上腺素 0.05–0.5，米力农 0.375–0.75，量纲差一两个数量级；同列取最大等于
-- 把强心药当加压药计（实测该做法下队列内"最大速率"最高值 32.5 来自多巴酚丁胺）。
-- 另有两种量纲问题：血管加压素记在 units/hour（全库 37,160 行），
-- 去氧肾上腺素有按 mcg/min 录入却标 mcg/kg/min 的记录（最大 379）。
--
-- 故统一换算为去甲肾上腺素当量（NEE，mcg/kg/min）——ICU 文献通用换算
-- （Gershengorn 等；NE + 肾上腺素 + 多巴胺/100 + 去氧肾上腺素/10 + 血管加压素×2.5）。
-- 换算前按**各药公认临床区间**剔除单位/录入错误：固定阈值，不依赖本数据的分布
-- （方案 4.1.3）；越界记录只从剂量中剔除，不影响上面 vasoactive_any 的"是否使用"。
--
-- 取最大值的口径：NEE 是**同时段多种药之和**，只在"各药最大速率之和"上取最大会
-- 高估（把从未同时使用的药加在一起）。故先在每个时点求和，再取时间轴上的最大值。
-- 只在各输注的 starttime 上求值即可取到最大值：任一取值区间的左端点若非 start
-- 就是 end，而 end 处之值必已在某个更早的 start 处取得。
CREATE OR REPLACE TABLE cov_vaso_nee AS
WITH v AS (
  SELECT ie.stay_id, ie.starttime, ie.endtime,
    CASE
      WHEN ie.itemid = 221906 AND ie.rateuom = 'mcg/kg/min'
           AND ie.rate > 0 AND ie.rate <= 3        THEN ie.rate
      WHEN ie.itemid = 221289 AND ie.rateuom = 'mcg/kg/min'
           AND ie.rate > 0 AND ie.rate <= 3        THEN ie.rate
      WHEN ie.itemid = 221662 AND ie.rateuom = 'mcg/kg/min'
           AND ie.rate >= 0.1 AND ie.rate <= 50    THEN ie.rate / 100.0
      WHEN ie.itemid = 221749 AND ie.rateuom = 'mcg/kg/min'
           AND ie.rate >= 0.01 AND ie.rate <= 10   THEN ie.rate / 10.0
      WHEN ie.itemid = 222315 AND ie.rateuom = 'units/hour'
           AND ie.rate >= 0.1 AND ie.rate <= 10    THEN ie.rate / 60.0 * 2.5
      WHEN ie.itemid = 222315 AND ie.rateuom = 'units/min'
           AND ie.rate >= 0.002 AND ie.rate <= 0.17 THEN ie.rate * 2.5
    END AS nee
  FROM cov_target t
  JOIN inputevents ie
    ON  ie.stay_id = t.stay_id
   AND ie.starttime < t.t24
   AND ie.endtime   > t.intime
  WHERE ie.itemid IN (221906, 221289, 221662, 221749, 222315)
),
pts AS (
  SELECT DISTINCT stay_id, starttime AS t FROM v WHERE nee IS NOT NULL
),
summed AS (
  SELECT p.stay_id, p.t, SUM(v.nee) AS nee_t
  FROM pts p
  JOIN v ON v.stay_id = p.stay_id
        AND v.starttime <= p.t
        AND (v.endtime IS NULL OR v.endtime > p.t)
  WHERE v.nee IS NOT NULL
  GROUP BY p.stay_id, p.t
)
SELECT stay_id, MAX(nee_t) AS nee_max FROM summed GROUP BY stay_id;

-- ---- 利尿剂（静脉部分）：inputevents 的实际给药量 ----
-- 221794 / 228340 是两种规格的呋塞米，229639 是布美他尼。只取 amountuom = 'mg'。
CREATE OR REPLACE TABLE cov_diur_iv AS
SELECT
  t.stay_id,
  MAX(CASE WHEN ie.itemid IN (221794, 228340, 229639) THEN 1 ELSE 0 END) AS diur_iv_any,
  COALESCE(SUM(CASE WHEN ie.itemid IN (221794, 228340) AND ie.amountuom = 'mg'
                    THEN ie.amount
                    WHEN ie.itemid = 229639 AND ie.amountuom = 'mg'
                    THEN ie.amount * 40
                    ELSE 0 END), 0) AS diur_iv_eq_mg
FROM cov_target t
LEFT JOIN inputevents ie
  ON  ie.stay_id = t.stay_id
  AND ie.starttime < t.t24
  AND ie.endtime   > t.intime
GROUP BY t.stay_id;

-- ---- 利尿剂（口服/肠内部分）：prescriptions 药名匹配（方案表 9 明文要求） ----
-- inputevents 只覆盖静脉给药，口服袢利尿剂必须从 prescriptions 取，否则出院前
-- 已转为口服的患者会被记为"未用利尿剂"。
--
-- **与静脉部分严格互斥**：此处只取非静脉途径（route 不在 IV / IV DRIP 内）。
-- 一次静脉给药会同时留下 prescriptions 医嘱与 inputevents 执行记录，两边都算会
-- 把这部分剂量的权重翻倍。
--
-- 折算为**呋塞米当量毫克**（口服）：呋塞米 ×1、布美他尼 ×40、托拉塞米 ×2、
-- 依他尼酸 ×0.8（40/50）。24 小时剂量 = 单次剂量 × doses_per_24_hrs。
-- 局限：静脉与口服的等效比并不相同（静脉约强一倍），此处按口服当量统一折算，
-- 会把静脉部分的暴露强度低估约一半；该变量只进扩展模型作稳健性分析、不作因果
-- 解释（方案 4.4），此局限在论文中说明。
-- 药名用词干匹配（furosem / bumetan / torsem / ethacryn），以覆盖库中
-- "Furosemi"、"Bumetanid" 一类截断写法与复方制剂名。
CREATE OR REPLACE TABLE cov_diur_po AS
SELECT
  t.stay_id,
  CASE WHEN COUNT(*) > 0 THEN 1 ELSE 0 END AS diur_po_any,
  SUM(p.dose_mg * p.eq_factor * COALESCE(p.doses_per_24_hrs, 1)) AS diur_po_eq_mg
FROM cov_target t
JOIN (
  SELECT
    hadm_id,
    starttime,
    stoptime,
    TRY_CAST(regexp_extract(dose_val_rx, '[0-9]+(\.[0-9]+)?', 0) AS DOUBLE) AS dose_mg,
    CASE WHEN regexp_matches(lower(drug), 'furosem')  THEN 1.0
         WHEN regexp_matches(lower(drug), 'bumetan')  THEN 40.0
         WHEN regexp_matches(lower(drug), 'torsem')   THEN 2.0
         WHEN regexp_matches(lower(drug), 'ethacryn') THEN 0.8
         ELSE 0.0 END AS eq_factor,
    doses_per_24_hrs
  FROM prescriptions
  WHERE regexp_matches(lower(drug), 'furosem|bumetan|torsem|ethacryn')
    AND dose_unit_rx = 'mg'
    -- 只取肠内途径（白名单，不用"非静脉"黑名单：后者会放进 IH 吸入等非肠内途径）
    AND upper(COALESCE(route, '')) IN ('PO', 'PO/NG', 'NG', 'ORAL')
) p
  ON  p.hadm_id   = t.hadm_id
  AND p.starttime < t.t24
  AND (p.stoptime IS NULL OR p.stoptime > t.intime)
GROUP BY t.stay_id;

-- ---- 利尿剂汇总 ----
CREATE OR REPLACE TABLE cov_meds AS
SELECT
  v.stay_id,
  v.vasoactive_any, v.vasopressor_any, v.inotrope_any, v.hypertonic_saline,
  COALESCE(n.nee_max, 0) AS nee_max,
  CASE WHEN COALESCE(iv.diur_iv_any, 0) = 1
         OR COALESCE(po.diur_po_any, 0) = 1 THEN 1 ELSE 0 END AS diuretic_any,
  COALESCE(iv.diur_iv_eq_mg, 0) AS diuretic_iv_eq_mg,
  COALESCE(po.diur_po_eq_mg, 0) AS diuretic_po_eq_mg,
  COALESCE(iv.diur_iv_eq_mg, 0) + COALESCE(po.diur_po_eq_mg, 0) AS diuretic_furosemide_eq_mg
FROM cov_vaso v
LEFT JOIN cov_vaso_nee n USING (stay_id)
LEFT JOIN cov_diur_iv iv USING (stay_id)
LEFT JOIN cov_diur_po po USING (stay_id);

-- ---- 机械通气（procedureevents） ----
CREATE OR REPLACE TABLE cov_vent AS
SELECT
  t.stay_id,
  MAX(CASE WHEN pe.itemid = 225792 THEN 1 ELSE 0 END) AS mech_vent_invasive,
  MAX(CASE WHEN pe.itemid = 225794 THEN 1 ELSE 0 END) AS mech_vent_noninvasive,
  MAX(CASE WHEN pe.itemid IN (225792, 225794) THEN 1 ELSE 0 END) AS mech_vent_any
FROM cov_target t
LEFT JOIN procedureevents pe
  ON  pe.stay_id = t.stay_id
  AND pe.starttime < t.t24
  AND pe.endtime   > t.intime
GROUP BY t.stay_id;

-- ---- 0–24h 液体入量（mL） ----
CREATE OR REPLACE TABLE cov_fluid_in AS
SELECT
  t.stay_id,
  COALESCE(SUM(CASE WHEN ie.amountuom = 'mL' THEN ie.amount ELSE 0 END), 0) AS fluid_in_24h_ml
FROM cov_target t
LEFT JOIN inputevents ie
  ON  ie.stay_id = t.stay_id
  AND ie.starttime < t.t24
  AND ie.endtime   > t.intime
GROUP BY t.stay_id;

-- ---- 0–24h 尿量（mL） ----
-- 尿路 itemid 清单来自方案 4.4，**必须显式列举**（清单与 21_urine_mimic.sql 一致）。
-- 其中 227489（GU Irrigant/Urine Volume Out）是膀胱冲洗时的**混合出量**，须扣除
-- 同一时段的冲洗液入量 227488 才是净尿量。方案原文写的 226566 在 MIMIC-IV v3.1
-- 中没有任何记录，已订正为 227489。
--
-- ⚠️ 此处早先误写成 `oe.itemid <> 227488`（"除冲洗液外的全部出量"）。outputevents
--    里胸管、纵隔/胸腔/腹腔引流、胃管、粪便、伤口负压等条目占了很大比重，
--    该写法把它们全算成了尿量，0–24h "尿量"中位数被抬到 2,480 mL。必须靠白名单。
CREATE OR REPLACE TABLE cov_urine_out AS
SELECT
  t.stay_id,
  COALESCE(SUM(CASE WHEN oe.itemid = 227488 THEN -oe.value ELSE oe.value END), 0)
    AS urine_out_24h_ml
FROM cov_target t
LEFT JOIN outputevents oe
  ON  oe.stay_id = t.stay_id
  AND oe.itemid IN (
        226559, 226560, 226561, 226563, 226564, 226565, 226567, 226584,
        226557, 226558, 226627, 226631, 227489, 227488)
  AND oe.charttime >= t.intime
  AND oe.charttime <= t.t24
GROUP BY t.stay_id;

-- ---- 0–24h 血糖与血红蛋白 ----
-- 血糖 50931（生化）∪ 50809（血气）；血红蛋白 51222（血液学）∪ 50811（血气）。
-- 方案 4.1.2 特别指出 50809 存在量级异常记录（最大 1276103 mg/dL），故限幅。
CREATE OR REPLACE TABLE cov_labs AS
SELECT
  t.stay_id,
  MAX(CASE WHEN le.itemid IN (50931, 50809)
            AND le.valuenum BETWEEN 20 AND 1000 THEN le.valuenum END) AS glucose_max_24h,
  MIN(CASE WHEN le.itemid IN (50931, 50809)
            AND le.valuenum BETWEEN 20 AND 1000 THEN le.valuenum END) AS glucose_min_24h,
  MIN(CASE WHEN le.itemid IN (51222, 50811)
            AND le.valuenum BETWEEN 3 AND 25 THEN le.valuenum END) AS hgb_min_24h
FROM cov_target t
LEFT JOIN labevents le
  ON  le.subject_id = t.subject_id
  AND le.itemid IN (50931, 50809, 51222, 50811)
  AND le.charttime >= t.intime
  AND le.charttime <= t.t24
  AND le.valuenum IS NOT NULL
GROUP BY t.stay_id;

-- ---- 汇总 ----
CREATE OR REPLACE TABLE cov_all AS
SELECT
  b.stay_id,
  b.hadm_id,
  b.subject_id,
  -- 合并症
  i.ckd, i.diabetes, i.hypertension, i.chd, i.atrial_fib, i.valvular,
  -- 生命体征与神经
  v.sbp_min_24h, v.hr_min_24h, v.hr_max_24h, v.rr_max_24h,
  g.gcs_min_24h,
  w.weight_kg,
  -- 治疗强度
  m.vasoactive_any, m.vasopressor_any, m.inotrope_any, m.nee_max,
  m.hypertonic_saline,
  -- 利尿剂：总量 + 静脉/口服两个分量（分量留作核查，主模型只用总量）
  m.diuretic_any, m.diuretic_furosemide_eq_mg,
  m.diuretic_iv_eq_mg, m.diuretic_po_eq_mg,
  ve.mech_vent_invasive, ve.mech_vent_noninvasive, ve.mech_vent_any,
  -- 出入量
  fi.fluid_in_24h_ml,
  uo.urine_out_24h_ml,
  fi.fluid_in_24h_ml - uo.urine_out_24h_ml AS fluid_balance_24h_ml,
  -- 化验
  lb.glucose_max_24h, lb.glucose_min_24h, lb.hgb_min_24h
FROM cov_target b
LEFT JOIN cov_icd       i  ON i.hadm_id = b.hadm_id
LEFT JOIN cov_vitals    v  USING (stay_id)
LEFT JOIN cov_gcs       g  USING (stay_id)
LEFT JOIN cov_weight    w  USING (stay_id)
LEFT JOIN cov_meds      m  USING (stay_id)
LEFT JOIN cov_vent      ve USING (stay_id)
LEFT JOIN cov_fluid_in  fi USING (stay_id)
LEFT JOIN cov_urine_out uo USING (stay_id)
LEFT JOIN cov_labs      lb USING (stay_id);
