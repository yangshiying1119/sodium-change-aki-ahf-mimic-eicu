-- =============================================================================
-- 40_covariates_eicu.sql —— 方案 4.4 规定的协变量（eICU 侧，与 MIMIC 同变量集）
--
-- 变量集合与 43_covariates_mimic.sql **逐一对应**，不做单变量筛选、不做逐步回归：
--   基础模型：年龄、性别、入科首次血钠、界标前最低肌酐、慢性肾病、糖尿病、
--             早期最低收缩压、高血压、冠心病、心房颤动、心脏瓣膜病
--   扩展模型：基础 + 机械通气、血管活性药、利尿剂、液体出入量、
--             血糖最差值、血红蛋白最低值、心率/呼吸频率最差值、GCS 最低值
-- 其中年龄/性别来自 30_cohort 的 base；首次血钠、界标前最低肌酐、血钠检测次数
-- 由长表在 R 侧算；本文件负责其余部分。
--
-- 口径：一律取与暴露窗 [0, 24h]（即 offset 0–1440 分钟）**有时间重叠**的事件。
--
-- ⚠️ 每个聚合各自成表，最后才 LEFT JOIN 汇总。**绝不能**把 intakeOutput 的入量与
--    出量放进同一条 SELECT 的两个 LEFT JOIN 里——会按行做笛卡尔积，SUM 成倍放大。
--
-- ⛔ 有意不纳入（方案 4.4）：任何含肌酐或血钠的整体严重程度评分。本库的
--    apacheApsVar / apachePatientResult 里有 apachescore 与 SOFA 类字段，但同一张
--    apacheApsVar 还直接含 sodium / creatinine / urine / dialysis —— 就是本研究的
--    暴露与结局本身。**整表不 join**，只从 00_schema_eicu.sql 里那张显式只投影
--    GCS 三列的 apache_gcs 视图取数。
-- =============================================================================

-- ---- 目标人群 = cohort_C（与 MIMIC 侧同：A ⊆ B ⊆ C，一份表覆盖三个口径） ----
CREATE OR REPLACE TABLE cov_target AS
SELECT
  patientunitstayid,
  patientunitstayid AS stay_id,
  hadm_id,
  subject_id,
  1440 AS t24                       -- 界标，分钟
FROM cohort_C;

-- ---- 合并症（诊断码，**按住院粒度**聚合） ----
-- 与 MIMIC 侧一样按**住院**而非 ICU 单元聚合：一个人在本次住院的其它单元里被编码的
-- 合并症，也应算作他的合并症。eICU 的 diagnosis 挂在 patientunitstayid 上，故先经
-- patient 表把它映射到 patienthealthsystemstayid（即 base 里的 hadm_id）。
--
-- ⚠️ 取数走 10_hf_phenotype_eicu.sql 建出的 **dx_code（一行一码）**，不直接在
--    icd9code 字段上匹配。原因见该文件的说明：icd9code 是一格多码的字符串，
--    在整串上写前缀通配符只能匹配到排在**第一位**的码，逗号之后的同码会被静默漏掉
--    ——不报错、不产生缺失值，只是把一批真患者算成"无该合并症"。拆行之后
--    `code LIKE '585%'` 才真正等价于"存在以 585 开头的码"。
-- ⚠️ 每个病种同时写 ICD-9 与 ICD-10 两版（本库实测 ICD-9 带小数点，如 '585.9'）。
CREATE OR REPLACE TABLE cov_icd AS
SELECT
  p.patienthealthsystemstayid AS hadm_id,
  MAX(CASE WHEN d.code LIKE '585%' OR d.code LIKE 'N18%'
                                                                 THEN 1 ELSE 0 END) AS ckd,
  MAX(CASE WHEN d.code LIKE '250%'
             OR d.code LIKE 'E1%'                                THEN 1 ELSE 0 END) AS diabetes,
  MAX(CASE WHEN d.code LIKE '401%'  OR d.code LIKE '402%'
             OR d.code LIKE '403%'  OR d.code LIKE '404%'
             OR d.code LIKE '405%'
             OR d.code LIKE 'I10%'  OR d.code LIKE 'I11%'
             OR d.code LIKE 'I12%'  OR d.code LIKE 'I13%'
             OR d.code LIKE 'I15%'                               THEN 1 ELSE 0 END) AS hypertension,
  MAX(CASE WHEN d.code LIKE '410%'  OR d.code LIKE '411%'
             OR d.code LIKE '412%'  OR d.code LIKE '413%'
             OR d.code LIKE '414%'
             OR d.code LIKE 'I20%'  OR d.code LIKE 'I21%'
             OR d.code LIKE 'I22%'  OR d.code LIKE 'I24%'
             OR d.code LIKE 'I25%'                               THEN 1 ELSE 0 END) AS chd,
  -- 心房颤动：ICD-9 427.31（本组里唯一带小数点的码），ICD-10 I48.x
  MAX(CASE WHEN d.code LIKE '427.31%' OR d.code LIKE '42731%'
             OR d.code LIKE 'I48%'                               THEN 1 ELSE 0 END) AS atrial_fib,
  MAX(CASE WHEN d.code LIKE '394%'  OR d.code LIKE '395%'
             OR d.code LIKE '396%'  OR d.code LIKE '397%'
             OR d.code LIKE '424%'
             OR d.code LIKE 'I05%'  OR d.code LIKE 'I06%'
             OR d.code LIKE 'I07%'  OR d.code LIKE 'I08%'
             OR d.code LIKE 'I34%'  OR d.code LIKE 'I35%'
             OR d.code LIKE 'I36%'  OR d.code LIKE 'I37%'
                                                                 THEN 1 ELSE 0 END) AS valvular
FROM dx_code d
JOIN patient p ON p.patientunitstayid = d.patientunitstayid
GROUP BY p.patienthealthsystemstayid;

-- ---- 生命体征：0–24h 最低收缩压 / 心率与呼吸频率最差值 ----
-- 收缩压取 vitalAperiodic 的无创压（覆盖 88%）∪ vitalPeriodic 的systemic 收缩压
-- （仅 19% 覆盖，作补充）。方案 4.1.3：生命体征按公认临床区间限幅，剔除伪值。
CREATE OR REPLACE TABLE cov_vitals AS
SELECT
  t.stay_id,
  MIN(CASE WHEN v.sbp BETWEEN 20 AND 300 THEN v.sbp END) AS sbp_min_24h
FROM cov_target t
LEFT JOIN (
  SELECT patientunitstayid, noninvasivesystolic AS sbp, observationoffset AS off
  FROM vitalAperiodic WHERE noninvasivesystolic IS NOT NULL
  UNION ALL
  SELECT patientunitstayid, systemicsystolic AS sbp, observationoffset AS off
  FROM vitalPeriodic WHERE systemicsystolic IS NOT NULL
) v ON v.patientunitstayid = t.patientunitstayid
   AND v.off >= 0 AND v.off <= t.t24
GROUP BY t.stay_id;

CREATE OR REPLACE TABLE cov_vitals2 AS
SELECT
  t.stay_id,
  MIN(CASE WHEN vp.heartrate   BETWEEN 20 AND 300 THEN vp.heartrate   END) AS hr_min_24h,
  MAX(CASE WHEN vp.heartrate   BETWEEN 20 AND 300 THEN vp.heartrate   END) AS hr_max_24h,
  MAX(CASE WHEN vp.respiration BETWEEN 1  AND 100 THEN vp.respiration END) AS rr_max_24h
FROM cov_target t
LEFT JOIN vitalPeriodic vp
  ON  vp.patientunitstayid = t.patientunitstayid
  AND vp.observationoffset >= 0
  AND vp.observationoffset <= t.t24
GROUP BY t.stay_id;

-- ---- GCS 最低值 ----
-- 方案 4.4 要求的是"0–24h **最低** GCS"，即一个**随时间而变的窗内最小值**。
-- 本库有两个 GCS 源，构念并不相同，按"哪个更贴近方案要求"排序取用：
--
--   ① nurseCharting（首选）：逐次记录的 GCS Total，全库 324.8 万行。取 0–24h 内的
--      最小值——这与 MIMIC 侧 chartevents 逐次记录取 24h 最小值的**构念完全一致**。
--      该表在两种 celltype 分类下各存一套 GCS：一路 valname 直接叫 'GCS Total'，
--      另一路 vallabel 为 'Score (Glasgow Coma Scale)' 而数值放在 valname='Value'。
--      两路都要取，只取一路会少掉其中一类（本条不报错、只是覆盖率变低）。
--   ② apacheApsVar（兜底）：每单元一行的**入科时点**GCS，全库仅 17.1 万行，且
--      是单次值而非窗内最小值。仅在 ① 完全没有该患者时使用，并在论文中说明
--      这一小部分人的 GCS 构念与主分析不同（敏感性分析可整段剔除后复核）。
--      该表缺失以 -1 哨兵编码（不是 NULL），故必须显式要求三项都 > 0。
--
-- 无论走哪一路，都限定 3–15 的合法总分区间，剔除录入错误与哨兵值。
CREATE OR REPLACE TABLE cov_gcs AS
SELECT
  t.stay_id,
  COALESCE(MIN(nc.gcs), MIN(ap.gcs)) AS gcs_min_24h
FROM cov_target t
LEFT JOIN (
  SELECT patientunitstayid, nursingchartoffset AS off,
         CASE WHEN TRY_CAST(REPLACE(nursingchartvalue, ' ', '') AS INTEGER) BETWEEN 3 AND 15
              THEN TRY_CAST(REPLACE(nursingchartvalue, ' ', '') AS INTEGER) END AS gcs
  FROM nurseCharting
  WHERE nursingchartcelltypevalname = 'GCS Total'
     OR (nursingchartcelltypevallabel = 'Score (Glasgow Coma Scale)'
         AND nursingchartcelltypevalname = 'Value')
) nc ON nc.patientunitstayid = t.patientunitstayid
    AND nc.off >= 0 AND nc.off <= t.t24
LEFT JOIN (
  SELECT patientunitstayid, eyes + motor + verbal AS gcs
  FROM apache_gcs
  WHERE eyes > 0 AND motor > 0 AND verbal > 0
) ap ON ap.patientunitstayid = t.patientunitstayid
GROUP BY t.stay_id;

-- ---- 体重（供液体出入量按 mL/(kg·h) 标化） ----
-- 主源 patient.admissionweight（本库标准的入科体重，单位为 kg；全库 91.7% 的单元有值）。
--
-- ⚠️ 回退源用 **intakeOutput 的 'Weight|Bodyweight (kg)'**，不是 nurseCharting。
--    实测 nurseCharting 里**根本没有体重条目**（按 '%weight%' 全表筛只返回体温与
--    GCS 两类），体重记在出入量表的 Weight 类目下（87.8 万行）。
-- ⚠️ 该处同样**必须限定 kg**：出入量表的体重条目并存 'Bodyweight (kg)' 与
--    'Bodyweight (lb)' 两条 cellpath，而 cellvaluenumeric 是**不带单位的数值**。
--    若只按"类目含 Weight"取数，一个 180 lb 的患者会被读成 180 kg，落在 20–400 的
--    合理区间内，**既不报错也不产生缺失值**，却让他的尿量速率（S8 的 mL/(kg·h)）
--    整体偏小一半以上，方向性地把 BMI 高的人推离 AKI。
--    故 cellpath 里显式要求 'Bodyweight (kg)'，把 lb 那条排除在外。
--    与 MIMIC 侧只认 226512（Admission Weight (Kg)）/ 224639（Daily Weight）
--    同一思路：**只认单位明确的条目**。
CREATE OR REPLACE TABLE cov_weight AS
SELECT
  t.stay_id,
  COALESCE(
    MIN(CASE WHEN pw.admissionweight BETWEEN 20 AND 400 THEN pw.admissionweight END),
    MIN(CASE WHEN io.w BETWEEN 20 AND 400 THEN io.w END)
  ) AS weight_kg
FROM cov_target t
LEFT JOIN patient pw ON pw.patientunitstayid = t.patientunitstayid
LEFT JOIN (
  SELECT patientunitstayid, TRY_CAST(cellvaluenumeric AS DOUBLE) AS w,
         intakeoutputoffset AS off
  FROM intakeOutput
  WHERE cellpath LIKE '%|Weight|Bodyweight (kg)'
    AND cellvaluenumeric IS NOT NULL
) io ON io.patientunitstayid = t.patientunitstayid
    AND io.off >= 0 AND io.off <= t.t24
GROUP BY t.stay_id;

-- ---- 血管活性药 + 高渗盐水：是否使用 ----
-- infusionDrug 为主源，medication 作补充（记忆已核实两者须取并集：
-- infusionDrug 26,234 单元、medication 21,722 单元且不完全重叠）。
-- drugname 带单位后缀（'Norepinephrine (ml/hr)'、'(mcg/kg/min)' 等），必须子串匹配。
-- ⚠️ 多巴酚丁胺与米力农是**正性肌力药**，不算血管加压药；分开统计，
--    与 MIMIC 侧的 vasopressor_any / inotrope_any 口径一致。
-- ⚠️ 高渗盐水在本库**只能从 infusionDrug 抽到**（medication 与 admissionDrug 全库 0 行），
--    静推/bolus 型基本漏掉，故该变量在 eICU 侧系统性低估，论文中说明。
--
-- ⚠️⚠️ is_hs 的判据**必须要求"高渗"本身的字样，绝不能匹配 'sodium chloride' 或
--    'nacl' 这两个泛词**。本库的输注药名普遍写成"药物 + 溶媒"的全称，
--    例如 'Insulin 250 Units Sodium Chloride 0.9% 250 ml'、
--    'FentaNYL (Sublimaze) 2500 mcg Sodium Chloride 0.9% 250 ml'、
--    'Norepinephrine MAX 32 mg Dextrose 5% 250 ml' —— 溶媒是**等渗**生理盐水，
--    与高渗盐水毫无关系。按泛词匹配会把几乎每一个 ICU 患者都判成"用过高渗盐水"：
--    全列恒为 1，既不报错也不缺值，而 S9（排除高渗盐水者）会因此把整个队列排空。
--    真正的高渗条目只有 '3% Sodium Chloride IVF'、'3% NS'、'3% Saline'、
--    '3% NaCl'、'Hypertonic Saline' 等少数几条（另有用 7.5% / 23.4% 的）。
--    check_covariates() 的"不同取值 ≤ 1"判据能拦住全 1 的结果，但那已是事后补救。
CREATE OR REPLACE TABLE cov_vaso AS
SELECT
  t.stay_id,
  MAX(CASE WHEN w.is_vaso THEN 1 ELSE 0 END)       AS vasoactive_any,
  MAX(CASE WHEN w.is_pressor THEN 1 ELSE 0 END)    AS vasopressor_any,
  MAX(CASE WHEN w.is_inotrope THEN 1 ELSE 0 END)   AS inotrope_any,
  MAX(CASE WHEN w.is_hs THEN 1 ELSE 0 END)         AS hypertonic_saline
FROM cov_target t
LEFT JOIN (
  SELECT patientunitstayid, infusionoffset AS off, drugname,
         LOWER(drugname) LIKE '%norepinephrine%' OR LOWER(drugname) LIKE '%epinephrine%'
      OR LOWER(drugname) LIKE '%dopamine%'       OR LOWER(drugname) LIKE '%vasopressin%'
      OR LOWER(drugname) LIKE '%phenylephrine%'  OR LOWER(drugname) LIKE '%neosynephrine%'
      OR LOWER(drugname) LIKE '%milrinone%'      OR LOWER(drugname) LIKE '%dobutamine%'
                                                          AS is_vaso,
         LOWER(drugname) LIKE '%norepinephrine%' OR LOWER(drugname) LIKE '%epinephrine%'
      OR LOWER(drugname) LIKE '%dopamine%'       OR LOWER(drugname) LIKE '%vasopressin%'
      OR LOWER(drugname) LIKE '%phenylephrine%'  OR LOWER(drugname) LIKE '%neosynephrine%'
                                                          AS is_pressor,
         LOWER(drugname) LIKE '%milrinone%'      OR LOWER(drugname) LIKE '%dobutamine%'
                                                          AS is_inotrope,
         LOWER(drugname) LIKE '%hypertonic%'
      OR LOWER(drugname) LIKE '%3% saline%'   OR LOWER(drugname) LIKE '%3% nacl%'
      OR LOWER(drugname) LIKE '%3% sodium chloride%'
      OR LOWER(drugname) LIKE '%3%ns%'        OR LOWER(drugname) LIKE '%3% ns%'
      OR LOWER(drugname) LIKE '%7.5% saline%' OR LOWER(drugname) LIKE '%23.4%'
                                                          AS is_hs
  FROM infusionDrug WHERE drugname IS NOT NULL
  UNION ALL
  SELECT patientunitstayid, drugstartoffset AS off, drugname,
         LOWER(drugname) LIKE '%norepinephrine%' OR LOWER(drugname) LIKE '%epinephrine%'
      OR LOWER(drugname) LIKE '%dopamine%'       OR LOWER(drugname) LIKE '%vasopressin%'
      OR LOWER(drugname) LIKE '%phenylephrine%'  OR LOWER(drugname) LIKE '%neosynephrine%'
      OR LOWER(drugname) LIKE '%milrinone%'      OR LOWER(drugname) LIKE '%dobutamine%',
         LOWER(drugname) LIKE '%norepinephrine%' OR LOWER(drugname) LIKE '%epinephrine%'
      OR LOWER(drugname) LIKE '%dopamine%'       OR LOWER(drugname) LIKE '%vasopressin%'
      OR LOWER(drugname) LIKE '%phenylephrine%'  OR LOWER(drugname) LIKE '%neosynephrine%',
         LOWER(drugname) LIKE '%milrinone%'      OR LOWER(drugname) LIKE '%dobutamine%',
         LOWER(drugname) LIKE '%hypertonic%'
      OR LOWER(drugname) LIKE '%3% saline%'   OR LOWER(drugname) LIKE '%3% nacl%'
      OR LOWER(drugname) LIKE '%3% sodium chloride%'
      OR LOWER(drugname) LIKE '%3%ns%'        OR LOWER(drugname) LIKE '%3% ns%'
      OR LOWER(drugname) LIKE '%7.5% saline%' OR LOWER(drugname) LIKE '%23.4%'
  -- 已取消的医嘱从未给药，计入会把治疗强度记高（全库 20.5 万行为 TRUE）
  FROM medication WHERE drugname IS NOT NULL
    AND COALESCE(drugordercancelled, FALSE) = FALSE
) w ON w.patientunitstayid = t.patientunitstayid
   AND w.off >= 0 AND w.off <= t.t24
GROUP BY t.stay_id;

-- ---- 血管活性药剂量：去甲肾上腺素当量最大速率 ----
-- 与 MIMIC 侧 cov_vaso_nee **同一套换算**（NE + 肾上腺素 + 多巴胺/100 +
-- 去氧肾上腺素/10 + 血管加压素×2.5），换算前按各药公认临床区间剔除单位/录入错误。
--
-- ⚠️ 两库在"时点"的取法上有一处**已记录的差异**：
--    MIMIC 有输注的 starttime/endtime 区间，在**每个 starttime 上对同时运行的多种药
--    求和**再取时间轴最大值；eICU 的 infusionDrug 没有结束时刻，只有 infusionoffset，
--    故改为**按 infusionoffset 分组求和**再取最大。泵速变更时 eICU 会新写一条记录，
--    因此这一近似与 MIMIC 的口径同向，但不是逐点可比。论文中说明。
--
-- ⚠️ 药名匹配一律"**药名子串 + 单位后缀**"两段式，不用行首精确匹配。
--    本库的药名把药物、浓度、溶媒连写成一条，如
--    'Norepinephrine MAX 32 mg Dextrose 5% 250 ml (mcg/min)'，
--    行首是药名但后面还跟着浓度与溶媒；按 'norepinephrine (mcg/kg/min%' 这种
--    "药名 + 空格 + 单位"的精确式去匹配，会把**全部预混制剂整类漏掉**
--    （实测仅 'Norepinephrine (mcg/min)' 一种就有 19.2 万行，另有 STD/MAX 预混
--    两种、肾上腺素与去氧肾上腺素同类若干）。漏掉不报错，只是把重患者的
--    升压药强度记成 0，方向性地削弱这个协变量。
--
-- ⚠️ 单位三分法，缺一不可：
--      (mcg/kg/min)  —— 已是目标单位，直接用
--      (mcg/min)     —— **须除以体重**换算成 /kg，故 join cov_weight
--      (ml/hr) / ()  —— 无浓度无法换算，与 MIMIC 侧一致地排除
--    体重缺失者的 mcg/min 记录换算不出结果，会被下面的 `nee IS NOT NULL` 自然丢弃，
--    这与"该患者这一条记录不可用"是同义的，不构成额外偏倚。
--
-- ⚠️ 顺序敏感：CASE 自上而下取第一个命中。
--      · '(mcg/kg/min)' 必须排在 '(mcg/min)' 之前；
--      · 'norepinephrine' 必须排在 'epinephrine' 之前——
--        'norepinephrine' 字面上**包含** 'epinephrine'，若先判肾上腺素，
--        所有去甲肾上腺素记录都会被当成肾上腺素（换算系数相同，但两者若将来
--        分开设系数就会静默出错）。
CREATE OR REPLACE TABLE cov_vaso_nee AS
WITH w AS (
  SELECT stay_id, weight_kg FROM cov_weight WHERE weight_kg IS NOT NULL
),
v AS (
  SELECT d.patientunitstayid, d.infusionoffset,
    CASE
      -- (mcg/kg/min)：已是目标单位
      WHEN LOWER(d.drugname) LIKE '%norepinephrine%' AND LOWER(d.drugname) LIKE '%(mcg/kg/min)%'
           AND TRY_CAST(d.drugrate AS DOUBLE) > 0
           AND TRY_CAST(d.drugrate AS DOUBLE) <= 3    THEN TRY_CAST(d.drugrate AS DOUBLE)
      WHEN LOWER(d.drugname) LIKE '%epinephrine%' AND LOWER(d.drugname) LIKE '%(mcg/kg/min)%'
           AND TRY_CAST(d.drugrate AS DOUBLE) > 0
           AND TRY_CAST(d.drugrate AS DOUBLE) <= 3    THEN TRY_CAST(d.drugrate AS DOUBLE)
      WHEN LOWER(d.drugname) LIKE '%dopamine%' AND LOWER(d.drugname) LIKE '%(mcg/kg/min)%'
           AND TRY_CAST(d.drugrate AS DOUBLE) >= 0.1
           AND TRY_CAST(d.drugrate AS DOUBLE) <= 50   THEN TRY_CAST(d.drugrate AS DOUBLE) / 100.0
      WHEN (LOWER(d.drugname) LIKE '%phenylephrine%' OR LOWER(d.drugname) LIKE '%neosynephrine%')
           AND LOWER(d.drugname) LIKE '%(mcg/kg/min)%'
           AND TRY_CAST(d.drugrate AS DOUBLE) >= 0.01
           AND TRY_CAST(d.drugrate AS DOUBLE) <= 10   THEN TRY_CAST(d.drugrate AS DOUBLE) / 10.0
      -- (mcg/min)：按该患者体重换算成 /kg（体重缺失或为 0 时结果为 NULL，自然剔除）
      WHEN LOWER(d.drugname) LIKE '%norepinephrine%' AND LOWER(d.drugname) LIKE '%(mcg/min)%'
           AND TRY_CAST(d.drugrate AS DOUBLE) > 0
           AND TRY_CAST(d.drugrate AS DOUBLE) <= 3 * w.weight_kg
                                                      THEN TRY_CAST(d.drugrate AS DOUBLE) / w.weight_kg
      WHEN LOWER(d.drugname) LIKE '%epinephrine%' AND LOWER(d.drugname) LIKE '%(mcg/min)%'
           AND TRY_CAST(d.drugrate AS DOUBLE) > 0
           AND TRY_CAST(d.drugrate AS DOUBLE) <= 3 * w.weight_kg
                                                      THEN TRY_CAST(d.drugrate AS DOUBLE) / w.weight_kg
      WHEN (LOWER(d.drugname) LIKE '%phenylephrine%' OR LOWER(d.drugname) LIKE '%neosynephrine%')
           AND LOWER(d.drugname) LIKE '%(mcg/min)%'
           AND TRY_CAST(d.drugrate AS DOUBLE) >= 0.01 * w.weight_kg
           AND TRY_CAST(d.drugrate AS DOUBLE) <= 10 * w.weight_kg
                                                      THEN TRY_CAST(d.drugrate AS DOUBLE) / w.weight_kg / 10.0
      -- 血管加压素：两种单位分别换算成 units/min 后再折 NEE
      WHEN LOWER(d.drugname) LIKE '%vasopressin%' AND LOWER(d.drugname) LIKE '%(units/hr)%'
           AND TRY_CAST(d.drugrate AS DOUBLE) >= 0.1
           AND TRY_CAST(d.drugrate AS DOUBLE) <= 10   THEN TRY_CAST(d.drugrate AS DOUBLE) / 60.0 * 2.5
      WHEN LOWER(d.drugname) LIKE '%vasopressin%' AND LOWER(d.drugname) LIKE '%(units/hour)%'
           AND TRY_CAST(d.drugrate AS DOUBLE) >= 0.1
           AND TRY_CAST(d.drugrate AS DOUBLE) <= 10   THEN TRY_CAST(d.drugrate AS DOUBLE) / 60.0 * 2.5
      WHEN LOWER(d.drugname) LIKE '%vasopressin%' AND LOWER(d.drugname) LIKE '%(units/min)%'
           AND TRY_CAST(d.drugrate AS DOUBLE) >= 0.002
           AND TRY_CAST(d.drugrate AS DOUBLE) <= 0.17 THEN TRY_CAST(d.drugrate AS DOUBLE) * 2.5
    END AS nee
  FROM infusionDrug d
  JOIN cov_target t ON t.patientunitstayid = d.patientunitstayid
  LEFT JOIN w ON w.stay_id = t.stay_id
  WHERE d.infusionoffset >= 0 AND d.infusionoffset <= 1440
),
summed AS (
  SELECT patientunitstayid, infusionoffset, SUM(nee) AS nee_t
  FROM v WHERE nee IS NOT NULL GROUP BY patientunitstayid, infusionoffset
)
SELECT t.stay_id, COALESCE(MAX(s.nee_t), 0) AS nee_max
FROM cov_target t
LEFT JOIN summed s ON s.patientunitstayid = t.patientunitstayid
GROUP BY t.stay_id;

-- ---- 利尿剂：袢利尿剂使用与呋塞米当量 ----
-- 主源 medication（口服与静脉都在里面，drugname 子串匹配）。
-- 目标构念与 MIMIC 侧一致：**24 小时呋塞米当量毫克 = 单次剂量 × 每日给药次数**，
-- 换算系数 呋塞米 ×1、布美他尼 ×40。
--
-- ⚠️ 剂量解析**不能用 TRY_CAST，必须用正则取首个数字**。本库 dosage 是"带单位的
--    字符串"，实测形态为 '40 mg' / '40 MG' / '40' / '1 mg'。直接 TRY_CAST 只有
--    **约 2%** 的行（纯数字那种）能转成功，其余静默变成 NULL，
--    求和时被 SUM 忽略——利尿剂当量于是被系统性低估，且不报错、不产生缺失标记。
--    regexp_extract(dosage, '[0-9]+(\.[0-9]+)?', 0) 能覆盖 98.9% 的行。
--
-- ⚠️ 每日次数由 frequency 折算：本库该字段是自由文本（'Once'、'BID'、'Q12H'、
--    'Once X1'、'BID M/E' …），无法直接比较。按"先特异后一般"的顺序映射，
--    未识别的取值保守地记 1 次／日。不折算的话 'BID' 与 'Once' 会被等同看待，
--    把维持治疗的患者记成与单次给药同强度。
--
-- ⚠️ 必须剔除已取消的医嘱（drugordercancelled，全库 20.5 万行为 TRUE）。
--    取消的医嘱从未给药，计入会抬高治疗强度。
-- ⚠️ 本库没有 torasemide/torsemide 记录（全库 0 行），故只用呋塞米与布美他尼。
-- **与 MIMIC 侧的局限一致**：静脉与口服的等效比不同，此处统一按口服当量折算；
-- 该变量只进扩展模型作稳健性分析、不作因果解释（方案 4.4）。
CREATE OR REPLACE TABLE cov_diur AS
SELECT
  t.stay_id,
  -- ⚠️ 必须 COUNT(m.medicationid) 而不是 COUNT(*)：LEFT JOIN 下没有匹配药物的单元
  --    也会留下**一行全 NULL**，COUNT(*) 恒 ≥ 1，于是每个人都被记成"用过利尿剂"。
  --    这类"聚合把左连接的空行也数进去"的错误不报错、结果看着也合理，最难发现。
  CASE WHEN COUNT(m.medicationid) > 0 THEN 1 ELSE 0 END AS diuretic_any,
  COALESCE(SUM(
    TRY_CAST(regexp_extract(m.dosage, '[0-9]+(\.[0-9]+)?', 0) AS DOUBLE)
    * CASE WHEN LOWER(m.drugname) LIKE '%bumetanide%'
             OR LOWER(m.drugname) LIKE '%bumex%'      THEN 40.0 ELSE 1.0 END
    * CASE
        WHEN LOWER(COALESCE(m.frequency, '')) LIKE '%q4h%'
          OR LOWER(COALESCE(m.frequency, '')) LIKE '%every 4 hour%'   THEN 6
        WHEN LOWER(COALESCE(m.frequency, '')) LIKE '%q6h%'
          OR LOWER(COALESCE(m.frequency, '')) LIKE '%qid%'
          OR LOWER(COALESCE(m.frequency, '')) LIKE '%every 6 hour%'   THEN 4
        WHEN LOWER(COALESCE(m.frequency, '')) LIKE '%q8h%'
          OR LOWER(COALESCE(m.frequency, '')) LIKE '%tid%'
          OR LOWER(COALESCE(m.frequency, '')) LIKE '%every 8 hour%'   THEN 3
        WHEN LOWER(COALESCE(m.frequency, '')) LIKE '%q12h%'
          OR LOWER(COALESCE(m.frequency, '')) LIKE '%bid%'
          OR LOWER(COALESCE(m.frequency, '')) LIKE '%every 12 hour%'  THEN 2
        ELSE 1 END
  ), 0) AS diuretic_furosemide_eq_mg
FROM cov_target t
LEFT JOIN medication m
  ON  m.patientunitstayid = t.patientunitstayid
  AND m.drugstartoffset >= 0
  AND m.drugstartoffset <= t.t24
  AND COALESCE(m.drugordercancelled, FALSE) = FALSE
  AND (LOWER(m.drugname) LIKE '%furosemide%' OR LOWER(m.drugname) LIKE '%lasix%'
    OR LOWER(m.drugname) LIKE '%bumetanide%' OR LOWER(m.drugname) LIKE '%bumex%')
GROUP BY t.stay_id;

-- ---- 机械通气 ----
-- 有创：respiratoryCare.airwaytype ∈ {Oral ETT, Nasal ETT, Tracheostomy}（有明确
--       起止时刻），并补 treatment 的 'mechanical ventilation'（覆盖更广但无时刻）。
-- 无创：treatment 的 'non-invasive ventilation' / 'CPAP' / 'BiPAP'。
-- ⚠️ treatment 只有时点没有时长，故一律按"该记录的时刻落在 0–24h 内"计；
--    respiratoryCare 有 ventstart/ventend，按区间与 [0,24h] 重叠计（与 MIMIC 一致）。
CREATE OR REPLACE TABLE cov_vent AS
SELECT
  t.stay_id,
  MAX(CASE WHEN rc.patientunitstayid IS NOT NULL THEN 1
           WHEN tr_inv.patientunitstayid IS NOT NULL THEN 1 ELSE 0 END) AS mech_vent_invasive,
  MAX(CASE WHEN tr_niv.patientunitstayid IS NOT NULL THEN 1 ELSE 0 END) AS mech_vent_noninvasive
FROM cov_target t
LEFT JOIN (
  SELECT DISTINCT patientunitstayid FROM respiratoryCare
  WHERE airwaytype IN ('Oral ETT', 'Nasal ETT', 'Tracheostomy')
    AND COALESCE(ventstartoffset, respcarestatusoffset) <= 1440
    AND COALESCE(ventendoffset, 999999) >= 0
) rc ON rc.patientunitstayid = t.patientunitstayid
LEFT JOIN (
  SELECT DISTINCT patientunitstayid FROM treatment
  WHERE LOWER(treatmentstring) LIKE '%mechanical ventilation%'
    AND LOWER(treatmentstring) NOT LIKE '%non-invasive%'
    AND treatmentoffset >= 0 AND treatmentoffset <= 1440
) tr_inv ON tr_inv.patientunitstayid = t.patientunitstayid
LEFT JOIN (
  SELECT DISTINCT patientunitstayid FROM treatment
  WHERE (LOWER(treatmentstring) LIKE '%non-invasive ventilation%'
      OR LOWER(treatmentstring) LIKE '%cpap%' OR LOWER(treatmentstring) LIKE '%bipap%')
    AND treatmentoffset >= 0 AND treatmentoffset <= 1440
) tr_niv ON tr_niv.patientunitstayid = t.patientunitstayid
GROUP BY t.stay_id;

-- ---- 0–24h 液体入量（mL） ----
-- 取 I&O 入量区的 cellvaluenumeric 求和。**cellpath 前缀必须限定在 Intake**，
-- 否则会把出量条目一并算进来。
CREATE OR REPLACE TABLE cov_fluid_in AS
SELECT
  t.stay_id,
  COALESCE(SUM(io.cellvaluenumeric), 0) AS fluid_in_24h_ml
FROM cov_target t
LEFT JOIN intakeOutput io
  ON  io.patientunitstayid = t.patientunitstayid
  AND io.intakeoutputoffset >= 0
  AND io.intakeoutputoffset <= t.t24
  AND io.cellvaluenumeric IS NOT NULL
  AND io.cellpath LIKE 'flowsheet|Flowsheet Cell Labels|I&O|Intake%'
GROUP BY t.stay_id;

-- ---- 0–24h 尿量（mL） ----
-- 与 21_urine_eicu.sql 用**同一套标签判据**（含 urine、且排除 count/occurrence/
-- incontinen/stool 四类非计量条目），保证协变量层的尿量与 S8 结局判定同源。
CREATE OR REPLACE TABLE cov_urine_out AS
SELECT
  t.stay_id,
  COALESCE(SUM(io.cellvaluenumeric), 0) AS urine_out_24h_ml
FROM cov_target t
LEFT JOIN intakeOutput io
  ON  io.patientunitstayid = t.patientunitstayid
  AND io.intakeoutputoffset >= 0
  AND io.intakeoutputoffset <= t.t24
  AND io.cellvaluenumeric IS NOT NULL
  AND io.cellpath LIKE 'flowsheet|Flowsheet Cell Labels|I&O|Output (ml)|%'
  AND LOWER(TRIM(io.celllabel)) LIKE '%urine%'
  AND LOWER(TRIM(io.celllabel)) NOT LIKE '%count%'
  AND LOWER(TRIM(io.celllabel)) NOT LIKE '%occurrence%'
  AND LOWER(TRIM(io.celllabel)) NOT LIKE '%incontinen%'
  AND LOWER(TRIM(io.celllabel)) NOT LIKE '%stool%'
GROUP BY t.stay_id;

-- ---- 0–24h 血糖与血红蛋白 ----
-- 血糖 labname ∈ {glucose, bedside glucose}；血红蛋白 labname ∈ {hgb, hemoglobin}。
-- 与 MIMIC 侧同样按公认临床区间限幅（方案 4.1.3）。
CREATE OR REPLACE TABLE cov_labs AS
SELECT
  t.stay_id,
  MAX(CASE WHEN LOWER(l.labname) IN ('glucose', 'bedside glucose')
            AND l.labresult BETWEEN 20 AND 1000 THEN l.labresult END) AS glucose_max_24h,
  MIN(CASE WHEN LOWER(l.labname) IN ('glucose', 'bedside glucose')
            AND l.labresult BETWEEN 20 AND 1000 THEN l.labresult END) AS glucose_min_24h,
  MIN(CASE WHEN LOWER(l.labname) IN ('hgb', 'hemoglobin', 'haemoglobin')
            AND l.labresult BETWEEN 3 AND 25 THEN l.labresult END)   AS hgb_min_24h
FROM cov_target t
LEFT JOIN lab l
  ON  l.patientunitstayid = t.patientunitstayid
  AND l.labresultoffset >= 0
  AND l.labresultoffset <= t.t24
  AND l.labresult IS NOT NULL
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
  v.sbp_min_24h, v2.hr_min_24h, v2.hr_max_24h, v2.rr_max_24h,
  g.gcs_min_24h,
  w.weight_kg,
  -- 治疗强度
  m.vasoactive_any, m.vasopressor_any, m.inotrope_any, n.nee_max,
  m.hypertonic_saline,
  -- 利尿剂
  d.diuretic_any, d.diuretic_furosemide_eq_mg,
  -- 机械通气
  ve.mech_vent_invasive, ve.mech_vent_noninvasive,
  CASE WHEN COALESCE(ve.mech_vent_invasive, 0) = 1
         OR COALESCE(ve.mech_vent_noninvasive, 0) = 1 THEN 1 ELSE 0 END AS mech_vent_any,
  -- 出入量
  fi.fluid_in_24h_ml,
  uo.urine_out_24h_ml,
  fi.fluid_in_24h_ml - uo.urine_out_24h_ml AS fluid_balance_24h_ml,
  -- 化验
  lb.glucose_max_24h, lb.glucose_min_24h, lb.hgb_min_24h
FROM cov_target b
LEFT JOIN cov_icd       i  ON i.hadm_id = b.hadm_id
LEFT JOIN cov_vitals    v  USING (stay_id)
LEFT JOIN cov_vitals2   v2 USING (stay_id)
LEFT JOIN cov_gcs       g  USING (stay_id)
LEFT JOIN cov_weight    w  USING (stay_id)
LEFT JOIN cov_vaso      m  USING (stay_id)
LEFT JOIN cov_vaso_nee  n  USING (stay_id)
LEFT JOIN cov_diur      d  USING (stay_id)
LEFT JOIN cov_vent      ve USING (stay_id)
LEFT JOIN cov_fluid_in  fi USING (stay_id)
LEFT JOIN cov_urine_out uo USING (stay_id)
LEFT JOIN cov_labs      lb USING (stay_id);
