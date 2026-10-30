-- =============================================================================
-- 42_rrt_eicu.sql —— 肾脏替代治疗事件（eICU 侧）
--
-- rrt 长表 (id, t_hours, value = 1)，t_hours 相对 ICU 入科小时数。用途与 MIMIC
-- 侧完全一致：
--   (a) 排除标准③：界标前（t ≤ 24h）已开始 RRT 者不进分析
--   (b) 次要结局：界标后（24 < t ≤ 72h）启动 RRT 者按 KDIGO 3 期计
--
-- ── 本库 RRT 的主源是 treatment，不是诊断码 ──────────────────────────────────
-- 记忆中已核实：diagnosis 的 V45.1 / V56 / Z99.2 / Z49 在本库**全库 0 行**，
-- 靠诊断码找 RRT 会全军覆没。treatment.treatmentstring 才是主源
-- （含 dialysis 者 8,298 单元），且有 treatmentoffset 可定位时点。
--
-- ── 必须排除的三类"看起来像 RRT 但不是 RRT"的记录 ────────────────────────────
--   insertion of venous catheter for hemodialysis  —— 置管操作（2,292 行）
--   insertion of venous catheter ... percutaneous / tunneled —— 同上
--   dialysis access surgery                        —— 造瘘
-- 这些是**为透析做准备**，不代表该时点已在接受 RRT；若计入，会把"计划透析"
-- 的患者整段算成"已在 RRT"，排除标准③会过度排除。
--
-- 另有 'renal|electrolyte correction|treatment of hyperkalemia|dialysis'（260 行）
-- 属治疗高钾的透析，**是** RRT，保留。
-- =============================================================================

CREATE OR REPLACE TABLE rrt AS
SELECT DISTINCT
  CAST(t.patientunitstayid AS VARCHAR)      AS id,
  t.treatmentoffset / 60.0                  AS t_hours,
  1                                         AS value
FROM treatment t
WHERE (
        LOWER(t.treatmentstring) LIKE '%dialysis%'
     OR LOWER(t.treatmentstring) LIKE '%c v v h%'          -- 'C V V H' / 'C V V H D'
     OR LOWER(t.treatmentstring) LIKE '%renal replacement%'
     OR LOWER(t.treatmentstring) LIKE '%ultrafiltration%'
      )
  AND LOWER(t.treatmentstring) NOT LIKE '%insertion of venous catheter%'
  AND LOWER(t.treatmentstring) NOT LIKE '%dialysis access surgery%'
  AND LOWER(t.treatmentstring) NOT LIKE '%dialysis catheter%'
  -- 'renal|dialysis|arteriovenous shunt for renal dialysis'（37 单元）同样是**造瘘**，
  -- 名里带 dialysis 只是说明用途，不代表该时点已在接受 RRT。不排掉会把造瘘当天
  -- 当成"已开始 RRT"，让这部分人被排除标准③误删。
  AND LOWER(t.treatmentstring) NOT LIKE '%arteriovenous shunt%'
  AND t.treatmentoffset IS NOT NULL;
