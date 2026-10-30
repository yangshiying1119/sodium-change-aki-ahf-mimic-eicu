## =============================================================================
## i18n.R —— 中英文两套呈现之间的唯一开关
##
## 设计原则：**数字只算一次，语言只换一层**。
## 表与图的全部数值仍由原来的模块产生，本文件不碰任何统计；
## 它只在落盘前把中文标签换成英文，因此中英文两版的数字不可能不一致。
##
## 用法：环境变量 OUT_LANG=en 时输出英文，否则（缺省）输出中文。
##   OUT_LANG=en Rscript _run6.R      → out/tables_en/
##   OUT_LANG=en Rscript _run5.R      → out/figures_en/
## 中文那一套完全不受影响：所有函数在 zh 模式下原样返回。
## =============================================================================

OUT_LANG <- tolower(Sys.getenv("OUT_LANG", "zh"))
is_en <- function() identical(OUT_LANG, "en")

## 英文图的字体：期刊通行 Arial，本机用其等宽度替代字 Liberation Sans。
## 中文那套仍用 Noto Sans SC（figures.R 里按语言选，见 FONT 的定义）。
FONT_EN <- "Liberation Sans"

## ---- 词表 --------------------------------------------------------------------
## 键是**生成时就已确定的字符串**（含全角标点），值是对应的英文。
## 不在表内的字符串原样返回，由 verify_en.R 扫描出来逐个补——不静默吞掉。
EN_TERMS <- c(
  ## 表头与结构性列名
  "变量" = "Variable", "整体" = "Overall", "缺失" = "Missing",
  "缺失数" = "No. missing", "模型" = "Model", "分析" = "Analysis",
  "常数项" = "(Intercept)",
  "人数" = "No. of patients", "事件数" = "Events", "亚组" = "Subgroup",
  "层" = "Level", "层内人数" = "No. in stratum", "值" = "Value",
  "项目" = "Item", "节点" = "Step", "排除原因" = "Reason for exclusion",
  "本层排除" = "Excluded, n", "留存" = "Remaining, n", "状态" = "Status",
  "累积发生率" = "Cumulative incidence",
  "累积发生率（72 h）" = "Cumulative incidence at 72 h",
  "备注" = "Note", "合计" = "Total", "上限" = "Upper limit",
  "下限" = "Lower limit", "交互 P" = "P for interaction",
  "分析人群" = "Analysis population", "界标后未复测" = "No repeat test after landmark",
  "可判定结局者" = "Evaluable for the outcome",
  "结局不可判定者" = "Outcome not evaluable",

  ## ΔNa 三档分组
  "ΔNa ≤ −3（下降）" = "ΔNa ≤ −3 (decrease)",
  "−3 < ΔNa < 3（稳定）" = "−3 < ΔNa < 3 (stable)",
  "ΔNa ≥ 3（升高）" = "ΔNa ≥ 3 (increase)",

  ## 基线表变量
  "年龄（岁）" = "Age, years", "男性" = "Male",
  "ICU 类型" = "ICU type", "入院类型" = "Admission type",
  "首次血钠（mmol/L）" = "First serum sodium (mmol/L)",
  "末次血钠（mmol/L）" = "Last serum sodium (mmol/L)",
  "ΔNa（mmol/L）" = "ΔNa (mmol/L)",
  "血钠检测次数" = "No. of sodium measurements",
  "0–24h 血钠检测次数" = "No. of sodium measurements, 0–24 h",
  "界标前最低肌酐（mg/dL）" = "Lowest creatinine before landmark (mg/dL)",
  "界标前肌酐检测次数" = "No. of creatinine measurements before landmark",
  "界标后肌酐检测次数" = "No. of creatinine measurements after landmark",
  "首末血钠间隔（h）" = "Interval between first and last sodium (h)",
  "慢性肾病" = "Chronic kidney disease", "糖尿病" = "Diabetes",
  "高血压" = "Hypertension", "冠心病" = "Coronary heart disease",
  "心房颤动" = "Atrial fibrillation", "心脏瓣膜病" = "Valvular heart disease",
  "0–24h 最低收缩压（mmHg）" = "Lowest systolic BP, 0–24 h (mmHg)",
  "0–24h 最高心率（次/分）" = "Highest heart rate, 0–24 h (bpm)",
  "0–24h 最高呼吸频率（次/分）" = "Highest respiratory rate, 0–24 h (breaths/min)",
  "0–24h 最低 GCS（分）" = "Lowest GCS, 0–24 h",
  "0–24h 最高血糖（mg/dL）" = "Highest glucose, 0–24 h (mg/dL)",
  "0–24h 最低血红蛋白（g/dL）" = "Lowest haemoglobin, 0–24 h (g/dL)",
  "0–24h 机械通气" = "Mechanical ventilation, 0–24 h",
  "0–24h 血管活性药" = "Vasoactive agents, 0–24 h",
  ## 这两条的层记号由 .tr1 的"，n（%）"分支统一补，词条里不再自带 "└ "
  "└ 其中血管加压药" = "Vasopressors",
  "└ 其中强心药" = "Inotropes",
  "0–24h 去甲肾上腺素当量最大速率，仅使用者（mcg/kg/min）" =
    "Maximum norepinephrine-equivalent rate, 0–24 h, users only (mcg/kg/min)",
  "0–24h 利尿剂" = "Diuretics, 0–24 h",
  "0–24h 利尿剂（呋塞米当量 mg）" =
    "Diuretic dose, 0–24 h (furosemide-equivalent mg)",
  "0–24h 液体入量（mL）" = "Fluid intake, 0–24 h (mL)",
  "0–24h 尿量（mL）" = "Urine output, 0–24 h (mL)",

  ## ICU 类型与入院类型的水平
  "冠脉监护病房（CCU）" = "Coronary care unit (CCU)",
  "内科监护病房（MICU）" = "Medical ICU (MICU)",
  "心血管外科监护病房（CVICU）" = "Cardiac vascular ICU (CVICU)",
  "内外科监护病房（MICU/SICU）" = "Medical/surgical ICU (MICU/SICU)",
  "外科监护病房（SICU）" = "Surgical ICU (SICU)",
  "创伤外科监护病房（TSICU）" = "Trauma SICU (TSICU)",
  "神经科中间监护" = "Neuro intermediate care",
  "神经外科监护病房" = "Neurosurgical ICU",
  "神经科降阶梯病房" = "Neuro stepdown",
  "麻醉后恢复室（PACU）" = "Post-anaesthesia care unit (PACU)",
  "急诊" = "Emergency", "紧急" = "Urgent", "观察入院" = "Observation",
  "直接急诊入院" = "Direct emergency",
  "手术当日入院" = "Surgical same-day admission", "择期" = "Elective",
  "急诊观察" = "Emergency observation",

  ## 三档模型
  "未调整" = "Unadjusted", "基础" = "Base", "扩展" = "Extended",

  ## 表 3c 判别力、校准与临床量级
  "范畴" = "Domain", "指标" = "Metric", "结果" = "Result",
  "判别力" = "Discrimination", "校准" = "Calibration",
  "临床量级" = "Clinical magnitude",
  "c 统计量（表观）" = "c-statistic (apparent)",
  "c 统计量（95% CI）" = "c-statistic (95% CI)",
  "去掉 ΔNa 后的 c 统计量" = "c-statistic without ΔNa",
  "ΔAUC（DeLong 检验）" = "ΔAUC (DeLong test)",
  "Hosmer-Lemeshow 检验" = "Hosmer–Lemeshow test",
  "校准斜率（bootstrap 校正）" = "Calibration slope (bootstrap-corrected)",
  "ΔNa 每升高 3 mmol/L 的 OR" = "OR per 3 mmol/L higher ΔNa",
  "似然比检验（ΔNa）" = "Likelihood ratio test (ΔNa)",

  ## 表 3d 校准曲线
  "十分位" = "Decile",
  "平均预测概率" = "Mean predicted probability",
  "实测发生率" = "Observed incidence", "差值" = "Difference",

  ## 表 14 两库并排
  "主队列 n" = "Primary cohort n",
  "主队列 OR（95% CI）" = "Primary cohort OR (95% CI)",
  "重复队列 n" = "Replication cohort n",
  "重复队列 OR（95% CI）" = "Replication cohort OR (95% CI)",

  ## 表 16 严格表型
  "口径" = "Definition",
  "未调整 OR（95% CI）" = "Unadjusted OR (95% CI)",
  "基础模型 OR（95% CI）" = "Base model OR (95% CI)",

  ## 亚组
  "首次血钠 <135 vs ≥135" = "First sodium <135 vs ≥135",
  "有无慢性肾病" = "Chronic kidney disease",
  "首次血钠 <135" = "First sodium <135",
  "首次血钠 ≥135" = "First sodium ≥135",
  "有慢性肾病" = "With CKD", "无慢性肾病" = "Without CKD",

  ## 敏感性分析
  "主分析" = "Main analysis",
  "主分析（表型 B）" = "Main analysis (phenotype B)",
  "主分析（表型 C）" = "Main analysis (phenotype C)",
  "S1 最早两次血钠之差" = "S1 ΔNa from the first two sodium measurements",
  "S2 暴露窗改为 0–18h" = "S2 Exposure window 0–18 h",
  "S3 剔除 |ΔNa| > 12 mmol/L" = "S3 Excluding |ΔNa| > 12 mmol/L",
  "S4 血气钠（50824）替换生化钠" =
    "S4 Blood-gas sodium (itemid 50824) instead of chemistry sodium",
  "S4b 主分析限制到有 ≥2 次血气钠者（仍用生化钠）" =
    "S4b Restricted to patients with ≥2 blood-gas sodium values (chemistry sodium still used)",
  "S5 血糖校正血钠" = "S5 Glucose-corrected sodium",
  "S6 表型定义 A（严格）" = "S6 Strict phenotype definition A",
  "S7 仅保留 48h ≥0.3 mg/dL 标准" = "S7 48-h ≥0.3 mg/dL criterion only",
  "S8 复合口径（肌酐 ∪ 尿量）" =
    "S8 Composite outcome (creatinine or urine output)",
  "S8 复合口径（肌酐 ∪ 尿量，表型 C）" =
    "S8 Composite outcome (creatinine or urine output), phenotype C",
  "S9 排除高渗盐水或 RRT 者" =
    "S9 Excluding patients given hypertonic saline or RRT",
  "S10 不限制首次住院" = "S10 No restriction to the first admission",
  "附加 调整血钠检测次数" =
    "Additional: adjusted for the number of sodium measurements",
  "附加 仅保留界标后留观 ≥48h 者" =
    "Additional: restricted to patients observed ≥48 h after the landmark",
  "附加 多重插补（基础模型，m=20）" =
    "Additional: multiple imputation (base model, m = 20)",
  "附加 多重插补（扩展模型，m=20）" =
    "Additional: multiple imputation (extended model, m = 20)",
  "附加 多重插补（mice）" = "Additional: multiple imputation (mice)",
  "Rubin 规则合并" = "Pooled by Rubin's rules",
  "复测概率加权（IPW）" = "Inverse probability of retesting weighting (IPW)",
  "权重截断 1%–99%。加权到“除复测外全部条件均满足”的目标人群" =
    "Weights trimmed at 1%–99%; weighted to the population meeting all inclusion criteria except retesting",

  ## 逐层留存
  "N0 数据库内全部 ICU 住院" = "N0 All ICU admissions in the database",
  "N1 成人（≥18 岁）" = "N1 Adults (≥18 years)",
  "N4 界标时仍在 ICU 且存活 — 24 h 内死亡" =
    "N4 In ICU and alive at landmark — death within 24 h",
  "N5 0–24 h 内 ≥2 次血钠" = "N5 ≥2 sodium measurements within 0–24 h",
  "N6 界标前可判定为无 AKI — 已有 AKI" =
    "N6 No AKI at baseline — AKI already present",
  "N7 界标后有可判定的肌酐比较值" =
    "N7 Evaluable creatinine comparator after landmark",
  "N8 结局窗内检出 AKI" = "N8 AKI detected within the outcome window",
  "最终分析人群" = "Final analysis population",
  "24 h 内转出" = "Transferred out within 24 h",
  "24 h 内死亡" = "Death within 24 h",
  "└ 24 h 内转出 ICU" = "└ Transferred out of ICU within 24 h",
  "└ 因非首次符合条件的住院而排除" = "└ Not the first eligible admission",
  "└ 排除①a 诊断码（N18.6/585.6/Z99.2/V45.1x）" =
    "└ Exclusion 1a: diagnosis code (N18.6/585.6/Z99.2/V45.1x)",
  "└ 排除①a 诊断码（N186/5856/Z992/V4511/V4512）" =
    "└ Exclusion 1a: diagnosis code (N186/5856/Z992/V4511/V4512)",
  "└ 排除①b 入科问诊'透析患者'（chartevents 225126）" =
    "└ Exclusion 1b: 'dialysis patient' in admission history (chartevents 225126)",
  "└ 排除①b 慢性肾衰透析治疗（treatment）" =
    "└ Exclusion 1b: chronic renal failure on dialysis (treatment)",
  "└ 排除③ 界标前已开始 RRT" = "└ Exclusion 3: RRT started before the landmark",
  "└ 界标前无法判定（肌酐 < 2 次）" =
    "└ Not evaluable before landmark (<2 creatinine values)",
  "不符合表型" = "Does not meet phenotype",
  "未满 18 岁" = "Aged <18 years",

  ## 偏倚分析
  "界标前已有 AKI" = "AKI present before landmark",
  "界标前不满足随访条件" = "Did not meet follow-up criteria before landmark",
  "界标前无法判定" = "Not evaluable before landmark",
  "界标前已行 RRT" = "RRT before landmark",
  "界标前已满足 KDIGO" = "Met KDIGO criteria before landmark",
  "缺乏比较基线" = "No comparator available",
  "终末期肾病/透析依赖" = "End-stage kidney disease / dialysis-dependent",
  "透析依赖" = "Dialysis-dependent",
  "非首次住院" = "Not the first admission",
  "血钠不足 2 次" = "Fewer than 2 sodium measurements",
  "界标后无肌酐复测" = "No creatinine retest after landmark",
  "界标后无肌酐复测以外的其他原因" =
    "Reasons other than a missing repeat creatinine",
  "界标后复测比例" = "Retested after landmark",
  "纳入分析比例" = "Included in analysis",
  "AKI 检出率" = "AKI detection rate",
  "2 次" = "2", "3 次" = "3", "4–5 次" = "4–5", "≥6 次" = "≥6",
  "界标后是否复测肌酐" = "Creatinine retested after landmark",
  "界标前是否已少尿" = "Oliguria before landmark",
  "肌酐标准" = "Creatinine criterion",
  "肌酐标准 ∪ 尿量标准" = "Creatinine or urine output criterion",
  "尿量标准检出" = "Detected by urine output criterion",
  "尿量标准未检出" = "Not detected by urine output criterion",
  "未检出" = "Not detected", "检出" = "Detected",
  "结局窗内 AKI" = "AKI within the outcome window",

  ## 剂量反应表
  "整体关联检验 P" = "P for overall association",
  "非线性检验 P" = "P for non-linearity",
  "节点位置（mmol/L）" = "Knot positions (mmol/L)",
  "参考点（mmol/L）" = "Reference point (mmol/L)",
  "0（参考）" = "0 (reference)", "1.00（参考）" = "1.00 (reference)",

  ## 图内文字
  "图 1  研究人群逐层留存" = "Figure 1. Flow of the study population",
  "图 2  早期血钠变化与 AKI 的剂量反应关系" =
    "Figure 2. Dose–response relationship between early change in serum sodium and AKI",
  "图 3 敏感性分析：ΔNa 每升高 1 mmol/L 的 OR" =
    "Figure 3. Sensitivity analyses: OR per 1 mmol/L higher ΔNa",
  "图 4 预设亚组分析（探索性）" =
    "Figure 4. Prespecified subgroup analyses (exploratory)",
  "图 5  界标后 AKI 与竞争事件的累积发生率" =
    "Figure 5. Cumulative incidence of AKI and of the competing event after the landmark",
  "图 6  早期血钠变化与 AKI 的剂量反应关系：两库对比" =
    "Figure 6. Dose–response relationship: MIMIC-IV versus eICU-CRD",
  "图 7  两库三档模型的效应量并排" =
    "Figure 7. Effect estimates across model tiers in both databases",
  "AKI 的比值比（95% CI）" = "Odds ratio for AKI (95% CI)",
  "累积发生率" = "Cumulative incidence (%)",
  "ICU 入科后小时数" = "Hours after ICU admission",
  "参考 ΔNa = 0" = "Reference ΔNa = 0",
  "整体 P = %s" = "Overall P = %s",
  "非线性 P = %s" = "Non-linear P = %s",
  "OR（95% CI）" = "OR (95% CI)",
  "比值比（OR，对数刻度）" = "Odds ratio (log scale)",
  "AKI（结局）" = "AKI (outcome)",
  "转出或死亡（竞争事件）" = "Discharge or death (competing event)",
  "转出或死亡" = "Discharge or death",
  "MIMIC-IV 主队列" = "MIMIC-IV primary cohort",
  "eICU-CRD 重复队列" = "eICU-CRD replication cohort",
  "亚组（各层独立建模）" = "Subgroup (model fitted within each level)",
  "各层内均为基础模型调整。" = "All estimates adjusted with the base model.",
  "全部采用基础模型调整；\n红色方块为主分析。" =
    "All estimates adjusted with the base model;\nred square denotes the main analysis.",
  ## 图 7 的脚注只能是单行，且从表格左边缘起排——写得比左列长就会横穿到刻度
  ## 标签上（中文原句 15 字正好卡住）。英文取"同一套分析代码"的短表述，完整的
  ## "判定与建模同源"放在图题与正文 6.4 节里说。
  "两库共用同一套判定与建模代码。" = "Same analysis code in both databases.",
  "模型（ΔNa 每 +1 mmol/L）" = "Model (OR per +1 mmol/L ΔNa)",
  ## eICU 侧的两条注脚：带表型字母的题名走下面的 EN_RULES
  "全部采用基础模型调整；红色方块为主分析。\n本库无血气钠，故无 S4/S4b。" =
    paste0("All estimates adjusted with the base model; red square denotes the main analysis.\n",
           "This database has no blood-gas sodium, so S4/S4b are not applicable."),

  ## 表 17 严重度梯度
  "结局" = "Outcome",
  "任何 AKI（≥1 期）" = "Any AKI (stage ≥1)",
  "中重度 AKI（≥2 期）" = "Moderate-to-severe AKI (stage ≥2)",
  "重度 AKI（3 期）" = "Severe AKI (stage 3)",
  "有序结局（无 / 1 期 / 2 期 / 3 期，比例优势）" =
    "Ordinal outcome (none / stage 1 / 2 / 3, proportional odds)",

  ## 表 18 绝对风险差与 E 值
  "对比" = "Contrast", "量" = "Metric", "估计值" = "Estimate",
  "ΔNa 每升高 1 mmol/L" = "Per 1 mmol/L higher ΔNa",
  "ΔNa 每升高 3 mmol/L" = "Per 3 mmol/L higher ΔNa",
  "OR（基础模型）" = "OR (base model)",
  "OR（基础模型，每 1 mmol/L）" = "OR (base model, per 1 mmol/L)",
  "OR（基础模型，每 3 mmol/L）" = "OR (base model, per 3 mmol/L)",
  "E 值（每 1 mmol/L）" = "E-value (per 1 mmol/L)",
  "E 值（每 3 mmol/L）" = "E-value (per 3 mmol/L)",
  "换算为近似 RR" = "Approximate RR",
  "E 值（点估计）" = "E-value (point estimate)",
  "E 值（95% CI 端点）" = "E-value (95% CI limit)",
  "标准化 72 h AKI 风险（ΔNa = 0）" = "Standardised 72-h AKI risk (ΔNa = 0)",
  "标准化 72 h AKI 风险（ΔNa = +3）" = "Standardised 72-h AKI risk (ΔNa = +3)",
  "绝对风险差（+3 vs 0 mmol/L）" = "Absolute risk difference (+3 vs 0 mmol/L)",
  "绝对风险差（+3 vs −3 mmol/L）" = "Absolute risk difference (+3 vs −3 mmol/L)",

  ## 表 19 阴性对照
  "暴露" = "Exposure",
  "ΔNa" = "ΔNa", "ΔK" = "ΔK", "ΔCl" = "ΔCl",
  "ΔNa（同时调整 ΔK）" = "ΔNa (mutually adjusted for ΔK)",
  "ΔNa（同时调整 ΔCl）" = "ΔNa (mutually adjusted for ΔCl)",
  "ΔCl（同时调整 ΔNa）" = "ΔCl (mutually adjusted for ΔNa)",
  "ΔNa（每 1 SD，同时调整 ΔCl）" = "ΔNa per 1 SD (mutually adjusted for ΔCl)",
  "ΔCl（每 1 SD，同时调整 ΔNa）" = "ΔCl per 1 SD (mutually adjusted for ΔNa)",
  "阴性对照暴露（0–24 h 末次减首次，基础模型）" =
    "Negative-control exposures (last minus first over 0–24 h, base model)"
)

## 表格里的 P 值在中文源串里写成 "P 0.125" / "P <0.001"（等号省、不等号紧贴数字）。
## 英文要补 "=" 并把不等号与数字分开，否则直接替换会得到 "P = <0.001"。
.p_en <- function(p) {
  p <- trimws(p)
  if (grepl("^[<≤>≥]", p)) paste0("P ", sub("^([<≤>≥])\\s*", "\\1 ", p))
  else paste0("P = ", p)
}

## ---- 带数字的串用规则处理 ----------------------------------------------------
## 这些串里的数字是运行时生成的，无法逐个写进词表，故用正则 + 替换函数。
## 每条规则只匹配**完整串**（^…$），避免误伤正文中的正常句子。
EN_RULES <- list(
  list(p = "^样本过小（n=([0-9]+)，事件=([0-9]+)），不足以拟合基础模型$",
       f = function(n, e) sprintf("Too few patients (n = %s, events = %s) to fit the base model", n, e)),
  list(p = "^被复测者的样本过小（n=([0-9]+)，事件=([0-9]+)），不足以拟合加权模型$",
       f = function(n, e) sprintf("Too few retested patients (n = %s, events = %s) to fit the weighted model", n, e)),
  list(p = "^其中仅尿量标准检出 ([0-9,]+) 例$",
       f = function(n) sprintf("Detected by the urine output criterion only: %s", n)),
  list(p = "^交互 ([0-9.]+)$", f = function(p) paste0("Interaction P = ", p)),
  list(p = "^N2 急性心衰表型（定义 ([A-C])）$",
       f = function(cc) sprintf("N2 Acute heart failure phenotype (definition %s)", cc)),
  list(p = "^N2 急性心衰表型 ([A-C])$",
       f = function(cc) sprintf("N2 Acute heart failure phenotype %s", cc)),
  list(p = "^附加 多重插补（([^）]*)）$",
       f = function(m) sprintf("Additional: multiple imputation (%s)", m)),
  ## 流程图方框里的"2,700 例"与"…；结局窗内检出 AKI\n… 例（26.6%）"
  list(p = "^([0-9,]+) 例$", f = function(n) paste0("n = ", n)),
  ## 最终人群框宽只有画布的四成多，中文 32 个"西文宽度"刚满；英文同类信息直译要
  ## 54 个，实测整行溢出框外。改用 the outcome window 的紧凑写法（40 个），
  ## 语义不变——框上是给出结局窗内的检出数，窗口本身在正文与图题里已交代。
  list(p = "^([0-9,]+) 例；结局窗内检出 AKI\n([0-9,]+) 例（([0-9.]+)%）$",
       f = function(n, e, p) sprintf("n = %s\nAKI in the outcome window: %s (%s%%)", n, e, p)),
  list(p = "^排除 ([0-9,]+) 例（余 ([0-9,]+) 例）$",
       f = function(n, l) sprintf("Excluded %s (remaining %s)", n, l)),
  ## eICU 侧题名与图例带表型字母（B / C），只能按规则拼
  list(p = "^敏感性分析：ΔNa 每升高 1 mmol/L 的 OR（eICU，表型 ([A-C])）$",
       f = function(cc) sprintf("Sensitivity analyses: OR per 1 mmol/L higher ΔNa (eICU, phenotype %s)", cc)),
  list(p = "^预设亚组分析（eICU，探索性）$",
       f = function() "Prespecified subgroup analyses (eICU, exploratory)"),
  list(p = "^MIMIC-IV 主队列（表型 ([A-C])）$",
       f = function(cc) sprintf("MIMIC-IV primary cohort (phenotype %s)", cc)),
  list(p = "^eICU-CRD 重复队列（表型 ([A-C])）$",
       f = function(cc) sprintf("eICU-CRD replication cohort (phenotype %s)", cc)),
  list(p = "^风险集 n = ([0-9,]+)；界标后 72 h 内 AKI ([0-9,]+) 例，转出或死亡 ([0-9,]+) 例$",
       f = function(n, a, d) sprintf("Risk set n = %s; AKI %s and discharge or death %s within 72 h after the landmark", n, a, d)),
  ## 表 3c/3d/14/16 的"结果"列与表头：数字本身是现成的呈现串，只有括号里的
  ## 中文说明要换，故按数字 + 固定结构写规则，不把它们拆开重排。
  list(p = "^([0-9.]+)（([0-9,]+) 例，([0-9,]+) 例事件）$",
       f = function(r, n, e) sprintf("%s (n = %s; %s events)", r, n, e)),
  list(p = "^c 统计量（bootstrap 校正，B=([0-9]+)）$",
       f = function(b) sprintf("c-statistic (bootstrap-corrected, B = %s)", b)),
  list(p = "^([0-9.]+)（P (.+)）$",
       f = function(v, p) sprintf("%s (%s)", v, .p_en(p))),
  list(p = "^χ² = (.+)，自由度 ([0-9]+)，P (.+)$",
       f = function(x, df, p) sprintf("χ² = %s, df = %s, %s", x, df, .p_en(p))),
  list(p = "^χ² = (.+)，P (.+)$",
       f = function(x, p) sprintf("χ² = %s, %s", x, .p_en(p))),
  list(p = "^原表型定义（表型 ([A-C])）$",
       f = function(cc) sprintf("Original phenotype (phenotype %s)", cc)),
  list(p = "^严格表型（定义 ([A-C])）$",
       f = function(cc) sprintf("Strict phenotype (definition %s)", cc)),
  list(p = "^CI 端点 ([0-9.]+)$",
       f = function(v) sprintf("CI limit %s", v))
)

## 行标签的通用后缀：多水平/二分类变量统一标注 n（%）
.tr_suffix <- function(s) {
  if (grepl("，n（%）$", s)) {
    base <- sub("，n（%）$", "", s)
    indent <- if (startsWith(base, "└ ")) "  " else ""
    base <- sub("^└ ", "", base)
    return(paste0(indent, "└ ", .tr1(base), ", n (%)"))
  }
  s
}

.tr1 <- function(s) {
  if (is.na(s) || !nzchar(s)) return(s)
  ## 表格里的层级靠缩进表达，而缩进是全角空格（"　冠脉监护病房（CCU）"、
  ## "　└ 其中血管加压药"）。全角空格只承担排版，本身无信息；若让它在查询前
  ## 留在串里，词典就得为每一层缩进各存一份词条，多一层缩进就静默漏翻一次。
  ## 所以先剥掉再查，译完按同一层级补两个半角空格回去——英文里同样是靠缩进
  ## 看从属关系，不能因为翻译就把层级抹平。
  pad <- if (grepl("^[　 ]", s)) "  " else ""
  s   <- sub("^[　 ]+", "", s)
  paste0(pad, .tr1_body(s))
}

## 不含缩进前缀的翻译本体（递归时也走这里，避免把缩进重复补一遍）
.tr1_body <- function(s) {
  if (!nzchar(s)) return(s)
  if (s %in% names(EN_TERMS)) return(unname(EN_TERMS[[s]]))
  ## "，n（%）"后缀：n 与百分比是数值，后缀是标签，只翻标签
  if (grepl("，n（%）$", s)) {
    base <- sub("，n（%）$", "", s)
    ## 先按带层记号的写法查（表 7 的词条就是这么写的），查不到再去掉记号重试
    b <- if (base %in% names(EN_TERMS)) unname(EN_TERMS[[base]])
         else .tr1_body(sub("^└ ", "", base))
    ## 层记号只跟着原文走：原文有 "└ " 才补，不然"男性，n（%）"也会被挂上一条横杆
    return(paste0(if (startsWith(base, "└ ")) "└ " else "", b, ", n (%)"))
  }
  ## 运行期生成的串
  for (r in EN_RULES) {
    m <- regmatches(s, regexec(r$p, s))[[1]]
    if (length(m)) return(do.call(r$f, as.list(m[-1])))
  }
  s
}

#' 中文 → 英文；zh 模式下恒等
tr <- function(x) {
  if (!is_en()) return(x)
  vapply(as.character(x), .tr1, character(1), USE.NAMES = FALSE)
}

#' 数值呈现串：翻译后再把全角括号换成半角（"0.93（0.90–0.96）" → "0.93 (0.90–0.96)"）
#'
#' 只对**已无中文字符**的串做标点转换：还有中文说明这个词条漏了，
#' 这时候硬换括号只会让漏译更难被发现。
tr_fmt <- function(x) {
  if (!is_en()) return(x)
  vapply(as.character(x), function(s) {
    s <- .tr1(s)
    if (!is.na(s) && !grepl("[一-鿿]", s)) {
      s <- gsub("（", " (", s, fixed = TRUE)
      s <- gsub("）", ")", s, fixed = TRUE)
      s <- gsub("，", ", ", s, fixed = TRUE)
    }
    s
  }, character(1), USE.NAMES = FALSE)
}

#' 整张表翻译：列名与字符型列都过一遍
tr_df <- function(d) {
  if (!is_en()) return(d)
  names(d) <- tr_fmt(names(d))
  for (j in seq_along(d)) if (is.character(d[[j]])) d[[j]] <- tr_fmt(d[[j]])
  d
}

## 表名（同时用作 CSV 文件名与 xlsx 页签名）
EN_TABLES <- c(
  "表1_基线特征"        = "Table1_Baseline",
  "表2_主分析"          = "Table2_Main",
  "表2b_扩展模型系数"   = "Table2b_FullCoefficients",
  "表3_剂量反应检验"    = "Table3_DoseResponseTests",
  "表3b_剂量反应节点"   = "Table3b_DoseResponsePoints",
  "表3c_每3mmol"        = "Table3c_Per3mmol",
  "表3d_判别与校准"     = "Table3d_Discrimination",
  "表4_预设亚组"        = "Table4_Subgroups",
  "表5_敏感性分析"      = "Table5_Sensitivity",
  "表5b_肌酐与尿量标准交叉" = "Table5b_CreatinineVsUrine",
  "表6_竞争风险_FineGray"  = "Table6_FineGray",
  "表6b_累积发生率"     = "Table6b_CumulativeIncidence",
  "表7_逐层留存"        = "Table7_Attrition",
  "表8_可判定与未纳入对比" = "Table8_EvaluableVsExcluded",
  "表9_监测强度分层"    = "Table9_MonitoringIntensity",
  "表10_复测概率模型"   = "Table10_RetestModel",
  "表17_绝对风险差与E值" = "Table17_AbsoluteRiskAndEValue",
  "表18_严重度梯度"     = "Table18_SeverityGradient",
  "表19_阴性对照"       = "Table19_NegativeControl"
)

tr_tabname <- function(nm) {
  if (!is_en()) return(nm)
  hit <- EN_TABLES[nm]
  ifelse(is.na(hit), nm, unname(hit))
}

#' 当前语言下的结果目录（中英文两套结果分开放，互不覆盖）
out_subdir <- function(zh_name, en_name) if (is_en()) en_name else zh_name
