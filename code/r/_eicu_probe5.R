## 定向补探：只问 40_covariates_eicu.sql 里**还没有把握**的几处谓词。
##
## 与 _eicu_recon4.R 的区别：recon4 是"把相关表都看一遍"，本脚本只回答几个
## 二选一的问题，且**全部避开 nurseCharting**（那是全库最大的表，recon4 正卡在
## 它上面）。每个查询都包在 try_q 里，一处失败不拖垮其余。
source("extract.R")
suppressMessages(library(DBI))

con <- dbConnect(duckdb())
on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)

say  <- function(...) cat(sprintf(...), "\n", sep = "")
try_q <- function(label, sql) {
  cat("\n## ", label, "\n", sep = "")
  r <- tryCatch(dbGetQuery(con, sql), error = function(e) {
    cat("  [查询失败] ", conditionMessage(e), "\n", sep = ""); NULL })
  if (!is.null(r)) print(r, row.names = FALSE)
  invisible(r)
}

## 建视图（只读本探针要用的几张表，不建全部，省去护士记录表的打开开销）
for (f in c("intakeOutput", "respiratoryCare", "treatment", "infusionDrug",
            "medication", "vitalAperiodic", "lab", "patient"))
  dbExecute(con, sprintf(
    "CREATE OR REPLACE VIEW %s AS SELECT * FROM read_csv_auto('%s/%s.csv.gz', header = true, escape = '\"')",
    f, EICU_ROOT, f))

## ---- 1. intakeOutput：入量侧的 cellpath 到底是什么前缀 ----
## 这一条最要紧：40_covariates_eicu.sql 的 cov_fluid_in 现在写死
## 'flowsheet|Flowsheet Cell Labels|I&O|Intake%'。若前缀不对，fluid_in_24h_ml
## 会全为 0——而 0 是**合法取值**，check_covariates() 的"不同取值 ≤ 1"判据能拦住，
## 但更早在这里问清代价更小。
try_q("1a intakeOutput：cellpath 的完整取值分布（前 30）", "
  SELECT cellpath, COUNT(*) AS n,
         SUM(CASE WHEN cellvaluenumeric IS NOT NULL THEN 1 ELSE 0 END) AS 有数值
  FROM intakeOutput GROUP BY cellpath ORDER BY n DESC LIMIT 30")

try_q("1b intakeOutput：入量相关 cellpath 的细目（含 Intake 字样）", "
  SELECT cellpath, celllabel, COUNT(*) AS n
  FROM intakeOutput WHERE LOWER(cellpath) LIKE '%intake%'
  GROUP BY cellpath, celllabel ORDER BY n DESC LIMIT 40")

try_q("1c intakeOutput：出量尿量条目抽样（核对 21_urine_eicu.sql 的判据）", "
  SELECT cellpath, celllabel, COUNT(*) AS n
  FROM intakeOutput
  WHERE cellpath LIKE 'flowsheet|Flowsheet Cell Labels|I&O|Output (ml)|%'
    AND LOWER(TRIM(celllabel)) LIKE '%urine%'
  GROUP BY cellpath, celllabel ORDER BY n DESC LIMIT 25")

try_q("1d intakeOutput：出量区里与 urine 相关但会被排除规则滤掉的条目", "
  SELECT celllabel, COUNT(*) AS n
  FROM intakeOutput
  WHERE cellpath LIKE 'flowsheet|Flowsheet Cell Labels|I&O|Output (ml)|%'
    AND (LOWER(TRIM(celllabel)) LIKE '%count%' OR LOWER(TRIM(celllabel)) LIKE '%occurrence%'
      OR LOWER(TRIM(celllabel)) LIKE '%incontinen%' OR LOWER(TRIM(celllabel)) LIKE '%stool%')
  GROUP BY celllabel ORDER BY n DESC LIMIT 15")

## ---- 2. respiratoryCare：airwaytype 的真实取值 ----
try_q("2 respiratoryCare：airwaytype 取值分布", "
  SELECT airwaytype, COUNT(*) AS n, COUNT(DISTINCT patientunitstayid) AS units
  FROM respiratoryCare GROUP BY airwaytype ORDER BY n DESC LIMIT 25")

## ---- 3. treatment：与通气相关的 treatmentstring ----
try_q("3 treatment：含 ventilat / cpap / bipap 的 treatmentstring", "
  SELECT treatmentstring, COUNT(DISTINCT patientunitstayid) AS units
  FROM treatment
  WHERE LOWER(treatmentstring) LIKE '%ventilat%'
     OR LOWER(treatmentstring) LIKE '%cpap%'
     OR LOWER(treatmentstring) LIKE '%bipap%'
  GROUP BY treatmentstring ORDER BY units DESC LIMIT 30")

## ---- 4. infusionDrug：血管活性药与高渗盐水的完整 drugname（带单位后缀） ----
try_q("4a infusionDrug：血管活性药相关 drugname 全量取值", "
  SELECT drugname, COUNT(*) AS n FROM infusionDrug
  WHERE LOWER(drugname) LIKE '%norepinephrine%' OR LOWER(drugname) LIKE '%epinephrine%'
     OR LOWER(drugname) LIKE '%dopamine%'       OR LOWER(drugname) LIKE '%vasopressin%'
     OR LOWER(drugname) LIKE '%phenylephrine%'  OR LOWER(drugname) LIKE '%neosynephrine%'
     OR LOWER(drugname) LIKE '%milrinone%'      OR LOWER(drugname) LIKE '%dobutamine%'
  GROUP BY drugname ORDER BY n DESC LIMIT 40")

try_q("4b infusionDrug：高渗盐水相关 drugname", "
  SELECT drugname, COUNT(*) AS n FROM infusionDrug
  WHERE LOWER(drugname) LIKE '%sodium chloride%' OR LOWER(drugname) LIKE '%nacl%'
     OR LOWER(drugname) LIKE '%hypertonic%' OR LOWER(drugname) LIKE '%3%%'
  GROUP BY drugname ORDER BY n DESC LIMIT 30")

## ---- 5. medication：利尿剂与血管活性药的药名形态 ----
try_q("5 medication：利尿剂相关 drugname 与给药途径", "
  SELECT drugname, routeadmin, COUNT(*) AS n FROM medication
  WHERE LOWER(drugname) LIKE '%furosemide%' OR LOWER(drugname) LIKE '%lasix%'
     OR LOWER(drugname) LIKE '%bumetanide%' OR LOWER(drugname) LIKE '%bumex%'
     OR LOWER(drugname) LIKE '%torsemide%'  OR LOWER(drugname) LIKE '%torasemide%'
  GROUP BY drugname, routeadmin ORDER BY n DESC LIMIT 30")

## ---- 6. lab：血糖 / 血红蛋白的 labmeasurenamesystem（核对单位判据） ----
try_q("6 lab：glucose 与 Hgb 的单位体系", "
  SELECT labname, labmeasurenamesystem, COUNT(*) AS n,
         MIN(labresult) AS lo, MAX(labresult) AS hi
  FROM lab WHERE labname IN ('glucose','bedside glucose','Hgb','hemoglobin')
  GROUP BY labname, labmeasurenamesystem ORDER BY n DESC LIMIT 20")

## ---- 7. patient：admissionweight 覆盖（决定是否需要 nurseCharting 兜底） ----
try_q("7 patient：admissionweight 覆盖与范围", "
  SELECT COUNT(*) AS 全部单元,
         SUM(CASE WHEN admissionweight IS NOT NULL AND admissionweight > 0 THEN 1 ELSE 0 END) AS 有体重,
         MIN(admissionweight) AS lo, MAX(admissionweight) AS hi
  FROM patient")

cat("\n已写入（探针结束）\n")
