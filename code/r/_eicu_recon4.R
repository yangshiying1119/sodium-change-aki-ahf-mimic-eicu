## eICU 第三轮勘察：只问写 40_covariates_eicu.sql 还缺的字段。
## 上一轮（_eicu_recon3.R）用了不存在的 nursechartcelltypelabel 列，且把大量
## 全表聚合串在一起跑太久；本轮按"一张表一次扫描"组织，只留必要的问题。
source("extract.R")

LOG <- file.path(OUT_DIR, "logs", "eicu_recon4.txt")
con <- dbConnect(duckdb())
on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)

out <- character(0)
say <- function(...) { s <- sprintf(...); out <<- c(out, s); cat(s, "\n", sep = "") }
T_ <- function(t, quote = FALSE) sprintf(
  "read_csv_auto('%s/%s.csv.gz', header=true, sample_size=-1%s)",
  EICU_ROOT, t, if (quote) ", quote='\"'" else "")
q <- function(sql) dbGetQuery(con, sql)
try_q <- function(sql) tryCatch(q(sql), error = function(e) {
  say("  [错误] %s", conditionMessage(e)); NULL })
show <- function(x, n = 30) if (!is.null(x)) print(utils::head(as.data.frame(x), n), row.names = FALSE)

say("eICU 第三轮勘察   %s", format(Sys.time()))

## ---- 1. patient 的分类字段（小表） --------------------------------------------
say("\n## 1. patient：unittype / unitadmitsource / 出院状态 / unitvisitnumber")
show(try_q(sprintf("SELECT unittype, COUNT(*) AS n FROM %s GROUP BY 1 ORDER BY n DESC",
                   T_("patient"))))
show(try_q(sprintf("SELECT unitadmitsource, COUNT(*) AS n FROM %s GROUP BY 1 ORDER BY n DESC",
                   T_("patient"))), 12)
show(try_q(sprintf(
  "SELECT unitdischargestatus, hospitaldischargestatus, COUNT(*) AS n
   FROM %s GROUP BY 1,2 ORDER BY n DESC", T_("patient"))))
show(try_q(sprintf(
  "SELECT unitvisitnumber, COUNT(*) AS n FROM %s GROUP BY 1 ORDER BY n DESC",
  T_("patient"))), 10)
say("  hospitaldischargeyear：")
show(try_q(sprintf("SELECT hospitaldischargeyear, COUNT(*) AS n FROM %s GROUP BY 1 ORDER BY 1",
                   T_("patient"))))

## ---- 2. admissionDx 的问诊名（小表） ------------------------------------------
say("\n## 2. admissionDx：含 CHF / heart failure 的问诊名")
show(try_q(sprintf(
  "SELECT admitdxname, COUNT(DISTINCT patientunitstayid) AS units, COUNT(*) AS n
   FROM %s
   WHERE LOWER(admitdxname) LIKE '%%chf%%' OR LOWER(admitdxname) LIKE '%%heart failure%%'
   GROUP BY 1 ORDER BY n DESC", T_("admissionDx"))), 40)

## ---- 3. lab：单位判据 + 血红蛋白名 --------------------------------------------
say("\n## 3. lab：sodium / creatinine 的 labmeasurenamesystem")
show(try_q(sprintf(
  "SELECT labname, labmeasurenamesystem, COUNT(*) AS n FROM %s
   WHERE labname IN ('sodium','creatinine') GROUP BY 1,2 ORDER BY 1, n DESC", T_("lab"))), 20)
say("\n  labname 中与血红蛋白/血色素相关的取值：")
show(try_q(sprintf(
  "SELECT labname, labmeasurenamesystem, COUNT(*) AS n, MIN(labresult) AS lo, MAX(labresult) AS hi
   FROM %s WHERE LOWER(labname) LIKE '%%hgb%%' OR LOWER(labname) LIKE '%%hemoglob%%'
      OR LOWER(labname) LIKE '%%haemoglob%%' OR LOWER(labname) LIKE '%%hct%%'
      OR LOWER(labname) LIKE '%%hematocrit%%'
   GROUP BY 1,2 ORDER BY n DESC", T_("lab"))), 20)
say("\n  labname 中与血糖相关的取值与单位：")
show(try_q(sprintf(
  "SELECT labname, labmeasurenamesystem, COUNT(*) AS n FROM %s
   WHERE LOWER(labname) LIKE '%%glucose%%' GROUP BY 1,2 ORDER BY n DESC", T_("lab"))), 20)

## ---- 4. diagnosis：diagnosispriority 分布（决定 A 口径能否落地） ---------------
say("\n## 4. diagnosis：diagnosispriority 分布与心衰码所在位次")
show(try_q(sprintf(
  "SELECT diagnosispriority, COUNT(*) AS n FROM %s GROUP BY 1 ORDER BY 1",
  T_("diagnosis"))), 15)
show(try_q(sprintf(
  "SELECT diagnosispriority, COUNT(DISTINCT patientunitstayid) AS units
   FROM %s WHERE icd9code LIKE '428%%' OR icd9code LIKE '%%I50%%' GROUP BY 1 ORDER BY 1",
  T_("diagnosis"))), 15)

## ---- 5. nurseCharting：体重 / 体温 / GCS 的标签（大表，一次扫描） -------------
say("\n## 5. nurseCharting：体重 / 体温 / GCS 标签")
show(try_q(sprintf(
  "SELECT nursingchartcelltypecat, nursingchartcelltypevallabel,
          nursingchartcelltypevalname, COUNT(*) AS n
   FROM %s
   WHERE LOWER(nursingchartcelltypevallabel) LIKE '%%weight%%'
      OR LOWER(nursingchartcelltypevalname) LIKE '%%weight%%'
      OR LOWER(nursingchartcelltypevalname) LIKE '%%temp%%'
      OR LOWER(nursingchartcelltypevallabel) LIKE '%%temp%%'
      OR LOWER(nursingchartcelltypevallabel) LIKE '%%gcs%%'
      OR LOWER(nursingchartcelltypevalname) LIKE '%%gcs%%'
      OR LOWER(nursingchartcelltypevallabel) LIKE '%%glasgow%%'
   GROUP BY 1,2,3 ORDER BY n DESC", T_("nurseCharting", quote = TRUE))), 30)

## ---- 6. intakeOutput：入量侧 cellpath（大表，一次扫描） -----------------------
say("\n## 6. intakeOutput：入量侧 cellpath 顶层分类")
show(try_q(sprintf(
  "SELECT SPLIT_PART(cellpath, '|', 4) AS lvl4, COUNT(*) AS n FROM %s
   GROUP BY 1 ORDER BY n DESC", T_("intakeOutput"))), 20)
say("\n  cellpath 含 Intake 的前 15 个完整路径：")
show(try_q(sprintf(
  "SELECT cellpath, COUNT(*) AS n FROM %s WHERE cellpath LIKE '%%Intake%%'
   GROUP BY 1 ORDER BY n DESC LIMIT 15", T_("intakeOutput"))), 15)

## ---- 7. 治疗与呼吸（各自一次扫描） --------------------------------------------
say("\n## 7. treatment：机械通气相关条目")
show(try_q(sprintf(
  "SELECT treatmentstring, COUNT(DISTINCT patientunitstayid) AS units, COUNT(*) AS n
   FROM %s WHERE LOWER(treatmentstring) LIKE '%%ventilation%%'
      OR LOWER(treatmentstring) LIKE '%%intubat%%' OR LOWER(treatmentstring) LIKE '%%cpap%%'
   GROUP BY 1 ORDER BY n DESC", T_("treatment"))), 25)
say("\n  透析类 treatmentoffset 范围：")
show(try_q(sprintf(
  "SELECT MIN(treatmentoffset) AS lo, MAX(treatmentoffset) AS hi, COUNT(*) AS n
   FROM %s WHERE LOWER(treatmentstring) LIKE '%%dialysis%%'", T_("treatment"))))

say("\n## 8. respiratoryCare：airwaytype 分布")
show(try_q(sprintf(
  "SELECT airwaytype, COUNT(DISTINCT patientunitstayid) AS units, COUNT(*) AS n
   FROM %s GROUP BY 1 ORDER BY n DESC", T_("respiratoryCare", quote = TRUE))), 20)

## ---- 9. vitalAperiodic 覆盖 ---------------------------------------------------
say("\n## 9. vitalAperiodic 血压覆盖")
show(try_q(sprintf(
  "SELECT COUNT(DISTINCT patientunitstayid) AS units,
          SUM(CASE WHEN noninvasivesystolic IS NOT NULL THEN 1 ELSE 0 END) AS 非空收缩,
          MIN(observationoffset) AS lo, MAX(observationoffset) AS hi
   FROM %s", T_("vitalAperiodic"))))

## ---- 10. 血管活性药 / 利尿剂 / 高渗盐水 ---------------------------------------
say("\n## 10. infusionDrug：苯肾上腺素 / 加压素 / 高渗盐水")
show(try_q(sprintf(
  "SELECT drugname, COUNT(DISTINCT patientunitstayid) AS units, COUNT(*) AS n FROM %s
   WHERE LOWER(drugname) LIKE '%%phenylephrine%%' OR LOWER(drugname) LIKE '%%neosynephrine%%'
      OR LOWER(drugname) LIKE '%%vasopressin%%' OR LOWER(drugname) LIKE '%%sodium chloride%%'
      OR LOWER(drugname) LIKE '%%nacl%%'
   GROUP BY 1 ORDER BY n DESC", T_("infusionDrug"))), 30)
say("\n  medication：呋塞米/布美他尼的 routeadmin 分布")
show(try_q(sprintf(
  "SELECT routeadmin, COUNT(*) AS n FROM %s
   WHERE LOWER(drugname) LIKE '%%furosemide%%' OR LOWER(drugname) LIKE '%%bumetanide%%'
   GROUP BY 1 ORDER BY n DESC", T_("medication"))), 20)
say("\n  medication：dosage 样例")
show(try_q(sprintf(
  "SELECT drugname, dosage, frequency, routeadmin, COUNT(*) AS n FROM %s
   WHERE LOWER(drugname) LIKE '%%furosemide%%' GROUP BY 1,2,3,4 ORDER BY n DESC",
  T_("medication"))), 15)

writeLines(out, LOG, useBytes = TRUE)
say("\n已写入 %s", LOG)
