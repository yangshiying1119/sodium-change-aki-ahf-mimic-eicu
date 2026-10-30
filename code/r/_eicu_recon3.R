## eICU 侧**第二轮**字段勘察：只核实 30_cohort_eicu.sql 与 40_covariates_eicu.sql
## 要写的每一列。第一轮（_eicu_recon2.R）确认了 labname / treatmentstring / 尿量
## 标签，但没打印 patient 的列名、admissionDx 的问诊名、lab 的单位判据。
## 结论落 out/logs/eicu_recon3.txt。
source("extract.R")

LOG <- file.path(OUT_DIR, "logs", "eicu_recon3.txt")
con <- dbConnect(duckdb())
on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)

out <- character(0)
say <- function(...) {
  s <- sprintf(...)
  out <<- c(out, s); cat(s, "\n", sep = "")
}
T_ <- function(t, quote = FALSE) sprintf(
  "read_csv_auto('%s/%s.csv.gz', header=true, sample_size=-1%s)",
  EICU_ROOT, t, if (quote) ", quote='\"'" else "")
q <- function(sql) dbGetQuery(con, sql)
try_q <- function(sql) tryCatch(q(sql), error = function(e) {
  say("  [错误] %s", conditionMessage(e)); NULL
})
show <- function(x, n = 30) {
  if (is.null(x)) return(invisible(NULL))
  print(utils::head(as.data.frame(x), n), row.names = FALSE)
}

say("eICU 第二轮勘察   %s", format(Sys.time()))

## ---- 1. patient 的列名与关键字段取值 ------------------------------------------
say("\n## 1. patient 列名")
print(names(q(sprintf("SELECT * FROM %s LIMIT 0", T_("patient")))))
say("\n  unitdischargestatus / hospitaldischargestatus 取值：")
show(try_q(sprintf(
  "SELECT unitdischargestatus, hospitaldischargestatus, COUNT(*) AS n
   FROM %s GROUP BY 1, 2 ORDER BY n DESC", T_("patient"))))
say("\n  unitdischargeoffset 分层（分钟）：")
show(try_q(sprintf(
  "SELECT
     COUNT(*) AS 全部,
     SUM(CASE WHEN unitdischargeoffset IS NULL OR unitdischargeoffset <= 0 THEN 1 ELSE 0 END) AS 非正,
     SUM(CASE WHEN unitdischargeoffset > 1440 THEN 1 ELSE 0 END) AS 超24h,
     SUM(CASE WHEN unitdischargeoffset > 4320 THEN 1 ELSE 0 END) AS 超72h
   FROM %s", T_("patient"))))
say("\n  unittype 取值：")
show(try_q(sprintf("SELECT unittype, COUNT(*) AS n FROM %s GROUP BY 1 ORDER BY n DESC",
                   T_("patient"))))
say("\n  admissionweight 缺失：")
show(try_q(sprintf(
  "SELECT COUNT(*) AS 全部, SUM(CASE WHEN admissionweight IS NULL THEN 1 ELSE 0 END) AS 缺失,
          MIN(admissionweight) AS 最小, MAX(admissionweight) AS 最大
   FROM %s", T_("patient"))))

## ---- 2. admissionDx 的问诊名 --------------------------------------------------
say("\n## 2. admissionDx：含 CHF / heart failure 的问诊名")
show(try_q(sprintf(
  "SELECT admitdxname, COUNT(DISTINCT patientunitstayid) AS units, COUNT(*) AS n
   FROM %s
   WHERE LOWER(admitdxname) LIKE '%%chf%%' OR LOWER(admitdxname) LIKE '%%heart failure%%'
      OR LOWER(admitdxname) LIKE '%%cardiac%%'
   GROUP BY 1 ORDER BY n DESC", T_("admissionDx"))), 40)

## ---- 3. lab 的单位判据 --------------------------------------------------------
say("\n## 3. lab：sodium / creatinine 的 labmeasurenamesystem 取值")
show(try_q(sprintf(
  "SELECT labname, labmeasurenamesystem, COUNT(*) AS n
   FROM %s WHERE labname IN ('sodium','creatinine')
   GROUP BY 1, 2 ORDER BY labname, n DESC", T_("lab"))), 20)

## ---- 4. treatment 的列名与 RRT 的时点字段 ------------------------------------
say("\n## 4. treatment 列名")
print(names(q(sprintf("SELECT * FROM %s LIMIT 0", T_("treatment")))))
say("\n  treatmentoffset 是否可用（透析类）：")
show(try_q(sprintf(
  "SELECT MIN(treatmentoffset) AS 最小, MAX(treatmentoffset) AS 最大, COUNT(*) AS n
   FROM %s WHERE LOWER(treatmentstring) LIKE '%%dialysis%%'", T_("treatment"))))

## ---- 5. intakeOutput / nurseCharting 列名 ------------------------------------
say("\n## 5. intakeOutput 列名")
print(names(q(sprintf("SELECT * FROM %s LIMIT 0", T_("intakeOutput")))))
say("\n## 6. nurseCharting 列名")
print(names(q(sprintf("SELECT * FROM %s LIMIT 0", T_("nurseCharting", quote = TRUE)))))
say("\n  nurseCharting 中体重与体温相关标签：")
show(try_q(sprintf(
  "SELECT nursingchartcelltypecat, nursingchartcelltypevallabel, nursingchartcelltypelabel,
          COUNT(*) AS n
   FROM %s
   WHERE LOWER(nursingchartcelltypelabel) LIKE '%%weight%%'
      OR LOWER(nursingchartcelltypevallabel) LIKE '%%weight%%'
      OR LOWER(nursingchartcelltypelabel) LIKE '%%temp%%'
   GROUP BY 1, 2, 3 ORDER BY n DESC", T_("nurseCharting", quote = TRUE))), 20)

## ---- 7. 血压：vitalAperiodic 具名列 + nurseCharting 覆盖 ----------------------
say("\n## 7. 血压覆盖")
show(try_q(sprintf(
  "SELECT COUNT(DISTINCT patientunitstayid) AS 单元数,
          SUM(CASE WHEN noninvasivesystolic IS NOT NULL THEN 1 ELSE 0 END) AS 非空收缩
   FROM %s", T_("vitalAperiodic"))))
show(try_q(sprintf(
  "SELECT nursingchartcelltypevallabel, COUNT(DISTINCT patientunitstayid) AS 单元数, COUNT(*) AS n
   FROM %s
   WHERE LOWER(nursingchartcelltypevallabel) LIKE '%%bp%%'
      OR LOWER(nursingchartcelltypevallabel) LIKE '%%blood pressure%%'
   GROUP BY 1 ORDER BY n DESC", T_("nurseCharting", quote = TRUE))), 20)

## ---- 8. pastHistory 结构（合并症） --------------------------------------------
say("\n## 8. pastHistory 列名与路径层级样例")
print(names(q(sprintf("SELECT * FROM %s LIMIT 0", T_("pastHistory", quote = TRUE)))))
show(try_q(sprintf(
  "SELECT pasthistorypath, COUNT(*) AS n FROM %s
   WHERE LOWER(pasthistorypath) LIKE '%%renal%%' OR LOWER(pasthistorypath) LIKE '%%diabet%%'
      OR LOWER(pasthistorypath) LIKE '%%hypertension%%' OR LOWER(pasthistorypath) LIKE '%%coronary%%'
   GROUP BY 1 ORDER BY n DESC", T_("pastHistory", quote = TRUE))), 30)

writeLines(out, LOG, useBytes = TRUE)
cat(sprintf("\n已写入 %s\n", LOG))
