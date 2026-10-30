## eICU-CRD 2.0 字段勘察（第二阶段提取适配层的依据）
##
## 目的只有一个：把 extract 适配层要用的**每一列**先核实存在、且取值符合预期。
## 勘察结论落 out/logs/eicu_recon.txt，并据此写 00/10/20/30/40_eicu_*.sql。
## 只读，不改数据。
source("extract.R")

LOG <- file.path(OUT_DIR, "logs", "eicu_recon.txt")
dir.create(dirname(LOG), showWarnings = FALSE, recursive = TRUE)
con <- dbConnect(duckdb())
on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)

out <- character(0)
say <- function(...) {
  s <- sprintf(...)
  out <<- c(out, s)
  cat(s, "\n", sep = "")
}
tbl <- function(t, quote = FALSE) sprintf(
  "read_csv_auto('%s/%s.csv.gz', header=true, sample_size=-1%s)",
  EICU_ROOT, t, if (quote) ", quote='\"'" else "")
q <- function(sql) dbGetQuery(con, sql)
try_q <- function(sql) tryCatch(q(sql), error = function(e) {
  say("  [错误] %s", conditionMessage(e)); NULL
})

say("eICU-CRD 2.0 字段勘察   %s", format(Sys.time()))
say("根目录：%s", EICU_ROOT)

## ---- 1. patient：年龄、停留时长、住院标识 -------------------------------------
say("\n## 1. patient")
say("  行数：%s", q(sprintf("SELECT COUNT(*) FROM %s", tbl("patient")))[[1]])
say("  age 的取值形态（前 20 个不同值）：%s",
    paste(q(sprintf("SELECT DISTINCT age FROM %s LIMIT 20", tbl("patient")))[[1]],
          collapse = " | "))
say("  含非数字年龄的行数：%s",
    q(sprintf("SELECT COUNT(*) FROM %s WHERE TRY_CAST(age AS INTEGER) IS NULL",
              tbl("patient")))[[1]])
say("  unitdischargeoffset 的分位：%s",
    paste(q(sprintf("SELECT MIN(unitdischargeoffset), MEDIAN(unitdischargeoffset),
                            MAX(unitdischargeoffset) FROM %s", tbl("patient")))[1, ],
          collapse = " / "))
say("  hospitaldischargeoffset 存在：%s",
    "hospitaldischargeoffset" %in% names(q(sprintf("SELECT * FROM %s LIMIT 0", tbl("patient")))))

## ---- 2. lab：血钠与肌酐 ---------------------------------------------------------
say("\n## 2. lab")
say("  行数：%s", q(sprintf("SELECT COUNT(*) FROM %s", tbl("lab")))[[1]])
say("  labname 中匹配 sodium 的取值：")
print(try_q(sprintf(
  "SELECT labname, COUNT(*) AS n, MIN(labresultoffset) AS off_min,
          MAX(labresultoffset) AS off_max
   FROM %s WHERE LOWER(labname) LIKE '%%sodium%%'
   GROUP BY labname ORDER BY n DESC", tbl("lab"))))
say("  labname 中匹配 creat 的取值：")
print(try_q(sprintf(
  "SELECT labname, COUNT(*) AS n, MIN(labresultoffset) AS off_min,
          MAX(labresultoffset) AS off_max
   FROM %s WHERE LOWER(labname) LIKE '%%creat%%'
   GROUP BY labname ORDER BY n DESC", tbl("lab"))))
say("  labresult 非数值行数（血钠）：%s",
    q(sprintf("SELECT COUNT(*) FROM %s WHERE LOWER(labname) LIKE '%%sodium%%'
               AND TRY_CAST(labresult AS DOUBLE) IS NULL", tbl("lab")))[[1]])
say("  labresult 非数值行数（肌酐）：%s",
    q(sprintf("SELECT COUNT(*) FROM %s WHERE LOWER(labname) LIKE '%%creat%%'
               AND TRY_CAST(labresult AS DOUBLE) IS NULL", tbl("lab")))[[1]])
say("  glucose 相关 labname：%s",
    paste(q(sprintf("SELECT DISTINCT labname FROM %s WHERE LOWER(labname) LIKE '%%glucose%%'",
                    tbl("lab")))[[1]], collapse = " | "))

## ---- 3. 尿量 ------------------------------------------------------------------
say("\n## 3. intakeOutput")
say("  行数：%s", q(sprintf("SELECT COUNT(*) FROM %s", tbl("intakeOutput")))[[1]])
say("  与尿相关的 celllabel（前 25）：")
print(try_q(sprintf(
  "SELECT celllabel, cellpath, COUNT(*) AS n
   FROM %s WHERE LOWER(celllabel) LIKE '%%urine%%' OR LOWER(cellpath) LIKE '%%urine%%'
   GROUP BY celllabel, cellpath ORDER BY n DESC LIMIT 25", tbl("intakeOutput"))))
say("  intakeoutputoffset 范围：%s",
    paste(q(sprintf("SELECT MIN(intakeoutputoffset), MAX(intakeoutputoffset) FROM %s",
                    tbl("intakeOutput")))[1, ], collapse = " / "))

## ---- 4. 血压 ------------------------------------------------------------------
say("\n## 4. 血压来源")
say("  vitalAperiodic 的 vitalparam 取值：")
print(try_q(sprintf(
  "SELECT vitalparam, COUNT(*) AS n FROM %s GROUP BY vitalparam ORDER BY n DESC LIMIT 15",
  tbl("vitalAperiodic"))))
say("  nurseCharting 中血压相关 celltypevallabel：")
print(try_q(sprintf(
  "SELECT nursingchartcelltypevallabel, COUNT(*) AS n FROM %s
   WHERE LOWER(nursingchartcelltypevallabel) LIKE '%%blood pressure%%'
      OR LOWER(nursingchartcelltypevallabel) LIKE '%%arterial%%'
   GROUP BY 1 ORDER BY n DESC LIMIT 12", tbl("nurseCharting", quote = TRUE))))

## ---- 5. RRT 与 HF 诊断 ---------------------------------------------------------
say("\n## 5. treatment 中的 RRT")
print(try_q(sprintf(
  "SELECT COUNT(DISTINCT patientunitstayid) AS units, COUNT(*) AS n FROM %s
   WHERE LOWER(treatmentstring) LIKE '%%dialysis%%'
      OR LOWER(treatmentstring) LIKE '%%crrt%%'
      OR LOWER(treatmentstring) LIKE '%%renal replacement%%'",
  tbl("treatment"))))
say("  透析相关 treatmentstring 取值（前 15）：")
print(try_q(sprintf(
  "SELECT treatmentstring, COUNT(*) AS n FROM %s
   WHERE LOWER(treatmentstring) LIKE '%%dialysis%%' OR LOWER(treatmentstring) LIKE '%%crrt%%'
   GROUP BY 1 ORDER BY n DESC LIMIT 15", tbl("treatment"))))

say("\n## 6. diagnosis 中的心衰码")
print(try_q(sprintf(
  "SELECT COUNT(DISTINCT patientunitstayid) AS units FROM %s
   WHERE icd9code LIKE '428%%' OR icd9code LIKE '%%428%%'", tbl("diagnosis"))))
say("  含 I50 的行：%s",
    q(sprintf("SELECT COUNT(*) FROM %s WHERE icd9code LIKE '%%I50%%'", tbl("diagnosis")))[[1]])

## ---- 7. apachePatientResult：确认含暴露/结局的列 -------------------------------
say("\n## 7. apachePatientResult / apacheApsVar（**确认这些列必须排除在模型外**）")
print(try_q(sprintf("SELECT apacheversion, COUNT(*) FROM %s GROUP BY 1",
                    tbl("apachePatientResult"))))

## ---- 8. infusionDrug / medication ---------------------------------------------
say("\n## 8. 血管活性药与利尿剂")
print(try_q(sprintf(
  "SELECT drugname, COUNT(*) AS n FROM %s
   WHERE LOWER(drugname) LIKE '%%norepinephrine%%' OR LOWER(drugname) LIKE '%%dopamine%%'
      OR LOWER(drugname) LIKE '%%dobutamine%%' OR LOWER(drugname) LIKE '%%epinephrine%%'
      OR LOWER(drugname) LIKE '%%vasopressin%%' OR LOWER(drugname) LIKE '%%milrinone%%'
   GROUP BY 1 ORDER BY n DESC LIMIT 20", tbl("infusionDrug"))))
print(try_q(sprintf(
  "SELECT drugname, COUNT(*) AS n FROM %s
   WHERE LOWER(drugname) LIKE '%%furosemide%%' OR LOWER(drugname) LIKE '%%bumetanide%%'
      OR LOWER(drugname) LIKE '%%torsemide%%'
   GROUP BY 1 ORDER BY n DESC LIMIT 15", tbl("medication"))))

writeLines(out, LOG, useBytes = TRUE)
cat(sprintf("\n已写入 %s\n", LOG))
