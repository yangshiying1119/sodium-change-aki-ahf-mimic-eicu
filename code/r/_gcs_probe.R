## 只问一件事：apacheApsVar 的 GCS 三项在 eICU 心衰队列内的覆盖率。
## 这直接决定 model_cases() 会不会把**基础**模型的样本量也拖下去——
## 因为 MODEL_VARS = BASIC_X + EXTEND_X，是整体取完整病例的。
source("extract.R"); suppressMessages(library(DBI))
con <- dbConnect(duckdb()); on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
dbExecute(con, sprintf("CREATE OR REPLACE VIEW apache_gcs AS SELECT patientunitstayid, eyes, motor, verbal FROM read_csv_auto('%s/apacheApsVar.csv.gz', header=true, escape='\"')", EICU_ROOT))
dbExecute(con, sprintf("CREATE OR REPLACE VIEW patient AS SELECT * FROM read_csv_auto('%s/patient.csv.gz', header=true, escape='\"')", EICU_ROOT))
cat("## A. 全库 apacheApsVar 的 GCS 三项分布（-1 = 缺失哨兵）\n")
print(dbGetQuery(con, "
  SELECT COUNT(*) AS 行数,
         SUM(CASE WHEN eyes > 0 AND motor > 0 AND verbal > 0 THEN 1 ELSE 0 END) AS 三项齐全,
         SUM(CASE WHEN eyes = -1 OR motor = -1 OR verbal = -1 THEN 1 ELSE 0 END) AS 含缺失哨兵,
         MIN(eyes) AS eyes_lo, MAX(eyes) AS eyes_hi,
         MIN(motor) AS motor_lo, MAX(motor) AS motor_hi,
         MIN(verbal) AS verbal_lo, MAX(verbal) AS verbal_hi
  FROM apache_gcs"))
cat("\n## B. 每个 patientunitstayid 的行数（是否一单元多行）\n")
print(dbGetQuery(con, "SELECT 行数, COUNT(*) AS 单元数 FROM (SELECT patientunitstayid, COUNT(*) AS 行数 FROM apache_gcs GROUP BY 1) GROUP BY 1 ORDER BY 1 LIMIT 10"))
cat("\n## C. 成人且停留>24h 的单元里，GCS 三项齐全的比例\n")
print(dbGetQuery(con, "
  WITH b AS (SELECT patientunitstayid FROM patient
             WHERE unitdischargeoffset IS NOT NULL AND unitdischargeoffset > 1440
               AND CASE WHEN TRY_CAST(age AS INTEGER) IS NOT NULL THEN TRY_CAST(age AS INTEGER)
                        WHEN TRIM(COALESCE(age,'')) LIKE '>%' THEN 90 ELSE NULL END >= 18)
  SELECT COUNT(*) AS 单元数,
         SUM(CASE WHEN g.patientunitstayid IS NOT NULL THEN 1 ELSE 0 END) AS 有GCS行,
         SUM(CASE WHEN g.eyes>0 AND g.motor>0 AND g.verbal>0 THEN 1 ELSE 0 END) AS 三项齐全
  FROM b LEFT JOIN apache_gcs g USING (patientunitstayid)"))
cat("\n已写入（探针结束）\n")
