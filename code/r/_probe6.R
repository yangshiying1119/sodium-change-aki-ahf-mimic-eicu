source("extract.R"); suppressMessages(library(DBI))
con <- dbConnect(duckdb()); on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
for (f in c("nurseCharting", "medication", "intakeOutput"))
  dbExecute(con, sprintf("CREATE OR REPLACE VIEW %s AS SELECT * FROM read_csv_auto('%s/%s.csv.gz', header=true, escape='\"')", f, EICU_ROOT, f))
q <- function(lab, sql) { cat("\n## ", lab, "\n", sep=""); print(tryCatch(dbGetQuery(con, sql), error=function(e) cat("[失败] ", conditionMessage(e), "\n", sep=""))) }

q("A nurseCharting GCS Total 的取值范围（是否 3-15 的合成分）", "
  SELECT nursingchartcelltypevalname AS 条目,
         COUNT(*) AS n,
         SUM(CASE WHEN TRY_CAST(nursingchartvalue AS INTEGER) IS NOT NULL THEN 1 ELSE 0 END) AS 可转整数,
         MIN(TRY_CAST(nursingchartvalue AS INTEGER)) AS lo,
         MAX(TRY_CAST(nursingchartvalue AS INTEGER)) AS hi
  FROM nurseCharting
  WHERE nursingchartcelltypevalname IN ('GCS Total','Eyes','Motor','Verbal')
    AND nursingchartcelltypecat = 'Scores'
  GROUP BY 1 ORDER BY n DESC")

q("B nurseCharting GCS 条目非整数取值的抽样（确认没有文本混入）", "
  SELECT nursingchartvalue, COUNT(*) AS n FROM nurseCharting
  WHERE nursingchartcelltypevalname IN ('GCS Total','Eyes','Motor','Verbal')
    AND nursingchartcelltypecat = 'Scores'
    AND TRY_CAST(nursingchartvalue AS INTEGER) IS NULL
  GROUP BY 1 ORDER BY n DESC LIMIT 10")

q("C medication.dosage 的取值形态（决定呋塞米当量能否解析）", "
  SELECT dosage, dosageunit, COUNT(*) AS n FROM medication
  WHERE LOWER(drugname) LIKE '%furosemide%' OR LOWER(drugname) LIKE '%lasix%'
     OR LOWER(drugname) LIKE '%bumetanide%' OR LOWER(drugname) LIKE '%bumex%'
  GROUP BY 1,2 ORDER BY n DESC LIMIT 25")

q("D intakeOutput 体重条目的单位与取值", "
  SELECT cellpath, celllabel, COUNT(*) AS n,
         MIN(TRY_CAST(cellvaluenumeric AS DOUBLE)) AS lo,
         MAX(TRY_CAST(cellvaluenumeric AS DOUBLE)) AS hi
  FROM intakeOutput WHERE cellpath LIKE '%|Weight|%'
  GROUP BY 1,2 ORDER BY n DESC LIMIT 10")

q("E 高渗盐水候选：infusionDrug 中真正含 3% / hypertonic 的药名", "
  SELECT drugname, COUNT(DISTINCT patientunitstayid) AS units FROM infusionDrug
  WHERE LOWER(drugname) LIKE '%hypertonic%' OR LOWER(drugname) LIKE '%3% saline%'
     OR LOWER(drugname) LIKE '%3% nacl%'     OR LOWER(drugname) LIKE '%3% sodium chloride%'
     OR LOWER(drugname) LIKE '%3%ns%'        OR LOWER(drugname) LIKE '%3% ns%'
  GROUP BY 1 ORDER BY units DESC LIMIT 20")
cat("\n已写入（探针结束）\n")
