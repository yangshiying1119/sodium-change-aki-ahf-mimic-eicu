## 尿路 itemid 勘察：outputevents 里到底有哪些出量类型
suppressMessages({ library(DBI); library(duckdb) })
source("extract.R")
con <- dbConnect(duckdb())
on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
run_sql_file(con, "00_schema_mimic.sql", list(MIMIC_ROOT = MIMIC_ROOT))

r <- dbGetQuery(con, "
  SELECT o.itemid, d.label, d.category, COUNT(*) AS n,
         ROUND(SUM(o.value)) AS tot_mL,
         ROUND(MEDIAN(o.value), 1) AS med_value
  FROM outputevents o JOIN d_items d USING (itemid)
  GROUP BY 1, 2, 3 ORDER BY n DESC")
print(r, row.names = FALSE)
