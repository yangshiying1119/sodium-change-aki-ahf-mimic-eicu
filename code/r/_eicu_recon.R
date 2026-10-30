source("extract.R")
root <- EICU_ROOT
con <- dbConnect(duckdb())
q <- function(t, extra="") dbGetQuery(con, sprintf(
  "SELECT * FROM read_csv_auto('%s/%s.csv.gz', header=true, sample_size=-1) %s LIMIT 0",
  root, t, extra))
for (t in c("patient","lab","intakeOutput","nurseCharting","medication",
            "treatment","diagnosis","infusionDrug","vitalAperiodic","hospital")) {
  cols <- tryCatch(names(q(t)), error=function(e) paste("ERR:", conditionMessage(e)))
  cat(sprintf("\n== %s ==\n%s\n", t, paste(cols, collapse=", ")))
}
dbDisconnect(con, shutdown=TRUE)
