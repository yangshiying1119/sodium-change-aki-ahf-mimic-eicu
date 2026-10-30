source("extract.R")
con <- dbConnect(duckdb())
T_ <- function(t, quote=FALSE) sprintf(
  "read_csv_auto('%s/%s.csv.gz', header=true, sample_size=2048%s)",
  EICU_ROOT, t, if (quote) ", quote='\"'" else "")
for (t in c("patient","admissionDx","diagnosis","lab","intakeOutput",
            "nurseCharting","vitalAperiodic","vitalPeriodic","medication",
            "infusionDrug","treatment","respiratoryCare","apacheApsVar")) {
  q <- if (t %in% c("nurseCharting","respiratoryCare")) TRUE else FALSE
  nm <- tryCatch(names(dbGetQuery(con, sprintf("SELECT * FROM %s LIMIT 0", T_(t, q)))),
                 error=function(e) paste("ERR:", conditionMessage(e)))
  cat("\n## ", t, "\n", paste(nm, collapse=", "), "\n", sep="")
}
