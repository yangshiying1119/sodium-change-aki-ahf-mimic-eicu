source("extract.R")
con <- dbConnect(duckdb()); on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
run_sql_file(con, "00_schema_mimic.sql", list(MIMIC_ROOT = MIMIC_ROOT))
run_sql_file(con, "30_cohort_mimic.sql")
q <- function(s) dbGetQuery(con, s)[[1]]
cat("\n=== 与勘察值对照 ===\n")
cat(sprintf("  base 成人+>24h      : %6d   (勘察 74827)\n", q("select count(*) from base")))
cat(sprintf("  A 主诊断急性心衰    : %6d   (勘察  1367)\n", q("select count(*) from cohort_A")))
cat(sprintf("  B 前3位+急性码      : %6d   (勘察  7411)  <- 主分析\n", q("select count(*) from cohort_B")))
cat(sprintf("  C 任意位 I50/428    : %6d   (勘察 21549)\n", q("select count(*) from cohort_C")))
cat(sprintf("  仅任意位含急性码    : %6d   (勘察 11620)\n",
    q("select count(*) from base b join hf_flags f using(stay_id) where f.acute_hf=1")))
