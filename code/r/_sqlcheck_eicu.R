## 空表干跑：把每张源表建成 **0 行** 的视图后按顺序执行全部 eICU SQL。
##
## 目的：在付出长扫描代价之前，把"列名写错 / 语法错 / join 键不存在"这类错误
## 一次性暴露出来。DuckDB 的 Binder 在**准备阶段**就会解析每一条语句引用的列，
## 与表里有没有数据无关——所以 0 行输入足以校验全部列引用与语法。
## 谓词命中率、覆盖度等**数据层面**的问题它查不出，那是随后 check_covariates()
## 与基线表的职责，两者互补。
##
## ⚠️ 上一版这里手写了一份 15 张表的清单，结果清单里有 respiratoryCare 而
##    00_schema_eicu.sql 里没有——干跑通过，真跑报 "Table respiratoryCare does
##    not exist"。手写清单本身就是 bug 源。现在改为：
##      (1) 视图清单从 00_schema_eicu.sql **推导**，不再手抄；
##      (2) 静态扫描全部 eICU SQL 的 FROM / JOIN 表名，逐一断言其来源可解释
##          （schema 视图，或同批更早脚本建出的表）。这一条能**静态**抓住
##          "引用了没人建的表"，而 (1) 只保证视图清单不脱节。
##
## 本脚本不落任何分析产物，仅用于提交长跑之前的自检。
source("extract.R")
suppressMessages(library(DBI))

SQL_ORDER <- c("10_hf_phenotype_eicu.sql", "20_labs_eicu.sql", "30_cohort_eicu.sql",
               "21_urine_eicu.sql", "42_rrt_eicu.sql", "40_covariates_eicu.sql")

.read <- function(f) paste(readLines(file.path(SQL_DIR, f), warn = FALSE), collapse = "\n")

## ---------- (1) 从 00_schema_eicu.sql 推导视图清单 ----------
schema <- .read("00_schema_eicu.sql")
blocks <- strsplit(schema, "CREATE OR REPLACE VIEW", fixed = TRUE)[[1]][-1]
views <- lapply(blocks, function(b) {
  nm  <- sub("^\\s*([A-Za-z_][A-Za-z0-9_]*).*", "\\1", b)
  csv <- regmatches(b, regexpr("[A-Za-z_][A-Za-z0-9_]*\\.csv\\.gz", b))
  if (!length(csv)) NULL else list(name = nm, csv = csv)
})
views <- Filter(Negate(is.null), views)
stopifnot(length(views) > 0)
## 视图名必须唯一，否则后面的 CREATE OR REPLACE 会互相覆盖而不报错
vnames <- vapply(views, `[[`, "", "name")
if (anyDuplicated(vnames)) stop("00_schema_eicu.sql 中有重复视图名：",
                                paste(vnames[duplicated(vnames)], collapse = ", "))

con <- dbConnect(duckdb())
on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
for (v in views)
  dbExecute(con, sprintf(
    "CREATE OR REPLACE VIEW %s AS SELECT * FROM read_csv_auto('%s/%s', header = true, escape = '\"') LIMIT 0",
    v$name, EICU_ROOT, v$csv))
cat(sprintf("依据 00_schema_eicu.sql 建了 %d 张 0 行视图\n", length(views)))

## ---------- (2) 静态扫描：每个 FROM / JOIN 的表名都必须有来源 ----------
## 逐文件、按语句先后累积"已建出的表"，同批更早脚本建的表对后续脚本可见。
created <- vnames
cte_re  <- "(?i)\\bwith\\s+([A-Za-z_][A-Za-z0-9_]*)\\s+as\\b|,\\s*([A-Za-z_][A-Za-z0-9_]*)\\s+as\\s*\\("
src_re  <- "(?i)\\b(?:from|join)\\s+([A-Za-z_][A-Za-z0-9_]*)"
made_re <- "(?i)create\\s+or\\s+replace\\s+(?:table|view)\\s+([A-Za-z_][A-Za-z0-9_]*)"

unknown <- list()
for (s in SQL_ORDER) {
  txt  <- .read(s)
  ## 本文件建出的表先并入：文件内部"先建后用"的顺序由第 (3) 步真跑兜住，
  ## 静态扫描只管**跨文件**的存在性。
  created <- unique(c(created, gsub(made_re, "\\1",
                 unlist(regmatches(txt, gregexpr(made_re, txt, perl = TRUE))), perl = TRUE)))
  ctes <- unlist(regmatches(txt, gregexpr(cte_re, txt, perl = TRUE)))
  ctes <- gsub(cte_re, "\\1\\2", ctes, perl = TRUE)
  ctes <- ctes[nzchar(ctes)]
  refs <- unique(gsub(src_re, "\\1", unlist(regmatches(txt, gregexpr(src_re, txt, perl = TRUE))), perl = TRUE))
  refs <- setdiff(refs, c(ctes, "read_csv_auto", "unnest", "values"))
  bad  <- setdiff(refs, created)
  if (length(bad)) unknown[[s]] <- bad
}
if (length(unknown)) {
  cat("\n✗ 以下表被引用但没有任何脚本建出：\n")
  for (s in names(unknown)) cat(sprintf("  %s：%s\n", s, paste(unknown[[s]], collapse = ", ")))
  cat("\n  ——若它是 eICU 源表，请在 00_schema_eicu.sql 里补视图；\n")
  cat("    若是中间表，请确认建表语句在引用它的脚本之前。\n")
  quit(status = 1)
}
cat("静态扫描通过：全部 FROM / JOIN 表名均可追溯到 schema 视图或更早的建表语句\n")

## ---------- (3) 真跑一遍解析与列名绑定 ----------
## 00_schema_eicu.sql 已被上面的推导取代（此处视图为 0 行，与 schema 同名同列），
## 故不重复执行它，只跑其余六个脚本。
for (s in SQL_ORDER) {
  ok <- tryCatch({ run_sql_file(con, s); TRUE },
                 error = function(e) { cat(sprintf("  ✗ %s\n    %s\n", s, conditionMessage(e))); FALSE })
  if (!ok) { cat("\n干跑失败，先修 SQL 再提交长跑。\n"); quit(status = 1) }
}
cat("\n全部 eICU SQL 在 0 行输入下均通过解析与列名绑定。\n")
