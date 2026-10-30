## =============================================================================
## extract.R —— 提取驱动（两库共用，只换根目录与 SQL 文件）
##
## 职责：把 sql/ 下的脚本在本库上按顺序执行，产出统一的 (id, t_hours, value)
## 长表与队列表，交给 common.R 的判定层。本文件不含任何判定逻辑。
## =============================================================================

suppressMessages({ library(DBI); library(duckdb) })

.here <- local({
  a <- commandArgs(trailingOnly = FALSE)
  s <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(s) && nzchar(s)) dirname(normalizePath(s)) else "."
})

MIMIC_ROOT <- "/path/to/mimic-iv-3.1"
EICU_ROOT  <- file.path(.here, "..", "..", "data", "eicu-crd-2.0")
SQL_DIR    <- file.path(.here, "..", "sql")
OUT_DIR    <- file.path(.here, "..", "out")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

#' 执行一个 .sql 文件；{VAR} 占位符由 vars 替换
#' 按行尾分号切分语句（脚本内不出现字符串内含分号的情形）
run_sql_file <- function(con, file, vars = list()) {
  path <- file.path(SQL_DIR, file)
  stopifnot(file.exists(path))
  sql <- paste(readLines(path, warn = FALSE), collapse = "\n")
  ## fixed = TRUE 时模式按字面处理，占位符直接写 {NAME}，不要再加反斜杠转义
  for (nm in names(vars)) sql <- gsub(paste0("{", nm, "}"), vars[[nm]], sql, fixed = TRUE)
  sql <- gsub("(?m)^\\s*--.*$", "", sql, perl = TRUE)      # 去整行注释（多行模式）
  stmts <- strsplit(sql, ";", fixed = TRUE)[[1]]
  n <- 0L
  for (s in stmts) {
    s <- trimws(s)
    if (!nzchar(s)) next
    dbExecute(con, s); n <- n + 1L
  }
  cat(sprintf("  %-28s 执行 %d 条语句\n", file, n))
  invisible(n)
}

#' 建库连接并跑完提取，返回长表与队列
#'
#' 41_chartevents 必须跑：RRT 事件（排除标准③与次要结局）要靠它取的 227290，
#' 且它是 chartevents 唯一的一次扫描，后续都从 ce_slim 取数。
#' 取出三个表型口径的人群表
#'
#' **必须三个都返回**：S6（改用表型定义 A）与 _runA.R 都只在组装层切换口径。
#' 若只带回被请求的那一个，`run_mimic_cohort(cohort = "A", ex = ex)` 会**静默**
#' 沿用 ex$cohort（即 B 口径）的人群，敏感性分析退化成主分析的复制品
#' ——实测 S6 曾输出与主分析逐位相同的 2700/718/0.933。
.all_cohorts <- function(con) {
  out <- list()
  for (cc in c("A", "B", "C")) {
    t <- tryCatch(dbGetQuery(con, sprintf("SELECT * FROM cohort_%s", cc)),
                  error = function(e) NULL)
    if (!is.null(t)) { t$id <- as.character(t$stay_id); out[[cc]] <- t }
  }
  stopifnot(length(out) >= 1L)
  out
}

extract_mimic <- function(cohort = "B", root = MIMIC_ROOT) {
  con <- dbConnect(duckdb())
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
  t0 <- Sys.time()

  run_sql_file(con, "00_schema_mimic.sql", list(MIMIC_ROOT = root))
  run_sql_file(con, "20_labs_mimic.sql")
  run_sql_file(con, "30_cohort_mimic.sql")
  run_sql_file(con, "41_chartevents_mimic.sql")
  run_sql_file(con, "42_rrt_mimic.sql")

  long <- function(tbl) {
    d <- dbGetQuery(con, sprintf("SELECT id, t_hours, value FROM %s", tbl))
    d$id <- as.character(d$id)
    d
  }
  coh <- dbGetQuery(con, sprintf("SELECT * FROM cohort_%s", cohort))
  coh$id <- as.character(coh$stay_id)

  ## 节点 N0–N2 的原始计数；N3–N8 依赖化验序列，在 R 侧完成
  cnt <- function(q) dbGetQuery(con, q)[[1]]
  attrition <- c(
    N0_icu_stays_all = cnt("SELECT COUNT(*) FROM icustays"),
    N1_adult         = cnt("SELECT COUNT(*) FROM base"),
    N2_hf_phenotype  = nrow(coh)
  )

  cat(sprintf("提取完成，用时 %.1f 秒\n", as.numeric(Sys.time() - t0, units = "secs")))
  list(db = "mimic",
       sodium = long("sodium"), sodium_blood = long("sodium_blood"),
       creatinine = long("creatinine"), rrt = long("rrt"),
       esrd_chart = dbGetQuery(con, "SELECT DISTINCT CAST(stay_id AS VARCHAR) AS id FROM esrd_chart")$id,
       cohort = coh, cohorts = .all_cohorts(con),
       attrition = attrition, con = NULL)
}

#' 单独跑一遍协变量层（扩展模型用）
#'
#' 复用 41 建出的 ce_slim：41→42→43 顺序不能乱，43 的汇总要 join 42 的产物，
#' 42 又要读 41 的 ce_slim。
extract_mimic_covariates <- function(root = MIMIC_ROOT) {
  con <- dbConnect(duckdb())
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
  t0 <- Sys.time()

  run_sql_file(con, "00_schema_mimic.sql", list(MIMIC_ROOT = root))
  run_sql_file(con, "30_cohort_mimic.sql")
  run_sql_file(con, "41_chartevents_mimic.sql")
  run_sql_file(con, "42_rrt_mimic.sql")
  run_sql_file(con, "43_covariates_mimic.sql")

  cv <- dbGetQuery(con, "SELECT * FROM cov_all")
  cv$id <- as.character(cv$stay_id)
  cat(sprintf("协变量提取完成，用时 %.1f 秒；%d 行 / 覆盖 cohort_C\n",
              as.numeric(Sys.time() - t0, units = "secs"), nrow(cv)))
  cv
}

#' 协变量体检：把"整列取不到数"这类静默失败变成硬停
#'
#' 本项目反复踩的坑是**同一类**：某个构造成分因为列名写错 / 标签不匹配 / join 不中，
#' 悄无声息地返回全 NA 或全 0，而后续代码不报错——
#'   · `glm` 遇 NA 按行删除，全 NA 的协变量会把**整个分析人群删空**；
#'   · 全 0 的二分类协变量只是让系数变成 NA，模型照跑，表里多一个"NA (NA–NA)"。
#' 两者都能一路混进成稿。故在建模之前统一体检，宁可在这里失败。
#'
#' @param cv   协变量表（cov_all 的返回值）
#' @param db   库名，只用于报错信息
#' @param expect 期望存在的列；缺一即停。默认取 BASIC_X/EXTEND_X 中**由 SQL 侧
#'        构造**的那些（age/gender 在 base、na_first/creat_baseline 在判定层）。
#' @return 不可见地返回覆盖率表（变量、非缺失数、不同取值数）
check_covariates <- function(cv, db = "mimic", expect = NULL) {
  if (is.null(expect))
    expect <- c("ckd", "diabetes", "hypertension", "chd", "atrial_fib", "valvular",
                "sbp_min_24h", "hr_max_24h", "rr_max_24h", "gcs_min_24h",
                "mech_vent_any", "vasoactive_any", "nee_max",
                "diuretic_furosemide_eq_mg", "fluid_in_24h_ml", "urine_out_24h_ml",
                "glucose_max_24h", "hgb_min_24h")
  expect <- expect[expect %in% c(names(cv), expect)]     # 保留原样，下面逐条判
  miss_col <- setdiff(expect, names(cv))
  if (length(miss_col))
    stop(sprintf("[%s] 协变量表缺少列：%s", db, paste(miss_col, collapse = ", ")))

  n <- nrow(cv)
  cov_tab <- data.frame(
    变量 = expect,
    非缺失 = vapply(expect, function(v) sum(!is.na(cv[[v]])), 0L),
    不同取值 = vapply(expect, function(v) length(unique(cv[[v]][!is.na(cv[[v]])])), 0L),
    row.names = NULL)
  cov_tab$覆盖率 <- sprintf("%.1f%%", 100 * cov_tab$非缺失 / n)

  bad <- cov_tab$不同取值 <= 1L
  if (any(bad))
    stop(sprintf("[%s] 下列协变量在全部 %d 行上取不到有效值（全 NA 或恒定），
说明其提取谓词没命中任何记录，模型会因此退化或静默删空人群：
%s", db, n,
      paste(sprintf("  %s（非缺失 %d，不同取值 %d）",
                    cov_tab$变量[bad], cov_tab$非缺失[bad], cov_tab$不同取值[bad]),
            collapse = "\n")))
  invisible(cov_tab)
}

#' 一次扫库同时拿到队列长表与协变量（正式跑批用）
#'
#' 与分别调用 extract_mimic() + extract_mimic_covariates() 等价，但 chartevents
#' 只扫一遍（那是全流程最贵的一步），并多返回一份 covariates。
#' 返回值多一个 covariates 元素，其余与 extract_mimic() 完全一致。
extract_mimic_all <- function(cohort = "B", root = MIMIC_ROOT) {
  con <- dbConnect(duckdb())
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
  t0 <- Sys.time()

  run_sql_file(con, "00_schema_mimic.sql", list(MIMIC_ROOT = root))
  run_sql_file(con, "20_labs_mimic.sql")
  run_sql_file(con, "30_cohort_mimic.sql")
  run_sql_file(con, "21_urine_mimic.sql")      # 尿量长表（S8）。须在 30 之后（用 cohort_C）
  run_sql_file(con, "41_chartevents_mimic.sql")
  run_sql_file(con, "42_rrt_mimic.sql")
  run_sql_file(con, "43_covariates_mimic.sql")

  long <- function(tbl) {
    d <- dbGetQuery(con, sprintf("SELECT id, t_hours, value FROM %s", tbl))
    d$id <- as.character(d$id)
    d
  }
  coh <- dbGetQuery(con, sprintf("SELECT * FROM cohort_%s", cohort))
  coh$id <- as.character(coh$stay_id)

  cnt <- function(q) dbGetQuery(con, q)[[1]]
  attrition <- c(
    N0_icu_stays_all = cnt("SELECT COUNT(*) FROM icustays"),
    N1_adult         = cnt("SELECT COUNT(*) FROM base"),
    N2_hf_phenotype  = nrow(coh)
  )

  cv <- dbGetQuery(con, "SELECT * FROM cov_all")
  cv$id <- as.character(cv$stay_id)
  ## 体重单独带出：S8 的尿量速率要按 kg 标化，而体重在协变量表里是可能缺失的单列，
  ## 不宜从 cov_all 反查后混进建模用的人群表
  wt <- dbGetQuery(con, "SELECT CAST(stay_id AS VARCHAR) AS id, weight_kg FROM cov_weight")
  weight <- setNames(wt$weight_kg, wt$id)

  cat(sprintf("提取完成（含协变量），用时 %.1f 秒；协变量 %d 行\n",
              as.numeric(Sys.time() - t0, units = "secs"), nrow(cv)))
  list(db = "mimic",
       sodium = long("sodium"), sodium_blood = long("sodium_blood"),
       creatinine = long("creatinine"), glucose = long("glucose"),
       ## 钾 / 氯只用于阴性对照暴露（advance.R），不进主分析
       potassium = long("potassium"), chloride = long("chloride"),
       rrt = long("rrt"),
       urine = long("urine_ts"),
       esrd_chart = dbGetQuery(con, "SELECT DISTINCT CAST(stay_id AS VARCHAR) AS id FROM esrd_chart")$id,
       cohort = coh, cohorts = .all_cohorts(con),
       attrition = attrition, covariates = cv, weight = weight, con = NULL)
}

## =============================================================================
## eICU-CRD 2.0 侧
##
## 与 MIMIC 侧的**唯一区别在提取**：SQL 文件不同、根目录不同。装配层的约定完全
## 一致——cohort_* 里 patientunitstayid 已别名为 stay_id，故 .all_cohorts() 与
## r/cohort.R 的判定层可原样复用，一行都不用改（方案 6.4 的"跨库口径漂移"对策）。
## =============================================================================

#' 一次扫库拿到 eICU 的队列长表与协变量
extract_eicu_all <- function(cohort = "B", root = EICU_ROOT) {
  con <- dbConnect(duckdb())
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
  t0 <- Sys.time()

  run_sql_file(con, "00_schema_eicu.sql", list(EICU_ROOT = root))
  run_sql_file(con, "10_hf_phenotype_eicu.sql")   # hf_admitdx / hf_dx_flags / excl_esrd
  run_sql_file(con, "20_labs_eicu.sql")           # sodium / creatinine / glucose
  run_sql_file(con, "30_cohort_eicu.sql")         # base / landmark / esrd_chart / cohort_A|B|C
  run_sql_file(con, "21_urine_eicu.sql")          # 尿量长表（S8）。须在 30 之后（用 cohort_C）
  run_sql_file(con, "42_rrt_eicu.sql")
  run_sql_file(con, "40_covariates_eicu.sql")

  long <- function(tbl) {
    d <- dbGetQuery(con, sprintf("SELECT id, t_hours, value FROM %s", tbl))
    d$id <- as.character(d$id)
    d
  }
  coh <- dbGetQuery(con, sprintf("SELECT * FROM cohort_%s", cohort))
  coh$id <- as.character(coh$stay_id)

  cnt <- function(q) dbGetQuery(con, q)[[1]]
  attrition <- c(
    N0_icu_stays_all = cnt("SELECT COUNT(*) FROM patient"),
    N1_adult         = cnt("SELECT COUNT(*) FROM base"),
    N2_hf_phenotype  = nrow(coh)
  )

  cv <- dbGetQuery(con, "SELECT * FROM cov_all")
  cv$id <- as.character(cv$stay_id)
  wt <- dbGetQuery(con, "SELECT CAST(stay_id AS VARCHAR) AS id, weight_kg FROM cov_weight")
  weight <- setNames(wt$weight_kg, wt$id)

  ## ⚠️ eICU 本库**没有血气钠**这一路（lab 里 sodium 只有血清/血浆一种，
  ##    urinary sodium 是尿钠，不是血气）。故作敏感性分析的 S4/S4b 在 eICU 侧
  ##    不适用；返回零行表而非 NULL，让 sensitivity.R 的 `nrow() > 0` 判据自然跳过，
  ##    而不是把一张空表当成"血气钠全是 NA"去做回归。
  sodium_blood <- data.frame(id = character(0), t_hours = numeric(0), value = numeric(0))

  cat(sprintf("eICU 提取完成，用时 %.1f 秒；协变量 %d 行\n",
              as.numeric(Sys.time() - t0, units = "secs"), nrow(cv)))
  list(db = "eicu",
       sodium = long("sodium"), sodium_blood = sodium_blood,
       creatinine = long("creatinine"), glucose = long("glucose"),
       ## 钾 / 氯只用于阴性对照暴露（advance.R），不进主分析
       potassium = long("potassium"), chloride = long("chloride"),
       rrt = long("rrt"),
       urine = long("urine_ts"),
       esrd_chart = dbGetQuery(con, "SELECT DISTINCT CAST(stay_id AS VARCHAR) AS id FROM esrd_chart")$id,
       cohort = coh, cohorts = .all_cohorts(con),
       attrition = attrition, covariates = cv, weight = weight, con = NULL)
}
