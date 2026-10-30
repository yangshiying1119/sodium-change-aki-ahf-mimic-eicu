## eICU-CRD 2.0 侧正式跑批：一次扫库 → 队列 + 协变量 → 基线表 → 落盘
##
## 与 MIMIC 侧 _run2.R 的唯一区别是提取函数（extract_eicu_all）与落盘文件名。
## 队列组装、逐层留存、基线表全部走同一套共享代码——这正是方案 6.4 要的
## "跨库口径漂移"对策：差异被结构性地限制在提取适配层。
source("cohort.R")
source("tables.R")

## 口径由命令行给定（A/B/C），缺省 B。三种口径共用同一次提取——
## 提取层一次就把 cohort_A/B/C 三张表都建好了，换口径不需要重扫库。
args <- commandArgs(trailingOnly = TRUE)
CC <- if (length(args) >= 1 && nzchar(args[1])) toupper(args[1]) else "B"
stopifnot(CC %in% c("A", "B", "C"))
cat(sprintf("eICU 口径：表型 %s\n", CC))

## 只统计 **eICU 侧**的 SQL 作缓存新鲜度判据，避免改 MIMIC 的脚本时白扫一遍 eICU 库。
## 缓存与口径无关，文件名固定为 eicu_B_extract.rds。
cache <- file.path(OUT_DIR, "eicu_B_extract.rds")
sqls  <- list.files(SQL_DIR, pattern = "_eicu\\.sql$", full.names = TRUE)
if (file.exists(cache) && file.mtime(cache) > max(file.mtime(sqls))) {
  cat(sprintf("复用提取缓存（%s）\n", format(file.mtime(cache), "%H:%M:%S")))
  ex <- readRDS(cache)
} else {
  ex <- extract_eicu_all(CC)
  saveRDS(ex, cache)
  cat(sprintf("提取结果已缓存至 %s\n", cache))
}

## 各表型口径的规模先摆出来，便于确认 A/B/C 的嵌套关系（A ⊆ B ⊆ C）
for (cc in c("A", "B", "C"))
  if (!is.null(ex$cohorts[[cc]]))
    cat(sprintf("  表型 %s：%d 个 ICU 单元\n", cc, nrow(ex$cohorts[[cc]])))

## 协变量体检：把"整列取不到数"这类静默失败挡在建模之前。
## 放在队列组装**之前**——若某个提取谓词没命中（标签写错、join 不中），
## 现在停下来改 SQL 只要重扫一次；等到模型跑出全 NA 系数再回头查，代价大得多。
cat("\n=== 协变量体检 ===\n")
print(as.data.frame(check_covariates(ex$covariates, db = "eicu")))

res <- run_cohort(CC, ex = ex)
print_attrition(res)

a <- attach_covariates(res$analysis, ex$covariates)
cat(sprintf("\n=== 协变量覆盖 ===\n  分析人群 %d 行，并上协变量后 %d 行\n",
            nrow(res$analysis), nrow(a)))
miss <- sapply(setdiff(names(ex$covariates),
                       c("stay_id", "hadm_id", "subject_id", "id")),
               function(v) sum(is.na(a[[v]])))
miss <- sort(miss[miss > 0], decreasing = TRUE)
if (length(miss)) {
  cat("  有缺失的协变量（缺失数）：\n")
  for (i in seq_along(miss)) cat(sprintf("    %-28s %5d（%.1f%%）\n",
      names(miss)[i], miss[i], 100 * miss[i] / nrow(a)))
} else cat("  无缺失\n")

cat("\n=== 基线表（按 ΔNa 三档）===\n")
print_baseline(baseline_table(a))

out_full <- file.path(OUT_DIR, sprintf("eicu_%s_full.rds", CC))
saveRDS(list(analysis = a, data = res$data, attrition = res$attrition,
             sodium = res$sodium, creatinine = res$creatinine,
             na_df = res$na_df, rrt = res$rrt,
             rrt_pre = res$rrt_pre, rrt_post = res$rrt_post,
             covariates = ex$covariates),
        out_full)
cat(sprintf("\n已存 %s\n", out_full))
