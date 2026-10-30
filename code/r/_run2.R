## 正式跑批：一次扫库 → 队列 + 协变量 → 基线表 → 落盘
## chartevents 只扫一遍（extract_mimic_all），比分开跑省一半时间。
source("cohort.R")
source("tables.R")

## 提取层缓存：扫库约 8 分钟，而后续建模、敏感性分析只需重新组装，不需重扫。
## 任一 .sql 比缓存新就重跑，避免改了 SQL 却读旧缓存。
## 只统计 **MIMIC 侧**的 SQL：eICU 的 .sql 与这份缓存无关，若一并计入，
## 每写一个 eICU 脚本都会把 MIMIC 的缓存判为过期，白扫一遍库。
cache <- file.path(OUT_DIR, "mimic_extract_B.rds")
sqls  <- list.files(SQL_DIR, pattern = "_mimic\\.sql$", full.names = TRUE)
if (file.exists(cache) && file.mtime(cache) > max(file.mtime(sqls))) {
  cat(sprintf("复用提取缓存（%s）\n", format(file.mtime(cache), "%H:%M:%S")))
  ex <- readRDS(cache)
} else {
  ex <- extract_mimic_all("B")
  saveRDS(ex, cache)
  cat(sprintf("提取结果已缓存至 %s\n", cache))
}

res <- run_mimic_cohort("B", ex = ex)
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

saveRDS(list(analysis = a, data = res$data, attrition = res$attrition,
             sodium = res$sodium, creatinine = res$creatinine,
             na_df = res$na_df, rrt = res$rrt,
             rrt_pre = res$rrt_pre, rrt_post = res$rrt_post,
             covariates = ex$covariates),
        file.path(OUT_DIR, "mimic_B_full.rds"))
cat(sprintf("\n已存 %s\n", file.path(OUT_DIR, "mimic_B_full.rds")))
