## =============================================================================
## _run11.R —— 补充分析：严重度梯度、绝对风险差、E 值（两库各跑一次）
##
## 只读已落盘的提取与队列结果，不动主流程，不重跑队列组装。
## 两者共用 export.R 之外的独立结果对象，出表由 _run12.R 负责。
## =============================================================================
source("extract.R")
source("common.R")
source("models.R")
source("advance.R")

DB <- list(
  list(名 = "MIMIC-IV", ex = "mimic_extract_B.rds",
       full = "mimic_B_full.rds", out = "mimic_advance.rds"),
  list(名 = "eICU-CRD", ex = "eicu_B_extract.rds",
       full = "eicu_C_full.rds", out = "eicu_C_advance.rds")
)

pct <- function(x, d = 1) sprintf(paste0("%.", d, "f%%"), 100 * x)
pt_ <- function(x, d = 1) sprintf(paste0("%.", d, "f 个百分点"), 100 * x)

for (s in DB) {
  cat("\n================ ", s$名, " ================\n", sep = "")
  a <- readRDS(file.path(OUT_DIR, s$full))$analysis
  ex <- readRDS(file.path(OUT_DIR, s$ex))

  ## 与主分析同一个分母：先取完整病例，再并纯肌酐分期
  a <- a[model_cases(a), , drop = FALSE]
  a <- attach_stage(a, ex)
  cat(sprintf("分析人群 %d，事件 %d（%s）\n", nrow(a), sum(a$aki), pct(mean(a$aki))))

  adv <- run_advance(a)

  cat("\n--- 严重度梯度（暴露：ΔNa 每升高 1 mmol/L，基础模型）---\n")
  for (i in seq_len(nrow(adv$严重度)))
    cat(sprintf("  %-42s %4d/%4d  OR %.3f（%.3f–%.3f）  P = %s\n",
                adv$严重度$结局[i], adv$严重度$事件数[i], adv$严重度$n[i],
                adv$严重度$OR[i], adv$严重度$下限[i], adv$严重度$上限[i],
                format.pval(adv$严重度$P[i], digits = 2, eps = 1e-4)))

  cat("\n--- 绝对风险差（边际标准化，bootstrap B = 200）---\n")
  for (nm in names(adv$风险)) {
    r <- adv$风险[[nm]]
    cat(sprintf("  ΔNa %-8s 风险 %s → %s，差 %s（%s – %s）\n", nm,
                pct(r$from), pct(r$to), pt_(r$diff),
                pt_(r$diff_ci[1]), pt_(r$diff_ci[2])))
  }

  cat("\n--- E 值 ---\n")
  for (i in seq_len(nrow(adv$E值)))
    cat(sprintf("  %-26s %-24s %s\n", adv$E值$对比[i], adv$E值$量[i], adv$E值$值[i]))

  saveRDS(adv, file.path(OUT_DIR, s$out))
  cat(sprintf("\n已存 %s\n", file.path(OUT_DIR, s$out)))
}
