## 外审稳线三项本地补充分析：间隔敏感性 + IPW 诊断 + eICU 校准
.args <- commandArgs(FALSE)
.sf <- sub("^--file=", "", grep("^--file=", .args, value = TRUE))
if (length(.sf)) setwd(dirname(normalizePath(.sf[1])))
suppressMessages({ source("common.R"); source("models.R"); source("cohort.R"); source("bias.R"); source("discrim.R") })
## 复现缓存目录：本地复现指向 repro 包的 code/out；换机器请改此处或设 SODIUM_OUT
OUT <- Sys.getenv("SODIUM_OUT", "D:/CC/repro_20261030/code/out")

fmt_or <- function(co, dig = 3) sprintf("%.3f (%.3f-%.3f), P=%s",
  exp(co[1]), exp(co[1]-1.96*co[2]), exp(co[1]+1.96*co[2]),
  format.pval(co[4], digits=3))

for (pair in list(c("MIMIC-IV", "mimic_B_full.rds", "mimic_B_models.rds"),
                  c("eICU-CRD", "eicu_C_full.rds", "eicu_C_models.rds"))) {
  tag <- pair[1]
  a <- readRDS(file.path(OUT, pair[2]))$analysis
  full_data <- readRDS(file.path(OUT, pair[2]))$data
  d <- a[model_cases(a), ]
  d$gap_h <- d$t_last - d$t_first
  cat(sprintf("\n========== %s (n=%d, AKI=%d) ==========\n", tag, nrow(d), sum(d$aki)))
  cat(sprintf("间隔 h: 中位 %.1f (IQR %.1f-%.1f)\n",
      median(d$gap_h), quantile(d$gap_h,.25), quantile(d$gap_h,.75)))

  .fx <- function(extra, xs, dat) stats::as.formula(
    paste("aki ~", paste(c(extra, xs), collapse = " + ")), env = parent.frame())

  ## ---- 分析1a: base + 间隔协变量 ----
  m <- glm(.fx(c("delta_na", "gap_h"), .pick(BASIC_X, d), d), data = d, family = binomial())
  cat(sprintf("[1a] base+gap_h:  ΔNa %s | gap_h %s\n",
      fmt_or(summary(m)$coeff["delta_na",]), fmt_or(summary(m)$coeff["gap_h",])))

  ## ---- 分析1b: 变化率暴露(每 0.1 mmol/L/h) ----
  d$rate10 <- d$delta_na / d$gap_h * 10
  m <- glm(.fx("rate10", .pick(BASIC_X, d), d), data = d, family = binomial())
  cat(sprintf("[1b] rate/0.1mmol/L/h: %s\n", fmt_or(summary(m)$coeff["rate10",])))

  ## ---- 分析1c: 间隔中位分层 ----
  med <- median(d$gap_h)
  for (lab in c("<=median", ">median")) {
    dd <- if (lab == "<=median") d[d$gap_h <= med, ] else d[d$gap_h > med, ]
    m <- glm(.fx("delta_na", .pick(BASIC_X, dd), dd), data = dd, family = binomial())
    cat(sprintf("[1c] gap %s (n=%d, ev=%d): %s\n", lab, nrow(dd), sum(dd$aki),
        fmt_or(summary(m)$coeff["delta_na",])))
  }

  ## ---- 分析2: IPW 权重诊断（须用全表型人群,同 _run7 口径）----
  ex_file <- if (tag == "eICU-CRD") "eicu_B_extract.rds" else "mimic_extract_B.rds"
  exx <- readRDS(file.path(OUT, ex_file))
  pop <- attach_covariates(full_data, exx$covariates)
  w <- fit_retest_ipw(pop)
  ww <- w$data$..w
  ess <- sum(ww)^2 / sum(ww^2)
  cat(sprintf("[2] IPW: n加权前=%d, n加权后=%.0f, ESS=%.0f (%.0f%%)\n",
      nrow(w$data), sum(ww), ess, 100*ess/nrow(w$data)))
  cat(sprintf("    截断点 %.3f-%.3f\n", w$截断点[1], w$截断点[2]))
  qq <- quantile(ww, c(.05, .25, .5, .75, .95))
  cat("    截断后权重分位数(5/25/50/75/95%):",
      paste(sprintf("%.3f", as.numeric(qq)), collapse = "/"), "\n")
  cat(sprintf("    截断前权重: min=%.4f 中位=%.4f max=%.4f\n",
      min(w$data$..w_raw), median(w$data$..w_raw), max(w$data$..w_raw)))
  cat(sprintf("    IPW: 复测 glm 收敛警告=%s（分离≈复测近乎普遍时属预期）\n",
      if (length(w$model$converged) && !w$model$converged) "是" else "否"))
  ip <- ips_retest_model(pop, w)
  cat(sprintf("    加权模型(稿件口径): OR %.3f (%.3f-%.3f), P=%s, n=%d, 事件=%d\n",
      ip$OR, ip$下限, ip$上限, format.pval(ip$P, digits=3), ip$n, ip$事件数))

  ## ---- 分析3: HL(基础模型, 与 discrim_all/稿件同口径) ----
  mb <- glm(.fx("delta_na", .pick(BASIC_X, d), d), data = d, family = binomial())
  hl <- hl_test(as.integer(d$aki), stats::fitted(mb))
  cat(sprintf("[3] HL(基础模型): chi2=%.2f, df=%d, P=%.3f (格子最小期望 %.1f)\n",
      hl$卡方, hl$自由度, hl$P, hl$格子最小期望))
}
