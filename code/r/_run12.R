## =============================================================================
## _run12.R —— 出表：表 17（严重度梯度）、表 18（绝对风险差与 E 值）、表 19（阴性对照）
##
## 与 export.R 同一套落盘约定：本文件不做统计，只把 _run11 / _run13 算好的
## 结果对象整理成表；中文与英文各写一套目录，靠 i18n.R 的 tr_df() 翻译。
## =============================================================================
source("extract.R")
source("common.R")
source("i18n.R")
source("export.R")

.fmt_num <- function(x, d = 1) {
  ## 负号统一用 U+2212，与正文一致（ASCII 连字符在 Word 里宽度不匹配）
  s <- sprintf(paste0("%.", d, "f"), x)
  gsub("-", "−", s, fixed = TRUE)
}

.fmt_pct <- function(x, d = 1) paste0(.fmt_num(100 * x, d), "%")

DB <- list(
  list(adv = "mimic_advance.rds",  neg = "mimic_negctl.rds",
       dir = file.path(OUT_DIR, out_subdir("tables", "tables_en"))),
  list(adv = "eicu_C_advance.rds", neg = "eicu_C_negctl.rds",
       dir = file.path(OUT_DIR, out_subdir("tables_eicu", "tables_eicu_en")))
)

for (s in DB) {
  adv <- readRDS(file.path(OUT_DIR, s$adv))

  ## ---- 表 17 严重度梯度 ----
  t17 <- data.frame(
    结局   = adv$严重度$结局,
    人数   = adv$严重度$n,
    事件数 = adv$严重度$事件数,
    `OR（95% CI）` = .fmt_or(adv$严重度$OR, adv$严重度$下限, adv$严重度$上限),
    P = .fmt_p(adv$严重度$P),
    check.names = FALSE, row.names = NULL)

  ## ---- 表 18 绝对风险差与 E 值 ----
  r03 <- adv$风险[["0 vs +3"]]
  r33 <- adv$风险[["-3 vs +3"]]
  e1  <- adv$E值[adv$E值$对比 == "ΔNa 每升高 1 mmol/L", ]
  e3  <- adv$E值[adv$E值$对比 == "ΔNa 每升高 3 mmol/L", ]
  ev  <- function(d, 量) d$值[d$量 == 量]

  t18 <- data.frame(
    项目 = c(
      "标准化 72 h AKI 风险（ΔNa = 0）",
      "标准化 72 h AKI 风险（ΔNa = +3）",
      "绝对风险差（+3 vs 0 mmol/L）",
      "绝对风险差（+3 vs −3 mmol/L）",
      "OR（基础模型，每 1 mmol/L）",
      "OR（基础模型，每 3 mmol/L）",
      "E 值（每 1 mmol/L）",
      "E 值（每 3 mmol/L）"),
    估计值 = c(
      .fmt_pct(r03$from),
      .fmt_pct(r03$to),
      .fmt_pct(r03$diff),
      .fmt_pct(r33$diff),
      ev(e1, "OR（基础模型）"),
      ev(e3, "OR（基础模型）"),
      ev(e1, "E 值（点估计）"),
      ev(e3, "E 值（点估计）")),
    `95% CI` = c(
      sprintf("%s – %s", .fmt_pct(r03$from_ci[1]), .fmt_pct(r03$from_ci[2])),
      sprintf("%s – %s", .fmt_pct(r03$to_ci[1]),   .fmt_pct(r03$to_ci[2])),
      sprintf("%s – %s", .fmt_pct(r03$diff_ci[1]), .fmt_pct(r03$diff_ci[2])),
      sprintf("%s – %s", .fmt_pct(r33$diff_ci[1]), .fmt_pct(r33$diff_ci[2])),
      "", "", "", ""),
    check.names = FALSE, row.names = NULL)

  ## OR 的区间直接取结果对象里的 lo/hi，不从格式化串反解析
  orci <- function(k) sprintf("%.3f – %.3f",
                              adv$基础效应$下限^k, adv$基础效应$上限^k)
  t18[["95% CI"]][5] <- orci(1)
  t18[["95% CI"]][6] <- orci(3)
  t18[["95% CI"]][7] <- sprintf("CI 端点 %s", ev(e1, "E 值（95% CI 端点）"))
  t18[["95% CI"]][8] <- sprintf("CI 端点 %s", ev(e3, "E 值（95% CI 端点）"))

  ## ---- 表 19 阴性对照 ----
  t19 <- NULL
  negf <- file.path(OUT_DIR, s$neg)
  if (file.exists(negf)) {
    neg <- readRDS(negf)
    ## 结果对象里的暴露是列名（delta_na 等），直接进表会露出代码口径；
    ## 这里换成读者看得懂的符号写法，与正文一致。
    disp <- c(delta_na = "ΔNa", delta_k = "ΔK", delta_cl = "ΔCl",
              "delta_na（同时调整 delta_k）" = "ΔNa（同时调整 ΔK）",
              "delta_na（同时调整 delta_cl）" = "ΔNa（同时调整 ΔCl）",
              "delta_cl（同时调整 delta_na）" = "ΔCl（同时调整 ΔNa）",
              "delta_na（每 1 SD，调整 delta_cl）" = "ΔNa（每 1 SD，同时调整 ΔCl）",
              "delta_cl（每 1 SD，调整 delta_na）" = "ΔCl（每 1 SD，同时调整 ΔNa）")
    t19 <- data.frame(
      暴露   = unname(disp[neg$效应$暴露]),
      人数   = neg$效应$n,
      事件数 = neg$效应$事件数,
      `OR（95% CI）` = .fmt_or(neg$效应$OR, neg$效应$下限, neg$效应$上限),
      P = .fmt_p(neg$效应$P),
      check.names = FALSE, row.names = NULL)
  } else {
    cat("提示：未找到", s$neg, "，表 19 跳过\n")
  }

  ## 表号按正文首次提到的顺序排：绝对风险差在 3.9 节先出现，严重度梯度随后。
  for (nm in c("表17_绝对风险差与E值", "表18_严重度梯度",
               if (!is.null(t19)) "表19_阴性对照")) {
    d <- switch(nm,
                表17_绝对风险差与E值 = t18,
                表18_严重度梯度 = t17,
                表19_阴性对照 = t19)
    .write_csv(tr_df(d), tr_tabname(nm), s$dir)
    cat(sprintf("  已写 %s/%s.csv\n", basename(s$dir), tr_tabname(nm)))
  }
}
