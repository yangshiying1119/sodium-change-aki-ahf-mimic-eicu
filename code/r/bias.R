## =============================================================================
## bias.R —— 选择与测量偏倚的专门分析（方案 5.6 第四步 / 5.7）
##
## 方案 5.7 要求两件事，本文件各实现一件：
##   ① 比较"可判定结局者"与"不可判定者"的基线特征、ΔNa 分布与血钠检测次数。
##      没有复测不等于没有 AKI，所以差异有多大、朝哪个方向，必须写清楚。
##   ② 复测概率权重：用界标前协变量建"界标后是否复测肌酐"的预测模型，以逆概率
##      加权重跑主模型，并报告权重截断规则与权重分布。
##
## 与 sensitivity.R 的分工：那里是"换一个口径看结论稳不稳"，这里是"看人群本身
## 被怎样筛选过"。两者的共同纪律是——**绝不另写一套筛选逻辑**，一律复用
## cohort.R 的产出。
## =============================================================================

source("cohort.R")
source("models.R")
source("tables.R")   # fmt_cont / fmt_cat / .skew

## ---- ① 可判定 vs 不可判定 ------------------------------------------------------

#' 分组：可判定结局者 / 界标后未复测者 / 界标前已有 AKI 者 / 界标前无法判定者
#'
#' 优先级按方案 3.4 的排除顺序：先"界标前已有 AKI"，再"无法判定"，最后才是
#' "界标后未复测"。同一人可能同时满足多条，只归入最先触发的那一类，各类互斥且
#' 合起来等于未纳入者全体。
.eligibility_class <- function(d) {
  cls <- rep(NA_character_, nrow(d))
  cls[d$evaluable] <- "分析人群"
  rest <- is.na(cls)
  cls[rest & d$aki_pre %in% TRUE]            <- "界标前已有 AKI"
  rest <- is.na(cls)
  cls[rest & d$not_judgeable]                <- "界标前无法判定"
  rest <- is.na(cls)
  cls[rest & !d$has_post_obs]                <- "界标后未复测"
  rest <- is.na(cls)
  cls[rest & d$has_delta_na]                 <- "界标后无肌酐复测以外的其他原因"
  rest <- is.na(cls)
  cls[rest]                                  <- "界标前不满足随访条件"
  factor(cls)
}

#' 可判定者与各类未纳入者的关键变量对比
#'
#' 只取与暴露、结局可检出性直接相关的变量：ΔNa 分布、血钠检测次数与间隔、
#' 界标前最低肌酐、以及反映病情轻重的几项。把全部基线变量都列出来反而会淹没重点。
#' @return data.frame(变量, 分析人群, 各类未纳入者...)
compare_evaluable <- function(d, cls = .eligibility_class(d)) {
  ## 首末血钠间隔：方案 4.3 预先要求保留，用于判断 ΔNa 是否被测量误差主导
  ## （间隔过短时两次检测的差异大半来自分析误差）。由 delta_na 的 t_first/t_last 算。
  if (!"na_gap_h" %in% names(d) && all(c("t_first", "t_last") %in% names(d)))
    d$na_gap_h <- d$t_last - d$t_first
  vars <- list(
    list("delta_na",        "ΔNa（mmol/L）",              "cont", 1),
    list("n_na",            "0–24h 血钠检测次数",          "cont", 0),
    list("na_gap_h",        "首末血钠间隔（h）",            "cont", 1),
    list("creat_baseline",  "界标前最低肌酐（mg/dL）",      "cont", 2),
    list("n_creat_before",  "界标前肌酐检测次数",           "cont", 0),
    list("n_creat_post",    "界标后肌酐检测次数",           "cont", 0),
    list("age",             "年龄（岁）",                  "cont", 0),
    list("gender",          "男性",                        "cat", "M"),
    list("ckd",             "慢性肾病",                    "cat"),
    list("vasoactive_any",  "0–24h 血管活性药",            "cat"),
    list("mech_vent_any",   "0–24h 机械通气",              "cat"),
    list("aki",             "结局窗内 AKI",                "cat")
  )
  lv <- levels(cls)
  lv <- c("分析人群", setdiff(lv, "分析人群"))
  rows <- lapply(vars, function(v) {
    nm <- v[[1]]
    if (!nm %in% names(d)) return(NULL)
    out <- data.frame(变量 = v[[2]], stringsAsFactors = FALSE)
    for (g in lv) {
      x <- d[[nm]][!is.na(cls) & cls == g]
      if (v[[3]] == "cont") {
        xv <- suppressWarnings(as.numeric(x)); xv <- xv[!is.na(xv)]
        st <- if (length(xv) && abs(.skew(xv)) > 1) "median" else "mean"
        val <- fmt_cont(x, if (length(v) >= 4) v[[4]] else 1L, st)
      } else {
        ## 分母固定为该列人数（下方"人数"行给出的就是它），不随变量缺失而变
        val <- fmt_cat(x, if (length(v) >= 4) v[[4]] else NULL,
                       denom = length(x))
      }
      out[[g]] <- val
    }
    out
  })
  out <- do.call(rbind, rows)
  ## 每列人数单独给一行，读者才知道百分比的分母
  hdr <- data.frame(变量 = "人数", stringsAsFactors = FALSE)
  for (g in lv) hdr[[g]] <- as.character(sum(!is.na(cls) & cls == g))
  rbind(hdr, out)
}

#' 复测率与血钠监测强度的分层描述（表 3 用的补充）
#'
#' @return data.frame(血钠检测次数分层, n, 界标后复测比例, AKI 检出率)
retest_by_monitoring <- function(d, breaks = c(-Inf, 2, 3, 5, Inf),
                                 labels = c("2 次", "3 次", "4–5 次", "≥6 次")) {
  g <- cut(d$n_na, breaks = breaks, labels = labels, right = TRUE)
  do.call(rbind, lapply(levels(g), function(k) {
    s <- d[!is.na(g) & g == k, , drop = FALSE]
    data.frame(血钠检测次数 = k, n = nrow(s),
               界标后复测比例 = sprintf("%.1f%%", 100 * mean(s$has_post_obs)),
               分析人群占比 = sprintf("%.1f%%", 100 * mean(s$evaluable)),
               AKI检出率 = if (any(s$evaluable))
                 sprintf("%.1f%%", 100 * mean(s$aki[s$evaluable])) else "—",
               row.names = NULL)
  }))
}

## ---- ② 复测概率权重（IPW） -----------------------------------------------------

#' "除界标后复测之外，其余入选条件都已满足"的分析人群
#'
#' 这是 IPW 的**目标人群**：对这些人而言，"界标后有没有复测肌酐"是唯一挡在
#' 分析人群之外的一道门。因此"复测选择"造成的偏倚只发生在这个人群里，
#' 权重模型也只该在这个人群里拟合。
#'
#' 若把权重模型拟合在全部表型人群（含界标前已患 AKI 者、没有 ΔNa 者）上，
#' 就等于用一批**根本不构成候选**的人去估计复测概率，权重随之失真。
.eligible_except_retest <- function(d) {
  ## 必须要求**已知为否**（%in% FALSE），不能用 `!(x %in% TRUE)`：
  ## 后者把 NA 也判成通过。NA 表示"界标前状态判不出来"，这种人即使补做了复测
  ## 也进不了分析人群，本就不属于目标人群。
  known_false <- function(x) x %in% FALSE
  d$has_delta_na %in% TRUE & known_false(d$aki_pre) & known_false(d$not_judgeable)
}

#' 建"界标后是否复测肌酐"的预测模型，返回逆概率权重
#'
#' 权重口径（方案 5.7 要求写明截断规则）：
#'   - 预测模型用**界标前**协变量（年龄、性别、首次血钠、界标前最低肌酐、CKD、
#'     糖尿病、高血压、冠心病、房颤、瓣膜病、界标前收缩压、血管活性药、机械通气、
#'     界标前血钠检测次数）。**不得用任何界标后变量**，否则权重会把结局信息带进来。
#'   - 模型在 .eligible_except_retest() 的目标人群上拟合。
#'   - 权重 = 1 / P(复测 | 界标前协变量)。
#'   - 截断在 1% 与 99% 分位数（Winsorize）。方案要求报告截断前后权重分布，
#'     故两者都返回。
#' @return list(data, model, 截断点, 分布_截断前, 分布_截断后)
fit_retest_ipw <- function(d, xs = NULL, pool = TRUE) {
  if (isTRUE(pool)) d <- d[.eligible_except_retest(d), , drop = FALSE]
  if (is.null(xs)) {
    xs <- c("age", "gender", "na_first", "creat_baseline", "ckd", "diabetes",
            "hypertension", "chd", "atrial_fib", "valvular", "sbp_min_24h",
            "vasoactive_any", "mech_vent_any", "n_na")
    xs <- .pick(xs, d)
  }
  dd <- d[stats::complete.cases(d[, c("has_post_obs", xs), drop = FALSE]), ,
          drop = FALSE]
  m  <- glm(.fml("has_post_obs", xs), data = dd, family = binomial())
  p  <- stats::predict(m, type = "response")
  w  <- 1 / p
  qs <- stats::quantile(w, c(0.01, 0.99), na.rm = TRUE)
  wq <- pmin(pmax(w, qs[1]), qs[2])
  dd$..w_raw <- w; dd$..w <- wq
  list(data = dd, model = m,
       截断点 = qs,
       分布_截断前 = summary(w), 分布_截断后 = summary(wq))
}

#' 用复测概率权重重跑基础模型（探索性，见方案 5.7）
#'
#' 回归只在**实际被复测者**身上拟合，权重取 1/P(复测 | 界标前协变量)，
#' 于是被复测者按各自被选中的概率还原到目标人群的份额上——
#' 这是标准的"选择偏倚逆概率校正"，估计量对应
#' "若目标人群中所有人都被复测"这一假设人群。
#'
#' ⚠️ 不能把未复测者也放进回归、把他们的结局按 0 处理：未复测者在 4.2.6 里
#' 明确归入"不可判定"，既不计入分母也不计入分子。一旦当作 0 放进来，
#' 他们与已被权重放大过的复测者会各算一遍，等于把同一批人重复计数，
#' 效应量被稀释向无效。这也是本函数与 fit_retest_ipw() 的分工：
#' 前者拟合的是**权重**（要在目标人群全体上拟合），后者拟合的是**加权结局模型**
#' （只在实际有结局的人身上拟合）。
#'
#' 解释边界：校正的是复测选择，不校正治疗混杂；"没有复测不等于没有 AKI"
#' 这一前提仍需在局限性中写明。
ips_retest_model <- function(d, w = NULL) {
  if (is.null(w)) w <- fit_retest_ipw(d)
  dd <- w$data[w$data$has_post_obs %in% TRUE, , drop = FALSE]
  if (nrow(dd) < 50 || sum(dd$aki) < 20)
    stop(sprintf("被复测者的样本过小（n=%d，事件=%d），不足以拟合加权模型",
                 nrow(dd), sum(dd$aki)))
  xs <- c("delta_na", .pick(BASIC_X, dd))
  m  <- suppressWarnings(
    glm(.fml("aki", xs), data = dd, family = binomial(), weights = dd$..w))
  co <- summary(m)$coefficients["delta_na", ]
  data.frame(
    分析 = "复测概率加权（IPW）", n = nrow(dd), 事件数 = sum(dd$aki),
    OR = exp(co[1]), 下限 = exp(co[1] - 1.96 * co[2]),
    上限 = exp(co[1] + 1.96 * co[2]), P = co[4],
    备注 = sprintf("权重截断 %.0f%%–%.0f%%。加权到“除复测外全部条件均满足”的目标人群",
                 100 * 0.01, 100 * 0.99),
    row.names = NULL)
}

#' 打印复测概率模型摘要
print_retest_ipw <- function(w) {
  cat("\n=== 界标后肌酐复测的预测模型（界标前协变量）===\n")
  co <- summary(w$model)$coefficients
  for (i in seq_len(nrow(co)))
    cat(sprintf("  %-16s OR %6.3f  95%%CI %6.3f–%6.3f  P = %s\n",
                rownames(co)[i], exp(co[i, 1]),
                exp(co[i, 1] - 1.96 * co[i, 2]), exp(co[i, 1] + 1.96 * co[i, 2]),
                format.pval(co[i, 4], digits = 2, eps = 1e-4)))
  cat(sprintf("\n  权重（截断前）：最小 %.3f  中位 %.3f  最大 %.3f\n",
              w$分布_截断前[1], w$分布_截断前[3], w$分布_截断前[6]))
  cat(sprintf("  权重（截断后）：最小 %.3f  中位 %.3f  最大 %.3f\n",
              w$分布_截断后[1], w$分布_截断后[3], w$分布_截断后[6]))
  invisible(w)
}
