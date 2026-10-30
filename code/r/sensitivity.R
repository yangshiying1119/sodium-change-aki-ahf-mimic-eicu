## =============================================================================
## sensitivity.R —— 方案表 10 的 S1–S11 敏感性分析
##
## 统一报告口径：每项都给出**基础模型**（BASIC_X）下 ΔNa 的 OR、95% CI、P、
## 分析人群与事件数，便于与主分析逐行对照。方向是否稳定、量级是否漂移，
## 一眼可比。
##
## 关键设计：所有变体都**沿用 cohort.R 的同一套逐层留存与排除逻辑**，只替换
## 暴露构造（na_fun）、结局构造（cr_fun）、提取表（ex）、表型或首次限制。
## 绝不在本文件里另写一份筛选逻辑——否则敏感性分析会悄悄偏离主分析的口径，
## 那正是审稿人会抓的地方。
## =============================================================================

source("cohort.R")
source("models.R")

#' 最早两次血钠之差（S1）
#' 与 delta_na 的区别只在取第 2 条而非末条，用于检验结果是否由全窗轨迹驱动。
delta_na_first2 <- function(dt, t_start = 0, t_end = LANDMARK_H) {
  dt <- merge_same_time(dt)
  d  <- dt[dt$t_hours >= t_start & dt$t_hours <= t_end, , drop = FALSE]
  if (nrow(d) == 0L) return(NULL)
  res <- do.call(rbind, lapply(split(d, d$id), function(x) {
    x <- x[order(x$t_hours), , drop = FALSE]
    if (nrow(x) < 2L) return(NULL)
    data.frame(id = x$id[1], n_na = nrow(x),
               t_first = x$t_hours[1], t_last = x$t_hours[2],
               na_first = x$value[1],  na_last = x$value[2],
               delta_na = x$value[2] - x$value[1], stringsAsFactors = FALSE)
  }))
  if (is.null(res)) return(NULL)
  rownames(res) <- NULL
  res[order(res$id), , drop = FALSE]
}

#' 血糖校正血钠（S5）
#'
#' 公式用临床通用的 Katz 校正（Hillier 简化式）：血钠每被 100 mg/dL 以上的血糖
#' 抬高 100 mg/dL 约 2.4 mmol/L，即 Na_c = Na + 0.024 × (Glu − 100)，血糖 ≤100 不校正。
#' 高血糖是急性心衰的常见合并状态，若不校正，"血钠变化"里会混进血糖波动带来的
#' 稀释性低钠成分。
#'
#' 配对规则：每条血钠取**同一次住院内、±1 小时内最接近**的一条血糖，不做插值、
#' 不跨时点搬运。1 小时与血钠/肌酐的"同时间"容差（15 分钟）不同——血糖与电解质
#' 常在同一管血但录入时刻略有先后，窗口过窄会大量配不上；而 1 小时内血糖的变化
#' 幅度远小于校正项本身。
#' @return 与输入同构的长表 (id, t_hours, value)，value 为校正血钠
correct_sodium <- function(sodium, glucose, tol_h = 1) {
  sod <- merge_same_time(sodium)
  glu <- merge_same_time(glucose)
  if (is.null(glu) || !nrow(glu)) return(sod)
  keep_ids <- unique(sod$id)
  glu <- glu[glu$id %in% keep_ids, , drop = FALSE]
  if (!nrow(glu)) return(sod)

  gl <- split(glu, glu$id)
  out <- lapply(split(sod, sod$id), function(x) {
    g <- gl[[x$id[1]]]
    if (is.null(g)) return(x)
    gt <- g$t_hours; gv <- g$value
    for (i in seq_len(nrow(x))) {
      d <- abs(gt - x$t_hours[i])
      j <- which(d <= tol_h)
      if (length(j)) {
        j <- j[which.min(d[j])]
        if (gv[j] > 100) x$value[i] <- x$value[i] + 0.024 * (gv[j] - 100)
      }
    }
    x
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out[order(out$id, out$t_hours), , drop = FALSE]
}

#' 跑一份变体并返回基础模型下 ΔNa 的效应量
#' @param label   显示名
#' @param ex_use  使用的提取结果（默认与主分析同一份；S4 换成血气钠）
#' @param extra_x 追加到 BASIC_X 之后的协变量（如 n_na）
#' @param drop    并上协变量**之后**追加的过滤函数，返回逻辑向量
#' @param ...     其余参数透传给 run_mimic_cohort
#' @return 一行 data.frame；失败时返回带 NA 与原因的行，不中断整批
one_sens <- function(ex, label, ex_use = ex, extra_x = character(0),
                     drop = NULL, ...) {
  out <- try({
    r <- run_mimic_cohort(ex = ex_use, ...)
    a <- attach_covariates(r$analysis, ex$covariates)
    ## 与主分析一样取完整病例，使每一项变体内部也只有一个分母；
    ## 各变体之间 n 可以不同（人群本就变了），但变体内不得静默按行删除。
    a <- a[model_cases(a), , drop = FALSE]
    if (!is.null(drop)) a <- a[drop(a), , drop = FALSE]
    if (nrow(a) < 50 || sum(a$aki) < 20)
      stop(sprintf("样本过小（n=%d，事件=%d），不足以拟合基础模型", nrow(a), sum(a$aki)))
    xs <- c("delta_na", .pick(c(BASIC_X, extra_x), a))
    m  <- glm(.fml("aki", xs), data = a, family = binomial())
    co <- summary(m)$coefficients["delta_na", ]
    data.frame(分析 = label, n = nrow(a), 事件数 = sum(a$aki),
               OR = exp(co[1]), 下限 = exp(co[1] - 1.96 * co[2]),
               上限 = exp(co[1] + 1.96 * co[2]), P = co[4], 备注 = "",
               row.names = NULL)
  }, silent = TRUE)
  if (inherits(out, "try-error")) {
    msg <- attr(out, "condition")
    data.frame(分析 = label, n = NA_integer_, 事件数 = NA_integer_,
               OR = NA_real_, 下限 = NA_real_, 上限 = NA_real_, P = NA_real_,
               备注 = if (is.null(msg)) "失败" else conditionMessage(msg),
               row.names = NULL)
  } else out
}

#' 多重插补敏感性分析：把"统一完整病例"丢掉的缺失协变量补回来
#'
#' 主分析为让三档模型同分母，剔除了协变量缺失者（本例 44 例 / 1.6%）。剔除本身
#' 是可辩护的，但审稿人会问"这 44 人是不是系统性地不同"。本项在**未剔除**的分析
#' 人群上做链式方程多重插补（mice），把基础与扩展模型各重跑一遍并按 Rubin 规则
#' 合并。若合并后的 OR 与主分析几乎一致，缺失就不构成威胁。
#'
#' 插补模型的变量集固定为 MODEL_VARS（与建模变量一致），不额外引入结局之外的
#' 信息；seed 固定，结果可复现。**注意插补是为了补齐协变量，不是补齐结局**——
#' `aki` 与 `delta_na` 在主分析人群里本就完整无缺，此处一并纳入插补矩阵只作为
#' 插补模型的预测变量。
#'
#' @param d 分析人群（covariate-attached，**未**按 model_cases 剔除）
#' @return 若干行 data.frame，列与 run_sensitivity 的主表一致
run_mi <- function(d, m = 20L) {
  if (!requireNamespace("mice", quietly = TRUE))
    return(data.frame(分析 = "附加 多重插补（mice）", n = NA_integer_,
                      事件数 = NA_integer_, OR = NA_real_, 下限 = NA_real_,
                      上限 = NA_real_, P = NA_real_,
                      备注 = "未安装 mice 包", row.names = NULL))
  suppressMessages(library(mice))   # 需要 attach 才能让 with() 派发到 with.mids
  vars <- .pick(MODEL_VARS, d)
  dd   <- d[stats::complete.cases(d[, c("aki", "delta_na"), drop = FALSE]),
            vars, drop = FALSE]
  out <- try({
    imp <- mice::mice(dd, m = m, seed = 20261001, printFlag = FALSE)
    rows <- lapply(list(基础 = .pick(BASIC_X, dd), 扩展 = .pick(EXTEND_X, dd)),
                   function(xs) {
      ## 公式必须在 with() 的表达式**内部**构造：with.mids 是把插补数据放进求值
      ## 环境，而不是传给 glm 的 data=，因此 model.frame 只会去 environment(formula)
      ## 找结局变量。若在函数外层先建好 fml，它的环境是本函数帧（没有 aki），
      ## 会直接报 "object 'aki' not found"。
      fml <- paste("aki ~", paste(c("delta_na", xs), collapse = " + "))
      fit <- with(imp, glm(as.formula(fml), family = binomial()))
      s   <- summary(mice::pool(fit), conf.int = TRUE)
      r   <- s[s$term == "delta_na", ]
      data.frame(OR = exp(r$estimate), 下限 = exp(r$conf.low),
                 上限 = exp(r$conf.high), P = r$p.value)
    })
    do.call(rbind, lapply(seq_along(rows), function(i)
      data.frame(分析 = sprintf("附加 多重插补（%s模型，m=%d）",
                                names(rows)[i], m),
                 n = nrow(dd), 事件数 = sum(dd$aki),
                 rows[[i]], 备注 = "Rubin 规则合并", row.names = NULL)))
  }, silent = TRUE)
  if (inherits(out, "try-error"))
    data.frame(分析 = "附加 多重插补（mice）", n = NA_integer_, 事件数 = NA_integer_,
               OR = NA_real_, 下限 = NA_real_, 上限 = NA_real_, P = NA_real_,
               备注 = conditionMessage(attr(out, "condition")), row.names = NULL)
  else out
}

#' S1–S10 + 附加分析
#' @param main_row 主分析的结果行，作为对照放在第一行
#' @param d_main   分析人群（covariate-attached）。给了才跑多重插补，否则跳过。
#' @param cohort   表型口径字母。**必须与主分析所用口径一致**：每一项变体都会
#'   重跑一遍队列组装，而 cohort 的默认值是 "B"，若不显式传入，在 eICU 侧
#'   （重复队列用表型 C）所有变体都会悄悄退回表型 B 的人群，
#'   得到一族"n≈187"的、与主分析并非同一人群的结果——表面上每一行都算出来了，
#'   口径却整个错了。S6 例外：它的定义就是换到严格表型，故显式覆盖。
run_sensitivity <- function(ex, main_row = NULL, d_main = NULL, cohort = "B") {
  rrt_post <- unique(ex$rrt$id[ex$rrt$t_hours > LANDMARK_H &
                               ex$rrt$t_hours <= HORIZON_H])
  ## S4：只把血钠长表换成血气钠，判定层与其余一切不变
  ex_bg <- ex
  if (!is.null(ex$sodium_blood) && nrow(ex$sodium_blood) > 0)
    ex_bg$sodium <- ex$sodium_blood

  rows <- list()
  if (!is.null(main_row)) rows[[length(rows) + 1L]] <- main_row

  ## 每一项都带上主分析的表型口径；调用方显式给了 cohort 就以调用方为准
  ## （S6 正是靠这一条切到严格表型），否则补上本函数的 cohort。
  add <- function(...) {
    dots <- list(...)
    if (is.null(dots$cohort)) dots$cohort <- cohort
    rows[[length(rows) + 1L]] <<- do.call(one_sens, c(list(ex = ex), dots))
  }

  add("S1 最早两次血钠之差",
      na_fun = function(dt) delta_na_first2(dt))
  add("S2 暴露窗改为 0–18h",
      na_fun = function(dt) delta_na(dt, t_end = 18))
  ## S3：剔除 |ΔNa| > 12。24h 内血钠变化超出临床可实现范围者多为录入/单位异常；
  ## 方案 4.1.3 要求回看原值、不采用依赖数据分布的截尾，故此处用固定阈值。
  add("S3 剔除 |ΔNa| > 12 mmol/L",
      drop = function(a) !is.na(a$delta_na) & abs(a$delta_na) <= 12)
  if (!is.null(ex$sodium_blood) && nrow(ex$sodium_blood) > 0) {
    add("S4 血气钠（50824）替换生化钠", ex_use = ex_bg)
    ## S4b：把主分析原样限制到"有 ≥2 次血气钠"的那批人，**仍用生化钠**。
    ## 这一行是 S4 归因的关键：血气一般只在较重、需要频繁查血气的人身上做，
    ## 若这批人本身的 OR 就已接近 1，S4 的零结果来自**选择**而非标本类型；
    ## 只有 S4b 仍接近主分析、独独 S4 掉到 1，才说明是血气钠的测量误差。
    bg <- ex$sodium_blood
    n_bg <- table(bg$id[bg$t_hours >= 0 & bg$t_hours <= LANDMARK_H])
    bg_ids <- names(n_bg)[n_bg >= 2L]
    add("S4b 主分析限制到有 ≥2 次血气钠者（仍用生化钠）",
        drop = function(a) a$id %in% bg_ids)
  }
  ## S5：血糖校正血钠。需要逐时刻血糖序列；未提取到就写明原因，不静默略过
  if (!is.null(ex$glucose) && nrow(ex$glucose) > 0) {
    add("S5 血糖校正血钠",
        na_fun = function(dt) delta_na(correct_sodium(dt, ex$glucose)))
  } else {
    rows[[length(rows) + 1L]] <- data.frame(
      分析 = "S5 血糖校正血钠", n = NA_integer_, 事件数 = NA_integer_,
      OR = NA_real_, 下限 = NA_real_, 上限 = NA_real_, P = NA_real_,
      备注 = "提取结果中无血糖长表", row.names = NULL)
  }
  add("S6 表型定义 A（严格）", cohort = "A")
  ## S7：只用 48h 内 ≥0.3 mg/dL 这一条
  add("S7 仅保留 48h ≥0.3 mg/dL 标准",
      cr_fun = function(dt) kdigo_first_event(dt, pre_from = 0, rrt = rrt_post,
                                              use_ratio = FALSE))
  ## S9：排除 0–24h 接受高渗盐水或已行 RRT 者（血钠被直接干预）
  add("S9 排除高渗盐水或 RRT 者",
      drop = function(a) a$hypertonic_saline == 0 & !a$rrt_pre & !a$rrt_post)
  ## S10：不施加首次住院限制
  add("S10 不限制首次住院", first_only = FALSE)
  ## S12：首末两次血钠至少相隔 6 h。
  ## 暴露是"末次减首次"，方案只要求 0–24h 内 ≥2 次，未规定最小间隔；若两次
  ## 抽血相隔一两小时，ΔNa 里很大一部分是分析前误差与生理波动，而不是 24 h
  ## 内的真实变化。6 h 是预先定下的阈值（覆盖半个班次，足以让变化有意义），
  ## 不按结果调整。
  add("S12 首末两次血钠间隔 ≥6 h",
      drop = function(a) !is.na(a$t_last - a$t_first) & (a$t_last - a$t_first) >= 6)
  ## 附加 1：调整血钠检测次数（监测强度的代理变量；审稿人必问）
  add("附加 调整血钠检测次数", extra_x = "n_na")
  ## 附加 2：剔除界标后 24h 内即转出/死亡者，检验结局 ascertainment 的影响
  add("附加 仅保留界标后留观 ≥48h 者",
      drop = function(a) !is.na(a$los_hours) & a$los_hours >= HORIZON_H)
  ## 附加 3：多重插补。与上面各项不同，它不用 one_sens——它要的是**未剔除**的
  ## 分析人群，而 one_sens 走的是同一套完整病例流水线。
  if (!is.null(d_main)) rows[[length(rows) + 1L]] <- run_mi(d_main)

  do.call(rbind, rows)
}

#' S8：结局改用"肌酐标准 ∪ 尿量标准"的复合 KDIGO 口径
#'
#' 主分析只用肌酐。肌酐是滞后指标且受肌肉量影响，单用肌酐会漏掉少尿型 AKI——
#' 这是审稿人对 AKI 研究最常见的质询。本项把尿量标准并列进来，界标前的"已满足
#' KDIGO"与界标后结局同步替换，其余环节（暴露、筛选、协变量）一律不动。
#'
#' @return list(行 = 森林图用的一行, 队列 = 复合口径的组装结果, 交叉 = 两套标准的四格表)
run_s8_uo <- function(ex, cohort = "B") {
  if (is.null(ex$urine) || nrow(ex$urine) == 0L)
    stop("提取结果中没有尿量长表（urine_ts），S8 无法进行")
  if (is.null(ex$weight) || !length(ex$weight))
    stop("提取结果中没有体重，尿量速率无法按 kg 标化")

  ## cohort 必须与主分析一致：run_mimic_cohort 的默认口径是 B，
  ## 在 eICU 侧（重复队列用表型 C）漏传会静默换掉整张 S8 表的人群。
  r <- run_mimic_cohort(cohort = cohort, ex = ex,
                        uo = list(urine = ex$urine, weight = ex$weight))
  a <- attach_covariates(r$analysis, ex$covariates)
  a <- a[model_cases(a), , drop = FALSE]   # 与主分析同一分母口径
  xs <- c("delta_na", .pick(BASIC_X, a))
  m  <- glm(.fml("aki", xs), data = a, family = binomial())
  co <- summary(m)$coefficients["delta_na", ]

  ## 两套标准的交叉表：谁贡献了额外事件，是判断复合口径有没有实质影响的关键。
  ## 落成 data.frame 而不是 table——write.csv 对 table 会走 matrix 分支，
  ## 表头变成 "FALSE"/"TRUE"，出到表件里读者看不懂。
  cross <- data.frame(
    肌酐标准 = c("检出", "未检出", "合计"),
    尿量标准检出   = c(sum(a$aki_creat & a$uo_aki),
                       sum(!a$aki_creat & a$uo_aki),
                       sum(a$uo_aki)),
    尿量标准未检出 = c(sum(a$aki_creat & !a$uo_aki),
                       sum(!a$aki_creat & !a$uo_aki),
                       sum(!a$uo_aki)),
    check.names = FALSE, row.names = NULL)
  cross[["合计"]] <- cross[["尿量标准检出"]] + cross[["尿量标准未检出"]]

  list(行 = data.frame(
         分析 = "S8 复合口径（肌酐 ∪ 尿量）", n = nrow(a), 事件数 = sum(a$aki),
         OR = exp(co[1]), 下限 = exp(co[1] - 1.96 * co[2]),
         上限 = exp(co[1] + 1.96 * co[2]), P = co[4],
         备注 = sprintf("其中仅尿量标准检出 %d 例", sum(!a$aki_creat & a$uo_aki)),
         row.names = NULL),
       队列 = r, 交叉 = cross,
       仅尿量 = a[!a$aki_creat & a$uo_aki, , drop = FALSE])
}

#' 打印敏感性分析表
print_sensitivity <- function(tab) {
  cat(sprintf("\n%-34s %6s %6s %-22s %s\n", "分析", "n", "事件", "OR（95% CI）", "P"))
  for (i in seq_len(nrow(tab))) {
    ci <- if (is.na(tab$OR[i])) tab$备注[i]
          else sprintf("%.3f（%.3f–%.3f）", tab$OR[i], tab$下限[i], tab$上限[i])
    cat(sprintf("%-34s %6s %6s %-22s %s\n", tab$分析[i],
                ifelse(is.na(tab$n[i]), "—", tab$n[i]),
                ifelse(is.na(tab$事件数[i]), "—", tab$事件数[i]), ci,
                ifelse(is.na(tab$P[i]), "—",
                       format.pval(tab$P[i], digits = 2, eps = 1e-4))))
  }
  invisible(tab)
}
