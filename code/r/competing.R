## =============================================================================
## competing.R —— S11 竞争风险（方案表 10）
##
## 问题：结局窗 (24, 72] 内，患者可能在检出 AKI 之前就转出 ICU 或死亡，此后
## 肌酐不再被测量，结局永远观察不到。普通 logistic 把这些人在分母里按"无事件"
## 计，若转出/死亡与暴露相关，效应量就是错的——而且方向不定。
##
## 两种处理：
##   1. Fine-Gray 次分布风险模型（cmprsk::crr）——协变量对**累积发生函数**的
##      影响，正是临床要解释的量（"发生 AKI 的绝对风险"）。
##   2. Aalen-Johansen 累积发生率（survival::survfit 多状态）——不做模型假设
##      地描述三条曲线：AKI、转出/死亡、仍在观察。
##
## 竞争事件定义：界标后、结局窗内，先于 AKI 发生的 ICU 转出或死亡。
## 注意转出 ICU 后会转入普通病房，肌酐仍可能被查——本方案结局窗止于 72h 且
## 要求界标后仍在 ICU，故把转出视为竞争事件是保守且与队列定义一致的做法；
## 该假设在论文中说明。
## =============================================================================

suppressMessages({ library(survival) })

#' 构造竞争风险数据集
#' @param a 分析人群（须含 intime, outtime, deathtime, t_event, aki）
#' @return data.frame(time, status, ...)  status: 1=AKI, 2=转出/死亡, 0=删失
build_competing <- function(a, horizon_h = HORIZON_H, landmark_h = LANDMARK_H) {
  need <- c("intime", "outtime", "deathtime", "t_event", "aki")
  stopifnot(all(need %in% names(a)))

  t_out <- as.numeric(difftime(a$outtime, a$intime, units = "hours"))
  t_die <- ifelse(is.na(a$deathtime), Inf,
                  as.numeric(difftime(a$deathtime, a$intime, units = "hours")))
  ## 转出与死亡取更早者作竞争事件时刻
  t_cmp <- pmin(t_out, t_die, na.rm = TRUE)

  n <- nrow(a); time <- rep(horizon_h, n); status <- rep(0L, n)

  ## 先填竞争事件：界标后、结局窗内发生转出/死亡者
  is_cmp <- !is.na(t_cmp) & t_cmp > landmark_h & t_cmp <= horizon_h
  time[is_cmp] <- t_cmp[is_cmp]; status[is_cmp] <- 2L
  ## AKI 覆盖竞争事件：AKI 先发生才算数（t_event 严格在 (24,72] 内）
  is_aki <- !is.na(a$aki) & a$aki & !is.na(a$t_event)
  earlier <- is_aki & (!is_cmp | a$t_event < time)
  time[earlier] <- a$t_event[earlier]; status[earlier] <- 1L

  ## 观察时长为 0 或为负者无法进入风险集
  ok <- time > landmark_h
  out <- a[ok, , drop = FALSE]
  out$time <- time[ok]; out$status <- status[ok]
  out
}

#' Fine-Gray 次分布风险模型
#' @return data.frame(项, sHR, 下限, 上限, P, n, 事件数)
fit_finegray <- function(d, xs = c("delta_na", BASIC_X)) {
  xs <- .pick(xs, d)
  d  <- d[stats::complete.cases(d[, xs, drop = FALSE]), , drop = FALSE]
  X  <- stats::model.matrix(stats::as.formula(paste("~", paste(xs, collapse = " + "))),
                            data = d)[, -1, drop = FALSE]
  fit <- cmprsk::crr(d$time, d$status, cov1 = X, failcode = 1, cencode = 0)
  s   <- summary(fit)$coef
  data.frame(项 = rownames(s), sHR = exp(s[, 1]),
             下限 = exp(s[, 1] - 1.96 * s[, 3]),
             上限 = exp(s[, 1] + 1.96 * s[, 3]),
             P = s[, 5], n = nrow(d),
             事件数 = sum(d$status == 1L), row.names = NULL)
}

#' Aalen-Johansen 累积发生率（不做模型假设）
#' @return list(fit, 汇总) —— 汇总为 t=72h 时各状态的累积概率
fit_aalen_johansen <- function(d) {
  ## 多状态：0 观察中 → 1 AKI / 2 转出或死亡
  st <- factor(d$status, levels = c(0, 1, 2),
               labels = c("censor", "AKI", "转出或死亡"))
  fit <- survfit(Surv(d$time, st) ~ 1, data = d)
  ## 72h 处的累积发生率。
  ## **不能写 `状态 = s$states`**：单一时点上 `s$states` 是长度 3 的因子（三个状态
  ## 各出现一次），`s$pstate` 却是 1×3 矩阵，两者按列回收后得到 3 行、每行三个
  ## 概率全都相同——表面上"有表"，实则每个状态都记成了同一组数。按列名取。
  s <- summary(fit, times = HORIZON_H, extend = TRUE)
  p <- s$pstate
  if (is.null(dim(p))) p <- matrix(p, nrow = 1)
  list(fit = fit,
       汇总 = data.frame(状态 = colnames(p),
                         累积发生率 = round(as.numeric(p[1, ]), 4),
                         row.names = NULL))
}
