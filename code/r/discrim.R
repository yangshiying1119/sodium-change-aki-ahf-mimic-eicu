## 判别力、校准与临床量级（方案 5.3 三处要求的落地）
##
## 方案 5.3 要求除了给出 OR，还要给出：按每 3 mmol/L 变化的 OR、按代表性初始血钠
## 与肾功能标准化的预测概率、以及判别力（AUC / c 统计量）与校准（校准曲线或
## Hosmer–Lemeshow）。这三项在原分析里没有落盘，本文件补齐。
##
## 为什么必须补：OR 0.93 这个量级，读者无法自行判断它在临床上算大还是算小。
## 判别力（加了 ΔNa 之后 c 统计量提高多少）与标准化预测概率，
## 是把"统计显著"翻译成"临床上有多少用"的唯一手段，也是审稿人一定会问的。

#' 每 3 mmol/L 变化的 OR
#'
#' 模型里 ΔNa 是连续线性项，每 1 mmol/L 的 OR 为 exp(b)。方案 5.3 要求按
#' 临床上更常见的 3 mmol/L 幅度报告，由 exp(3b) 直接换算，不重新拟合。
or_per_unit <- function(m, term = "delta_na", step = 3) {
  co <- summary(m)$coefficients
  if (!term %in% rownames(co)) return(NULL)
  b  <- co[term, 1]; se <- co[term, 2]
  data.frame(步长 = step,
             OR = exp(step * b),
             下限 = exp(step * b - 1.96 * step * se),
             上限 = exp(step * b + 1.96 * step * se),
             P = co[term, 4])
}

#' Hosmer–Lemeshow 拟合优度检验（十分位）
hl_test <- function(y, p, g = 10) {
  br <- unique(quantile(p, probs = seq(0, 1, length.out = g + 1), na.rm = TRUE))
  if (length(br) < 3) return(NULL)
  grp <- cut(p, breaks = br, include.lowest = TRUE)
  obs <- tapply(y, grp, sum); exp_ <- tapply(p, grp, sum)
  n   <- tapply(y, grp, length)
  ok  <- !is.na(obs) & !is.na(exp_)
  obs <- obs[ok]; exp_ <- exp_[ok]; n <- n[ok]
  ## 期望频数过小的格子会让卡方近似失效，直接告知调用方，不静默出数
  chi <- sum((obs - exp_)^2 / (exp_ * (1 - exp_ / n)))
  df  <- length(obs) - 2
  list(卡方 = chi, 自由度 = df,
       P = pchisq(chi, df = df, lower.tail = FALSE),
       格子最小期望 = min(exp_))
}

#' 校准斜率与截距（logistic 重校准）
#'
#' **注意：在同一批数据上把线性预测子当唯一自变量重拟合，斜率恒等于 1、截距恒
#' 等于 0**——这是原模型得分方程的直接推论（Σ(y−p)=0 与 Σlp(y−p)=0），
#' 不是拟合良好的证据。所以这里一并给出 bootstrap 乐观校正后的斜率，
#' 表里只报校正值；表观值单列出来仅用于说明两者差多少（差得越多说明越乐观）。
cal_slope_intercept <- function(y, lp) {
  m <- glm(y ~ lp, family = binomial())
  co <- summary(m)$coefficients
  list(表观斜率 = co["lp", 1], 表观截距 = co["(Intercept)", 1])
}

#' Bootstrap 乐观校正后的判别力与校准斜率
#'
#' 方案 5.3 要的"判别力（AUC 或 c 统计量）与校准"若只给表观值，等于把过拟合
#' 当成性能。这里按标准的 optimism bootstrap 做校正：每次重抽样都重新拟合，
#' 用"在重抽样样本上的表现减去在原样本上的表现"作为该次的乐观量，
#' 再取平均从表观值里扣掉。同时给出校正值的 bootstrap 区间。
#'
#' 不用 rms::validate 的原因：它只返回汇总的五个数，不保留每次重抽样的结果，
#' 拿不到区间；而且它会连带加载 ggplot2，与 R 4.1.2 下的版本有告警。
#'
#' @return list(c表观, c校正, c校正下限, c校正上限, 斜率校正, B)
boot_validate <- function(d, adj, B = 200, seed = 20260930) {
  xs <- c("delta_na", adj)
  f  <- stats::as.formula(paste("aki ~", paste(xs, collapse = " + ")),
                          env = environment())
  y  <- as.integer(d$aki)
  c_of <- function(m, dd, yy)
    as.numeric(pROC::auc(pROC::roc(yy, stats::predict(m, newdata = dd,
                                                      type = "response"),
                                   quiet = TRUE)))
  slope_of <- function(m, dd) {
    lp <- stats::predict(m, newdata = dd, type = "link")
    co <- summary(glm(as.integer(dd$aki) ~ lp, family = binomial()))$coefficients
    co["lp", 1]
  }
  fit0  <- glm(f, data = d, family = binomial())
  c_app <- c_of(fit0, d, y)
  s_app <- slope_of(fit0, d)
  set.seed(seed)                         # 校正量本身随机，固定种子以便复现
  n <- nrow(d); dc <- ds <- numeric(B)
  for (b in seq_len(B)) {
    i  <- sample.int(n, n, replace = TRUE)
    db <- d[i, , drop = FALSE]
    mb <- tryCatch(glm(f, data = db, family = binomial()),
                   error = function(e) NULL)
    if (is.null(mb)) { dc[b] <- NA; ds[b] <- NA; next }
    dc[b] <- c_of(mb, db, as.integer(db$aki)) - c_of(mb, d, y)
    ds[b] <- slope_of(mb, db) - slope_of(mb, d)
  }
  ck <- c_app - dc; sk <- s_app - ds
  list(c表观 = c_app, 乐观量 = mean(dc, na.rm = TRUE),
       c校正 = c_app - mean(dc, na.rm = TRUE),
       c校正下限 = stats::quantile(ck, 0.025, na.rm = TRUE, names = FALSE),
       c校正上限 = stats::quantile(ck, 0.975, na.rm = TRUE, names = FALSE),
       斜率校正 = s_app - mean(ds, na.rm = TRUE),
       B = sum(!is.na(dc)))
}

#' 校准曲线用的分位分组：每组的平均预测概率与实测发生率
cal_curve <- function(y, p, g = 10) {
  br  <- unique(stats::quantile(p, probs = seq(0, 1, length.out = g + 1)))
  grp <- cut(p, breaks = br, include.lowest = TRUE)
  agg <- do.call(rbind, lapply(levels(grp), function(k) {
    i <- which(grp == k)
    if (!length(i)) return(NULL)
    data.frame(组 = k, n = length(i),
               平均预测 = mean(p[i]), 实测发生率 = mean(y[i]))
  }))
  agg
}

#' 判别力：基础模型 vs 去掉 ΔNa 的模型
#'
#' ΔAUC 的区间与检验走 pROC 的 DeLong 法。**只报加了 ΔNa 之后的 AUC，
#' 不报增量，等于回避了"这个指标到底有没有增量价值"这个问题**，
#' 故两个模型一并报，并给出似然比检验。
discrim_auc <- function(y, p_full, p_null) {
  r1 <- pROC::roc(y, p_full, quiet = TRUE)
  r0 <- pROC::roc(y, p_null, quiet = TRUE)
  ci <- as.numeric(pROC::ci.auc(r1, method = "delong"))
  tt <- tryCatch(pROC::roc.test(r1, r0, method = "delong"),
                 error = function(e) NULL)
  list(AUC = as.numeric(pROC::auc(r1)),
       AUC下限 = ci[1], AUC上限 = ci[3],
       AUC无暴露 = as.numeric(pROC::auc(r0)),
       差值 = as.numeric(pROC::auc(r1)) - as.numeric(pROC::auc(r0)),
       差值P = if (is.null(tt)) NA_real_ else tt$p.value)
}

#' 标准化预测概率
#'
#' 方案 5.3 要求"按代表性初始血钠与肾功能标准化"。做法是在**该层的实际人群**上
#' 把 ΔNa 设为给定值、其余协变量保持各自实测值，再对全层取平均——
#' 即对协变量分布作标准化，而不是把协变量都钉在均值上。
#' 后者在含分类变量与非线性项时会给出人群里并不存在的"平均人"，
#' 预测概率随之失去意义。
#'
#' @param d 已切好的完整病例集（含协变量与结局）
#' @param m 基础模型
#' @param deltas 要评估的 ΔNa 值
#' @return data.frame(ΔNa, 各层的标准化预测概率)
std_pred_prob <- function(d, m, deltas = c(-5, -3, 0, 3, 5)) {
  na0  <- d$na_first                     # 首次血钠，见 models.R 的 BASIC_X
  if (is.null(na0)) return(NULL)
  cr0  <- d$creat_baseline
  cr_med <- stats::median(cr0, na.rm = TRUE)
  lay <- list(
    "首次血钠 <135 且基线肌酐 ≤中位" = na0 <  135 & cr0 <= cr_med,
    "首次血钠 <135 且基线肌酐 >中位" = na0 <  135 & cr0 >  cr_med,
    "首次血钠 ≥135 且基线肌酐 ≤中位" = na0 >= 135 & cr0 <= cr_med,
    "首次血钠 ≥135 且基线肌酐 >中位" = na0 >= 135 & cr0 >  cr_med)
  out <- data.frame(ΔNa = deltas, check.names = FALSE)
  for (k in names(lay)) {
    idx <- which(lay[[k]])
    out[[k]] <- vapply(deltas, function(v) {
      dd <- d[idx, , drop = FALSE]
      dd$delta_na <- v
      mean(stats::predict(m, newdata = dd, type = "response"))
    }, numeric(1))
  }
  attr(out, "层内人数") <- vapply(lay, sum, integer(1))
  out
}

#' 一次算齐上述四项
#'
#' @param d 完整病例集
#' @return list(每3mmol, HL, 校准, 判别, 标准化概率, 层内人数, 基础模型系数)
discrim_all <- function(d, adj = .pick(BASIC_X, d)) {
  y <- as.integer(d$aki)
  m  <- glm(.fml("aki", c("delta_na", adj)), data = d, family = binomial())
  m0 <- glm(.fml("aki", adj),                data = d, family = binomial())
  pf <- stats::fitted(m); p0 <- stats::fitted(m0)
  lp <- stats::predict(m, type = "link")
  ll1 <- as.numeric(stats::logLik(m)); ll0 <- as.numeric(stats::logLik(m0))
  list(
    n = nobs(m), 事件数 = sum(y),
    每3mmol = or_per_unit(m, "delta_na", 3),
    LRT = list(卡方 = 2 * (ll1 - ll0), 自由度 = 1,
               P = pchisq(2 * (ll1 - ll0), df = 1, lower.tail = FALSE)),
    HL = hl_test(y, pf),
    校准 = cal_slope_intercept(y, lp),
    校准曲线 = cal_curve(y, pf),
    判别 = discrim_auc(y, pf, p0),
    校正 = boot_validate(d, adj),
    标准化概率 = std_pred_prob(d, m),
    层内人数 = attr(std_pred_prob(d, m), "层内人数"),
    模型 = m, 数据 = d)
}
