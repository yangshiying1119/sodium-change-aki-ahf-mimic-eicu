## =============================================================================
## common.R —— 两库共享的判定层
##
## 设计要点：KDIGO 判定与 ΔNa 计算只在本文件实现一次。MIMIC-IV 与 eICU-CRD
## 各自把数据归一化成同一张长表后调用这里的函数，"两库口径一致"由代码结构
## 保证，不靠人工对齐。两库的差异只允许留在各自的提取适配层（extract.R）。
##
## 长表约定（三个库一律遵守）：
##   id      字符  患者-停留唯一键（MIMIC 用 stay_id，eICU 用 patientunitstayid）
##   t_hours 数值  相对 t0（ICU 入科时刻）的小时数；**可为负**，用于容纳院前实测值
##   value   数值  指标值
##
## 本文件不连接任何数据库，纯函数，可独立单元测试（见 test_common.R）。
## =============================================================================

## ---- 0. 语言层（OUT_LANG=en 时表与图输出英文，zh 时原样） --------------------
## 只换标签不换数：数值仍由各模块算一次，翻译发生在落盘前。
.here_i18n <- local({
  a <- commandArgs(trailingOnly = FALSE)
  s <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(s) && nzchar(s)) dirname(normalizePath(s)) else "."
})
source(file.path(.here_i18n, "i18n.R"))

## ---- 1. 清洗边界（方案 4.1.2 / 4.2.2 预先规定，不得依数据分布调整） ----------

SODIUM_RANGE <- c(90, 200)    # mmol/L（= mEq/L）；方案：<90 或 >200 视为录入错误
CREAT_RANGE  <- c(0.1, 50)    # mg/dL；方案：剔除 0.0 与 808 这类录入错误

SAME_TIME_TOL_H <- 0.25       # 同时间重复检测容差：15 分钟 = 0.25 小时

## 界标与观察窗（小时）
LANDMARK_H <- 24              # 暴露窗 [0, 24]；界标 = 24
HORIZON_H  <- 72              # 结局窗 (24, 72]
PRE_H      <- 168             # 比较基线可回溯至 t0 - 7 天（院前实测值）

clean_sodium <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(!is.na(x) & x >= SODIUM_RANGE[1] & x <= SODIUM_RANGE[2], x, NA_real_)
}

clean_creat <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(!is.na(x) & x >= CREAT_RANGE[1] & x <= CREAT_RANGE[2], x, NA_real_)
}

## ---- 2. 同时间重复检测合并（方案 4.1.2） -------------------------------------
## 规则预先固定：同一时刻（或相差 ≤15 分钟）多条时取中位数；若为偶数条则取
## **较低者**（保守方向）。注意官方 mimic-code 此处用 AVG，本方案按方案原文
## 用中位数/较低者，刻意不跟官方，以免与预注册口径不符。

.merge_one_id <- function(t, v, tol_h = SAME_TIME_TOL_H) {
  o <- order(t, v)
  t <- t[o]; v <- v[o]
  ## 相邻间隔超过容差即开新簇
  grp <- cumsum(c(TRUE, diff(t) > tol_h))
  t_rep <- tapply(t, grp, min)                    # 簇的代表时刻 = 簇内最早时刻
  v_rep <- tapply(seq_along(v), grp, function(i) {
    s <- sort(v[i]); n <- length(s)
    if (n %% 2L == 1L) s[(n + 1L) %/% 2L] else s[n %/% 2L]   # 偶数条取较低者
  })
  data.frame(t_hours = as.numeric(t_rep), value = as.numeric(v_rep))
}

#' @param dt data.frame(id, t_hours, value)，value 须已清洗（NA 行应先剔除）
merge_same_time <- function(dt) {
  dt <- dt[!is.na(dt$value), , drop = FALSE]
  if (nrow(dt) == 0L) return(dt)
  parts <- lapply(split(dt, dt$id), function(d) .merge_one_id(d$t_hours, d$value))
  out <- do.call(rbind, lapply(names(parts), function(k)
    data.frame(id = k, parts[[k]], stringsAsFactors = FALSE)))
  rownames(out) <- NULL
  out[order(out$id, out$t_hours), , drop = FALSE]
}

## ---- 3. 暴露：ΔNa = 0–24h 末次 − 首次（方案 4.3） ---------------------------

#' @return data.frame(id, n_na, t_first, t_last, na_first, na_last, delta_na)
#'         0–24h 内血钠不足 2 次者不出现在结果中（该层即"可算 ΔNa"人群）
delta_na <- function(dt, t_start = 0, t_end = LANDMARK_H) {
  dt <- merge_same_time(dt)
  in_win <- dt$t_hours >= t_start & dt$t_hours <= t_end
  d <- dt[in_win, , drop = FALSE]
  if (nrow(d) == 0L) return(NULL)

  res <- do.call(rbind, lapply(split(d, d$id), function(x) {
    x <- x[order(x$t_hours), , drop = FALSE]
    if (nrow(x) < 2L) return(NULL)
    data.frame(
      id       = x$id[1],
      n_na     = nrow(x),
      t_first  = x$t_hours[1],
      t_last   = x$t_hours[nrow(x)],
      na_first = x$value[1],
      na_last  = x$value[nrow(x)],
      delta_na = x$value[nrow(x)] - x$value[1],
      stringsAsFactors = FALSE
    )
  }))
  if (is.null(res)) return(NULL)
  rownames(res) <- NULL
  res[order(res$id), , drop = FALSE]
}

## ---- 4. KDIGO 肌酐分期 --------------------------------------------------------
## 严格只用实测值。比较值一律取自**严格早于**当前时点的检测（t' < t），
## 绝不使用事件之后的更低值——这是本方案相对宽松做法的核心区别，也是
## 文献中 AKI 发生率报得异常高的常见成因。禁止 MDRD / eGFR=75 假定基线。
##
## 与官方 mimic-code/kdigo_stages.sql 一致：3 期 = 3× 基线；或 ≥4.0 且伴急性
## 升高；2 期 = 2×；1 期 = 48h 内升 ≥0.3 或 1.5×。尿量标准不在本函数内。

## KDIGO 用的是 ≥（"≥0.3 mg/dL"、"≥1.5 倍"），等号必须计入。但 `1.60 + 0.3` 在双
## 精度下得到 1.90000000000000013，而 1.90 的字面值是 1.89999999999999991——直接
## 比较会把"恰好升高 0.30"判成未达标（实测 38777702 就中招：t=41.7 的 0.30 升高被
## 漏检，事件被推迟到 t=54.4）。加 1e-9 容差消掉纯浮点误差；肌酐报告精度为
## 0.01 mg/dL，容差比它小 7 个数量级，不会制造假阳性。
KDIGO_EPS <- 1e-9

#' @param use_48h   是否启用"48 小时内升高 ≥0.3 mg/dL"这一条
#' @param use_ratio 是否启用"升至基线 1.5/2/3 倍"这一组（7 天窗）
#'   两个开关只在 S7 敏感性分析里关掉其一；主分析一律 Both TRUE。
.kdigo_stage_scalar <- function(creat, low48, low7, eps = KDIGO_EPS,
                                use_48h = TRUE, use_ratio = TRUE) {
  if (is.na(creat)) return(NA_integer_)
  ge <- function(a, b) !is.na(b) && a >= b - eps
  a_abs <- use_48h   && ge(creat, low48 + 0.3)   # 48h 内绝对升高
  r1    <- use_ratio && ge(creat, 1.5 * low7)    # 1.5 倍基线
  r2    <- use_ratio && ge(creat, 2.0 * low7)
  r3    <- use_ratio && ge(creat, 3.0 * low7)
  ## 3 期：升至基线 7 天最低值的 3 倍
  if (r3) return(3L)
  ## 3 期：≥4.0 mg/dL，且此前确有急性升高
  if (creat >= 4.0 && (a_abs || r1)) return(3L)
  if (r2) return(2L)
  if (a_abs || r1) return(1L)
  0L
}

#' 为单个 id 的肌酐序列计算每个时点的分期
#' @param t 数值向量，小时（可含负值 = 院前）
#' @param v 数值向量，mg/dL
#' @return data.frame(t_hours, creat, low48, low7, stage)
.kdigo_series <- function(t, v, use_48h = TRUE, use_ratio = TRUE) {
  o <- order(t); t <- t[o]; v <- v[o]
  n <- length(t)
  low48 <- rep(NA_real_, n); low7 <- rep(NA_real_, n); stage <- rep(NA_integer_, n)
  for (i in seq_len(n)) {
    ## 严格早于当前时点：idx < i
    prev <- seq_len(i - 1L)
    if (length(prev) > 0L) {
      p48 <- prev[t[prev] >= t[i] - 48]
      p7  <- prev[t[prev] >= t[i] - 168]
      if (length(p48) > 0L) low48[i] <- min(v[p48])
      if (length(p7)  > 0L) low7[i]  <- min(v[p7])
    }
    stage[i] <- .kdigo_stage_scalar(v[i], low48[i], low7[i],
                                    use_48h = use_48h, use_ratio = use_ratio)
  }
  data.frame(t_hours = t, creat = v, low48 = low48, low7 = low7, stage = stage)
}

#' 检出界标后首次 KDIGO 肌酐事件
#'
#' @param dt data.frame(id, t_hours, value)，value 为血清肌酐 mg/dL
#' @param landmark_h 界标（小时），默认 24
#' @param horizon_h  结局窗上界（小时），默认 72
#' @param pre_from   界标前判定的起点（小时）。方案 4.2.6 原文为"界标前（0–24h）已满足
#'                   KDIGO 者"，故默认 0；设为 -168 即把院前实测值也纳入界标前判定，
#'                   作为敏感性分析（见 S 系列）。注意院前实测值**始终**参与 48h/7d
#'                   滚动比较基线，只用于构造比较值，与这个参数无关。
#' @param rrt ids 向量；这些 id 在界标后接受肾脏替代治疗者按 3 期计（次要结局）
#' @return data.frame(id, aki_pre, n_creat_pre, n_creat_post, aki, t_event,
#'                    stage_first, stage_max, low48_at_event, low7_at_event,
#'                    creat_first, creat_min_post)
#'   aki_pre = 界标前窗内已满足 KDIGO —— 主分析排除
#'   aki     = 界标之后、观察窗之内首次由肌酐检出 KDIGO AKI
kdigo_first_event <- function(dt, landmark_h = LANDMARK_H, horizon_h = HORIZON_H,
                              pre_from = 0, rrt = NULL,
                              use_48h = TRUE, use_ratio = TRUE) {
  dt <- merge_same_time(dt)
  if (nrow(dt) == 0L) return(NULL)

  res <- do.call(rbind, lapply(split(dt, dt$id), function(x) {
    s <- .kdigo_series(x$t_hours, x$value,
                       use_48h = use_48h, use_ratio = use_ratio)
    id <- x$id[1]

    pre  <- s$t_hours >= pre_from & s$t_hours <= landmark_h   # 界标前窗
    post <- s$t_hours >  landmark_h & s$t_hours <= horizon_h   # (24, 72]

    aki_pre <- any(pre & !is.na(s$stage) & s$stage >= 1L, na.rm = TRUE)
    ## 界标前（含院前）肌酐检测总次数。方案 3.4④ 要求"可依据实测肌酐判断界标前
    ## 是否已存在 AKI"——少于 2 次则构造不出任何比较值，属"无法判定"，须与
    ## "已有 AKI"分列。注意这与 pre_from 无关：院前实测值同样是有效比较基线。
    n_creat_before <- sum(s$t_hours <= landmark_h)

    idx_post <- which(post & !is.na(s$stage) & s$stage >= 1L)
    if (length(idx_post) > 0L) {
      i <- idx_post[1]
      aki <- TRUE; t_event <- s$t_hours[i]
      stage_first <- s$stage[i]; low48e <- s$low48[i]; low7e <- s$low7[i]
      creat_first <- s$creat[i]
    } else {
      aki <- FALSE; t_event <- NA_real_
      stage_first <- 0L; low48e <- NA_real_; low7e <- NA_real_; creat_first <- NA_real_
    }
    stage_max <- suppressWarnings(max(s$stage[post], na.rm = TRUE))
    if (!is.finite(stage_max)) stage_max <- 0L

    ## 界标后启动 RRT 者按 3 期计（次要结局口径）
    if (!is.null(rrt) && id %in% rrt) stage_max <- max(stage_max, 3L)

    data.frame(
      id = id, aki_pre = aki_pre,
      n_creat_pre    = sum(pre),
      n_creat_before = n_creat_before,
      n_creat_post   = sum(post),
      aki = aki, t_event = t_event,
      stage_first = stage_first, stage_max = stage_max,
      low48_at_event = low48e, low7_at_event = low7e,
      creat_first = creat_first,
      creat_min_post = if (any(post)) min(s$creat[post]) else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
  if (is.null(res)) return(NULL)
  rownames(res) <- NULL
  res[order(res$id), , drop = FALSE]
}

## ---- 4b. KDIGO 尿量标准（S8） -------------------------------------------------
## 主分析只用肌酐标准。尿量标准的独立价值在于：肌酐是滞后指标，且受肌肉量影响；
## 单用肌酐会漏掉相当一部分少尿型 AKI，这是审稿人对 AKI 研究最常见的质询之一。
## 本节的实现刻意对齐官方 mimic-code/urine_output_rate.sql 的口径：
##
##   1. 在每个**尿量记录时刻** t 上评估（不在固定网格上评估）。记录值被视为
##      "自上一记录以来"的累计尿量。
##   2. w 小时窗内取 (t−w+1, t]（官方用 5/11/23，因为 t 这一条本身覆盖 1 小时）。
##   3. 分母不是名义窗长 w，而是**实测覆盖时长**——窗内各条记录的间隔之和。
##      记录稀疏时名义窗长会严重高估覆盖、进而高估尿量速率；用真实覆盖时长
##      可避免把"没记录"误判成"尿得多"。覆盖时长不足 w 小时则本期不评估（NA）。
##   4. 阈值按 KDIGO 2012：6–12h <0.5、≥12h <0.5、≥24h <0.3（mL/kg/h），
##      以及 ≥12h 无尿。本函数只用 (24, 72] 内的评估点判结局。

UO_WS <- c(6, 12, 24)

#' 单个 id 的尿量速率序列（mL/kg/h），窗 6/12/24h
#'
#' @param t 数值向量，小时；@param v 数值向量，mL（可为负，冲洗液入量已记负）
#' @param weight 体重 kg（正数；缺失则全部返回 NA）
#' @return data.frame(t_hours, vol6, vol12, vol24, cov6, cov12, cov24,
#'                    rate6, rate12, rate24)
.uo_rate_series <- function(t, v, weight, ws = UO_WS) {
  o <- order(t); t <- t[o]; v <- v[o]
  n <- length(t)
  out <- data.frame(t_hours = t)
  if (!is.finite(weight) || weight <= 0) weight <- NA_real_
  ## 每条记录代表的时长：自上一记录起算；首条自 ICU 入科（t=0）起算。
  ## 与官方 tm_since_last_uo 同义（官方首条用入科至首次心率记录的时间）。
  gap <- c(t[1], diff(t))

  for (w in ws) {
    vol <- numeric(n); cov <- numeric(n)
    for (i in seq_len(n)) {
      ## 窗为 (t−w, t]，**两端都要限**。只写下界会把 t 之后的记录也算进来：
      ## 尿量恒定时总量与覆盖时长同比例放大，速率照样对，看不出错；一旦窗前
      ## 窗后尿量不同（少尿→恢复），就会把未来的尿灌进当前窗，把少尿抹平。
      idx <- which(t >= t[i] - (w - 1) & t <= t[i])
      vol[i] <- sum(v[idx])
      cov[i] <- sum(gap[idx])
    }
    out[[paste0("vol", w)]]  <- vol
    out[[paste0("cov", w)]]  <- cov
    ## 覆盖时长不足名义窗长 → 不评估（NA），而不是拿偏小的分母去算高速率
    out[[paste0("rate", w)]] <- ifelse(cov >= w, vol / (weight * cov), NA_real_)
  }
  out
}

#' KDIGO 尿量分期（单时点）
#' @return 0/1/2/3；三期判据不足时按可达的最高期计
.uo_stage_scalar <- function(rate6, rate12, rate24, vol12, cov12) {
  anuria12 <- !is.na(vol12) && !is.na(cov12) && cov12 >= 12 && vol12 <= 0
  if (anuria12) return(3L)
  if (!is.na(rate24) && rate24 < 0.3) return(3L)
  if (!is.na(rate12) && rate12 < 0.5) return(2L)
  if (!is.na(rate6)  && rate6  < 0.5) return(1L)
  0L
}

#' 检出界标后首次由尿量标准检出的 KDIGO AKI
#'
#' @param urine data.frame(id, t_hours, value)（value 为净尿量 mL）
#' @param weight 具名数值向量 id → kg
#' @return data.frame(id, uo_pre, uo_aki, uo_t_event, uo_stage_first,
#'                    uo_stage_max, n_uo_pre, n_uo_post)
kdigo_uo_first <- function(urine, weight, landmark_h = LANDMARK_H,
                           horizon_h = HORIZON_H, pre_from = 0) {
  if (is.null(urine) || nrow(urine) == 0L) return(NULL)
  ids <- unique(urine$id)
  res <- do.call(rbind, lapply(ids, function(k) {
    x <- urine[urine$id == k, , drop = FALSE]
    s <- .uo_rate_series(x$t_hours, x$value, unname(weight[k]))
    st <- vapply(seq_len(nrow(s)), function(i)
      .uo_stage_scalar(s$rate6[i], s$rate12[i], s$rate24[i],
                       s$vol12[i], s$cov12[i]), integer(1))
    t  <- s$t_hours
    pre  <- t >= pre_from & t <= landmark_h
    post <- t >  landmark_h & t <= horizon_h
    ## 可评估 = 该时点至少有一个窗的覆盖时长达到名义窗长。
    ## **必须逐时点区分"评了但没达标"与"压根没数据"**：覆盖不足时三个速率全为 NA、
    ## 分期给 0，若只看分期就会把"没有尿量记录"误当成"尿量正常"，
    ## 于是 uo_pre 全为 FALSE，"界标前是否已少尿"这件事被悄悄跳过。
    eval <- !is.na(s$rate6) | !is.na(s$rate12) | !is.na(s$rate24)
    idx_post <- which(post & st >= 1L)
    data.frame(
      id = k,
      uo_pre  = any(pre & st >= 1L),
      uo_aki  = length(idx_post) > 0L,
      uo_t_event = if (length(idx_post)) t[idx_post[1]] else NA_real_,
      uo_stage_first = if (length(idx_post)) st[idx_post[1]] else 0L,
      uo_stage_max   = if (any(post)) max(st[post]) else 0L,
      n_uo_pre  = sum(pre), n_uo_post = sum(post),
      n_uo_pre_eval  = sum(pre & eval),
      n_uo_post_eval = sum(post & eval),
      stringsAsFactors = FALSE)
  }))
  rownames(res) <- NULL
  res[order(res$id), , drop = FALSE]
}

## ---- 5. 组装分析人群（方案表 6 的逐层留存） ---------------------------------

#' 把队列、暴露、结局合成分析表
#' @param cohort data.frame(id, ...) 界标队列（已含表型与随访条件）
#' @param sodium 长表；@param creat 长表
#' @param pre_from 传给 kdigo_first_event 的界标前窗起点，默认 0（方案原文）
assemble <- function(cohort, sodium, creat, rrt = NULL, pre_from = 0) {
  na_df  <- delta_na(sodium)
  cr_df  <- kdigo_first_event(creat, rrt = rrt, pre_from = pre_from)

  d <- cohort
  d$has_delta_na <- d$id %in% if (is.null(na_df)) character(0) else na_df$id
  d <- merge(d, na_df, by = "id", all.x = TRUE, sort = FALSE)

  if (!is.null(cr_df)) {
    d <- merge(d, cr_df, by = "id", all.x = TRUE, sort = FALSE)
  } else {
    d$aki_pre <- NA; d$aki <- NA
  }
  ## 结局可判定：界标前未满足 KDIGO，且界标后确有肌酐复测、可构造无时序倒置的比较值
  d$evaluable <- !is.na(d$aki) & !d$aki_pre & d$n_creat_post > 0L
  d
}

## ---- 6. 数值呈现 --------------------------------------------------------------
## 放在共享层而不是 export.R：图与表报的是同一个 P、同一个 OR，格式化成两副面孔
## 就是自找麻烦——审稿人对着图 3 的 "0.00067" 去表 5 找 "0.0007" 是找不到的。
## 这两个函数必须在任何取用方（export.R 出表、figures.R 出图）之前定义。

#' 数值列 → "OR（95% CI）"字符串；OR 缺失时给'—'（如亚组的交互 P 行）
#'
#' 括号随输出语言走：英文版写成半角 "0.933 (0.903–0.965)"。图 3/图 4 的格子
#' 直接取这里的串，不经过 export.R 的 tr_fmt —— 若只在那里换标点，表格是半角、
#' 图里还留着全角括号，同一篇稿子里两种括号并存。
.fmt_or <- function(or, lo, hi) {
  fmt <- if (is_en()) "%.3f (%.3f–%.3f)" else "%.3f（%.3f–%.3f）"
  ifelse(is.na(or), "—", sprintf(fmt, or, lo, hi))
}

#' P 值格式化：小到 1e-4 以下记 <0.0001，避免正文出现 2.86e-05 这种写法
#'
#' **必须逐元素格式化**。整列交给 format.pval() 会按列里最小的那个数对齐小数位
#' ——一个 P=0.00067 会把整列都撑成 5 位小数（实测 "0.24281" "0.07744"）。
.fmt_p <- function(p) {
  ifelse(is.na(p), "—",
         ifelse(p < 1e-4, "<0.0001", sprintf("%.4f", p)))
}
