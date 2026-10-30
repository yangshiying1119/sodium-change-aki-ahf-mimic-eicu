## =============================================================================
## advance.R —— 三项补充分析：严重度梯度、绝对风险差、E 值
##
## 三者回答的是同一类问题：主分析给出"每升 1 mmol/L，OR 0.93"之后，
##   ① 效应是否随 AKI 严重度递增（剂量-反应的又一种形式，支持因果解释）；
##   ② 翻译成人话是多少个百分点的绝对风险（临床读者读得懂的量级）；
##   ③ 需要多强的未测混杂才能把这点关联抹平（结论有多脆）。
##
## 与主分析共用 BASIC_X 与统一完整病例，不另设变量集合、不按 P 值增删。
## =============================================================================

source("models.R")

## ---- 1. 纯肌酐的 AKI 分期 ----------------------------------------------------
## 主流程的 stage_max 在 common.R 里被 RRT 抬到 3 期（"界标后启动 RRT 者按 3 期计"），
## 于是存在 aki = FALSE 但 stage_max = 3 的病例（MIMIC 94 例）。拿它做严重度分析
## 会把未发生 AKI 的人算进重度组。这里以 rrt = NULL 重算一份**只由肌酐定义**的
## 最大分期，与主结局 aki 严格等价：stage_creat >= 1 当且仅当 aki = TRUE。
## 不改动 common.R 与既有产物，避免牵动整条流水线。
stage_creat_only <- function(ex, ids = NULL) {
  ## **必须先按 id 过滤再进判定**：提取层的肌酐长表覆盖全部 ICU 住院
  ## （MIMIC 94,209 个 id、eICU 190,331 个），而 kdigo_first_event 是逐 id 的
  ## R 循环，全量跑要十几分钟；分析人群只有两三千人，先切再算不到一秒。
  cr <- if (is.null(ids)) ex$creatinine
        else ex$creatinine[ex$creatinine$id %in% ids, , drop = FALSE]
  s <- kdigo_first_event(cr, rrt = NULL)
  s[, c("id", "stage_max"), drop = FALSE]
}

#' 把纯肌酐分期并进分析人群，并断言与 aki 完全一致
#'
#' 断言不是形式主义：一旦两者对不上（比如日后改了 kdigo_first_event 的窗口），
#' 下面三档严重度模型的分母就会与主分析不同，而表面上一切正常。
attach_stage <- function(a, ex) {
  sc <- stage_creat_only(ex, a$id)
  a$stage_creat <- sc$stage_max[match(a$id, sc$id)]
  a$stage_creat[is.na(a$stage_creat)] <- 0L
  bad <- sum((a$stage_creat >= 1L) != a$aki)
  if (bad > 0)
    stop(sprintf("纯肌酐分期与主结局不一致：%d 例", bad), call. = FALSE)
  a
}

## ---- 2. 严重度梯度 ------------------------------------------------------------

SEVERITY_CUTS <- list(
  list(名 = "任何 AKI（≥1 期）",   cut = 1L),
  list(名 = "中重度 AKI（≥2 期）", cut = 2L),
  list(名 = "重度 AKI（3 期）",    cut = 3L)
)

#' 三个累积阈值各拟合一次基础模型
#'
#' 用累积阈值而不是有序 logistic：① 无须比例优势假设；② 三行可直接读出
#' "效应是否随严重度变强"，读者不需要理解比例优势模型。
#'
#' @return data.frame(结局, n, 事件数, OR, 下限, 上限, P)
fit_severity <- function(d, adj = .pick(BASIC_X, d)) {
  out <- lapply(SEVERITY_CUTS, function(s) {
    d$..y <- as.integer(d$stage_creat >= s$cut)
    m  <- glm(.fml("..y", c("delta_na", adj)), data = d, family = binomial())
    co <- summary(m)$coefficients["delta_na", ]
    data.frame(结局 = s$名, n = nobs(m), 事件数 = sum(d$..y),
               OR = exp(co[1]), 下限 = exp(co[1] - 1.96 * co[2]),
               上限 = exp(co[1] + 1.96 * co[2]), P = co[4],
               row.names = NULL)
  })
  do.call(rbind, out)
}

#' 有序 logistic（4 级：无 / 1 期 / 2 期 / 3 期），单个 OR 概括整条梯度
#'
#' 只在各阈值事件数都够时拟合；否则返回 NULL，由调用方跳过。
fit_ordinal <- function(d, adj = .pick(BASIC_X, d)) {
  d$..sev <- factor(pmin(d$stage_creat, 3L), levels = 0:3, ordered = TRUE)
  if (any(table(d$..sev) < 10L)) return(NULL)
  m <- try(MASS::polr(.fml("..sev", c("delta_na", adj)), data = d, Hess = TRUE),
           silent = TRUE)
  if (inherits(m, "try-error")) return(NULL)
  co <- summary(m)$coefficients
  est <- co["delta_na", "Value"]; se <- co["delta_na", "Std. Error"]
  data.frame(结局 = "有序结局（无 / 1 期 / 2 期 / 3 期，比例优势）",
             n = nobs(m), 事件数 = sum(d$stage_creat >= 1L),
             OR = exp(est), 下限 = exp(est - 1.96 * se),
             上限 = exp(est + 1.96 * se),
             P = 2 * pnorm(-abs(est / se)), row.names = NULL)
}

## ---- 3. 绝对风险差（边际标准化 / g-computation）-------------------------------

#' 把给定的 ΔNa 值代入每个个体，取预测概率的平均
#'
#' 用模型框 m$model 而不是原始数据：模型框里的分类变量已是带完整水平的 factor，
#' 且已按 glm 的规则去掉缺失行，复制后不会掉水平、也不会混进未参与拟合的人。
std_risk <- function(m, at) {
  nd <- m$model
  nd$delta_na <- at
  mean(stats::predict(m, newdata = nd, type = "response"))
}

#' 两个 ΔNa 取值之间的标准化绝对风险差及其 bootstrap 区间
#'
#' 非参数 bootstrap（按行重抽后重新拟合）：边际标准化是两个平均预测值之差，
#' 没有现成的解析方差可用；B = 200 在本样本量下对三位有效数字足够。
#' bootstrap 同时给出两个端点的标准化风险区间，便于表中并列。
#'
#' @param from,to ΔNa 的两个对比取值
#' @return list(风险_from, 风险_to, 差值, 区间, B)
risk_contrast <- function(d, from, to, adj = .pick(BASIC_X, d), B = 200L) {
  vars <- c("aki", "delta_na", adj)
  dd <- d[stats::complete.cases(d[, vars, drop = FALSE]), vars, drop = FALSE]
  fit <- function(x) glm(.fml("aki", c("delta_na", adj)), data = x, family = binomial())
  m  <- fit(dd)
  pt <- c(from = std_risk(m, from), to = std_risk(m, to))

  set.seed(20261025)
  bs <- vapply(seq_len(B), function(i) {
    x <- fit(dd[sample.int(nrow(dd), replace = TRUE), , drop = FALSE])
    c(from = std_risk(x, from), to = std_risk(x, to))
  }, numeric(2))

  q <- function(v) unname(stats::quantile(v, c(0.025, 0.975)))
  list(from = pt[["from"]], to = pt[["to"]],
       from_ci = q(bs["from", ]), to_ci = q(bs["to", ]),
       diff = pt[["to"]] - pt[["from"]],
       diff_ci = q(bs["to", ] - bs["from", ]),
       B = B)
}

## ---- 4. E 值 ------------------------------------------------------------------

#' OR 按基线风险换算成近似 RR（VanderWeele & Ding 的通用换算）
#'
#' 本结局发生率 26.6%，属常见结局，OR 会明显放大 RR；直接把 OR 当 RR 代入
#' E 值公式会**低估**所需混杂强度（方向不利于自曝）。故先换算。
or_to_rr <- function(or, p0) or / (1 - p0 + p0 * or)

#' E 值：把 RR 折向 1 后按 VW-Ding 公式计算
evalue_rr <- function(rr) {
  rr <- ifelse(rr < 1, 1 / rr, rr)
  rr + sqrt(rr * (rr - 1))
}

#' 由 OR 与基线风险给出点估计与 CI 的 E 值
#'
#' 连续暴露的 E 值**依赖所取的对比尺度**：同一条模型，每 1 mmol/L 与每 3 mmol/L
#' 的 E 值完全不同。只报 1 mmol/L 会低估稳健性（因为没人认为 1 mmol/L 是有临床
#' 意义的对比），只报 3 mmol/L 又像是在挑好看的尺度。故两者并列，并在表注里
#' 写明 E 值对应的是哪一个对比。
#'
#' @param k 对比的 ΔNa 跨度（mmol/L），点估计取 OR^k、区间取 lo^k / hi^k
#' @return data.frame(对比, 量, 值) 供出表
evalue_summary <- function(or, lo, hi, p0, k = 1L) {
  o  <- or^k
  lk <- lo^k; hk <- hi^k
  ## CI 的 E 值取**离 1 更近**的那一端：OR < 1 时上限更靠近 1，OR > 1 时下限更靠近。
  ## 方向写反会让 CI 的 E 值反而大于点估计——看起来"区间比点估计更稳健"，
  ## 是明显的内部矛盾。
  ci_lim <- if (o < 1) hk else lk
  ## 区间跨 1 时，按惯例 E 值记为 1（数据本身不排除无效应）
  e_ci <- if (lk < 1 && hk > 1) 1 else evalue_rr(or_to_rr(ci_lim, p0))
  rr   <- or_to_rr(o, p0)
  data.frame(
    对比 = sprintf("ΔNa 每升高 %d mmol/L", k),
    量 = c("OR（基础模型）", "换算为近似 RR", "E 值（点估计）", "E 值（95% CI 端点）"),
    值 = c(sprintf("%.3f", o),
           sprintf("%.3f", rr),
           sprintf("%.2f", evalue_rr(rr)),
           sprintf("%.2f", e_ci)),
    row.names = NULL)
}

## ---- 5. 打包 ------------------------------------------------------------------

#' 基础模型下 ΔNa 的效应量（与 fit_tiers 的"基础"行同一条公式，用于 E 值）
base_effect <- function(d, adj = .pick(BASIC_X, d)) {
  m  <- glm(.fml("aki", c("delta_na", adj)), data = d, family = binomial())
  co <- summary(m)$coefficients["delta_na", ]
  list(OR = exp(co[1]), 下限 = exp(co[1] - 1.96 * co[2]),
       上限 = exp(co[1] + 1.96 * co[2]), P = co[4], n = nobs(m))
}

#' 一次算完三块，返回可直接出表的结果对象
#' @param a  已并协变量、已取完整病例的分析人群（应已 attach_stage）
#' @param ex 提取结果（用于 stage_creat_only，若 a 已含 stage_creat 可省）
run_advance <- function(a, ex = NULL) {
  if (is.null(a$stage_creat)) {
    if (is.null(ex)) stop("a 中没有 stage_creat，且未提供 ex")
    a <- attach_stage(a, ex)
  }
  sev <- do.call(rbind, Filter(Negate(is.null),
                               list(fit_severity(a), fit_ordinal(a))))
  rownames(sev) <- NULL

  p0 <- mean(a$aki)
  be <- base_effect(a)

  list(严重度 = sev,
       风险 = list(`0 vs +3`  = risk_contrast(a, 0, 3),
                   `-3 vs +3` = risk_contrast(a, -3, 3),
                   `0 vs +6`  = risk_contrast(a, 0, 6)),
       基础效应 = be,
       E值 = rbind(evalue_summary(be$OR, be$下限, be$上限, p0, k = 1L),
                   evalue_summary(be$OR, be$下限, be$上限, p0, k = 3L)),
       基线风险 = p0, n = nrow(a), 事件数 = sum(a$aki))
}
