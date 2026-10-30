## =============================================================================
## models.R —— 主模型、剂量反应、预设亚组（严格按方案 4.4 的变量集合）
##
## 输入：cohort.R 产出的分析人群 data.frame（一行一 stay）。
## 暴露：delta_na（mmol/L，连续）。结局：aki（界标后首次由肌酐检出的 KDIGO AKI）。
##
## 三档模型（方案 4.4 / 5.3）：
##   未调整  delta_na
##   基础    + 年龄、性别、首次血钠、界标前最低肌酐、慢性肾病、糖尿病、
##             早期最低收缩压、高血压、冠心病、心房颤动、心脏瓣膜病
##   扩展    + 机械通气、血管活性药、利尿剂、液体出入量、血糖最差值、
##             血红蛋白最低值、心率与呼吸频率最差值、GCS 最低值
##
## 三条刻意的取舍：
##   (1) 血管活性药按方案表 9 同时纳入"是否使用"（vasoactive_any，二分类）与
##       "最大速率"（nee_max，连续）。二者**不共线**：nee_max 在使用者中取值连续
##       有别，在未使用者中恒为 0，两者并非彼此的线性函数（已用 qr 秩检验确认
##       设计阵满秩）。且 nee_max 覆盖 5 种血管加压药，vasoactive_any 覆盖含强心药
##       在内的 7 种，确有"用过强心药但加压药当量为 0"的正当取值。
##   (2) 剂量用**去甲肾上腺素当量**而非原样最大速率。7 种药量纲不可比
##       （多巴酚丁胺 20 vs 去甲肾上腺素 0.1 mcg/kg/min），同列取最大是错的；
##       详见 43_covariates_mimic.sql 中 cov_vaso_nee 的说明。
##   (3) 方案 4.4 明令不纳入含肌酐或血钠的整体严重程度评分（SAPS/APACHE/
##       OASIS/SOFA）。Charlson 指数亦不在方案变量集合内，故不进主模型。
##       治疗变量只进扩展模型，不作因果解释（方案 4.4）。
## =============================================================================

suppressMessages({ library(stats) })

## ---- 变量集合（提取前固定，不得按单变量 P 值或逐步回归增删） ----
BASIC_X <- c("age", "gender", "na_first", "creat_baseline", "ckd", "diabetes",
             "sbp_min_24h", "hypertension", "chd", "atrial_fib", "valvular")

EXTEND_X <- c("mech_vent_any", "vasoactive_any", "nee_max",
              "diuretic_furosemide_eq_mg",
              "fluid_in_24h_ml", "urine_out_24h_ml",
              "glucose_max_24h", "hgb_min_24h",
              "hr_max_24h", "rr_max_24h", "gcs_min_24h")

#' 由变量名拼出公式
#'
#' **环境必须显式设为调用者的帧**（parent.frame），不能用 as.formula 的默认值：
#' 默认会把公式的环境设成 .fml 自己的求值帧，那里只有 y 与 xs 两个形参。
#' 一旦某个参数不是从 data= 里取，而是从调用环境里取——例如
#' `glm(..., data = dd, weights = dd$..w)`——glm 会按
#' `eval(extras, data, env)` 先在模型框里找、再到 environment(formula) 里找，
#' 两处都没有 dd，于是报出与病因相隔很远的 "object 'dd' not found"。
#' 默认情况下 data= 里的变量仍然优先生效，所以这一改动只补上缺失的回退路径。
.fml <- function(y, xs)
  stats::as.formula(paste(y, "~", paste(xs, collapse = " + ")), env = parent.frame())

.pick <- function(xs, d) xs[xs %in% names(d)]

## ---- 统一建模人群 -------------------------------------------------------------

#' 参与建模的变量全集（暴露 + 结局 + 三档模型全部协变量）
MODEL_VARS <- c("aki", "delta_na", BASIC_X, EXTEND_X)

#' 三档模型共用的完整病例集合
#'
#' `glm` 遇缺失会**静默按行删除**，于是三档模型各用各的分母（未调整 2700、
#' 基础 2683、扩展 2656），同一张主表里 n 各不相同，读者无法横向比较，审稿人也
#' 会直接质疑"是不是挑着人跑的"。这里先算一次完整病例，由调用方**只切一次**，
#' 三档、RCS、亚组全部落在同一批人上，表里 n 只有一个。
#'
#' 代价是丢弃协变量缺失者（本例 1.6%）。缺失比例小且与结局无关的假设由多重插补
#' 敏感性分析（sensitivity.R 的 run_mi）单独检验，不在此处弥补。
#'
#' @return 逻辑向量，长度等于 nrow(d)
model_cases <- function(d) {
  vs <- .pick(MODEL_VARS, d)
  stats::complete.cases(d[, vs, drop = FALSE])
}

## ---- 三档 logistic -----------------------------------------------------------

#' @param d 分析人群（应已由 model_cases() 切过，保证三档同分母）
#' @return data.frame(模型, 项, OR, 下限, 上限, P, n, 事件数)
fit_tiers <- function(d) {
  tiers <- list(
    "未调整" = c("delta_na"),
    "基础"   = c("delta_na", .pick(BASIC_X, d)),
    "扩展"   = c("delta_na", .pick(BASIC_X, d), .pick(EXTEND_X, d))
  )
  out <- list()
  for (nm in names(tiers)) {
    xs <- tiers[[nm]]
    m  <- glm(.fml("aki", xs), data = d, family = binomial())
    co <- summary(m)$coefficients
    out[[nm]] <- data.frame(
      模型 = nm,
      项   = rownames(co),
      OR   = exp(co[, 1]),
      下限 = exp(co[, 1] - 1.96 * co[, 2]),
      上限 = exp(co[, 1] + 1.96 * co[, 2]),
      P    = co[, 4],
      n    = nobs(m), 事件数 = sum(d$aki),
      row.names = NULL)
  }
  res <- do.call(rbind, out)
  ## 兜底：若调用方忘了先用 model_cases() 切数据，这里必须吵出来——
  ## glm 的按行删除是静默的，不吵就会带着三个不同的分母进论文表。
  ns <- unique(res$n)
  if (length(ns) > 1L)
    warning(sprintf(
      "三档模型的分母不一致（%s）；应先以 model_cases() 取完整病例再拟合",
      paste(ns, collapse = " / ")), call. = FALSE)
  res
}

## ---- 限制性立方样条（剂量反应，方案 5.4）--------------------------------------
## 预设 4 节点，参考点 ΔNa = 0，节点位置按预先规定的分位数取，不得依结果调整。
RCS_KNOTS_Q <- c(0.05, 0.35, 0.65, 0.95)
RCS_REF     <- 0

#' @return list(model, knots, grid, p_overall, p_nonlinear, data, adj)
fit_rcs <- function(d, adj = .pick(BASIC_X, d)) {
  dd <- d[!is.na(d$delta_na) & !is.na(d$aki), , drop = FALSE]
  knots <- unname(quantile(dd$delta_na, RCS_KNOTS_Q))
  dd$delta_na_r <- dd$delta_na
  ## 字符型分类变量（MIMIC 的 gender 是 "M"/"F" 字符串）必须在拟合前转成 factor，
  ## 否则模型框里仍是 character；rcs_curve 复制单行构造预测数据时，character 会在
  ## model.matrix 里被转成只有 1 个水平的 factor 而报错。转成 factor 后水平随
  ## 数据整体固定，复制单行不掉水平。
  for (v in adj) if (is.character(dd[[v]])) dd[[v]] <- factor(dd[[v]])

  ## 显式构造限制性立方样条基，再把基矩阵当普通列放进数据框拟合。
  ##
  ## **为什么不用 rms::rcs() 写进公式**：那样 `model.matrix(terms, newdata)` 会在
  ## 新数据上**重新求值** `rcs()`，基函数与拟合时的基不一致（实测曲线塌成常数
  ## OR=1.000）。显式建基后，预测只是对已存系数做线性组合，结果可复现。
  ##
  ## `splines::ns()` 就是自然立方样条，即限制性立方样条（边界外线性），与
  ## rms::rcs 是同一个模型，只是参数化不同；节点由上面预先规定的分位数给定。
  inner <- knots[-c(1, length(knots))]
  bnd   <- range(knots)
  B <- splines::ns(dd$delta_na, knots = inner, Boundary.knots = bnd)
  bs <- paste0("s", seq_len(ncol(B)))
  for (j in seq_along(bs)) dd[[bs[j]]] <- B[, j]

  m     <- glm(.fml("aki", c(bs, adj)),       data = dd, family = binomial())
  m_lin <- glm(.fml("aki", c("delta_na", adj)), data = dd, family = binomial())
  m_no  <- glm(.fml("aki", adj),                data = dd, family = binomial())
  p_overall   <- anova(m_no,  m,     test = "LRT")$`Pr(>Chi)`[2]
  p_nonlinear <- anova(m_lin, m,     test = "LRT")$`Pr(>Chi)`[2]

  grid <- seq(quantile(dd$delta_na, 0.01), quantile(dd$delta_na, 0.99), length.out = 200)
  list(model = m, knots = knots, inner = inner, bnd = bnd, basis = bs,
       grid = grid, p_overall = p_overall, p_nonlinear = p_nonlinear,
       data = dd, adj = adj)
}

#' 用预测值构造相对参考点的 OR 曲线（RCS 作图用）
#'
#' 预测数据必须从**模型框**（glm 的 `$model`）构造，不能用原始 data.frame 复制
#' 单行：原始表里 gender 一类分类变量是字符型，复制成 200 行同值后经 model.matrix
#' 转 factor 只剩 1 个水平，直接报 "contrasts can be applied only to factors with
#' 2 or more levels"。模型框里的分类变量已是带完整水平的 factor，复制不会掉水平。
#'
#' @return data.frame(delta_na, OR, 下限, 上限)
rcs_curve <- function(fit, conf = 0.95) {
  m  <- fit$model
  bs <- fit$basis
  ## 用**同一组节点**在网格上重建基，故与拟合时完全一致
  B  <- splines::ns(fit$grid, knots = fit$inner, Boundary.knots = fit$bnd)
  B0 <- splines::ns(RCS_REF,  knots = fit$inner, Boundary.knots = fit$bnd)
  ## 模型里样条项与 adj 之间没有交互，adj 在"相对参考点"的差值里精确抵消，
  ## 因此只需取样条项系数，不必构造协变量矩阵。
  bco <- coef(m)[bs]
  V   <- vcov(m)[bs, bs, drop = FALSE]
  lp  <- as.numeric(B  %*% bco)
  lp0 <- as.numeric(B0 %*% bco)
  d   <- sweep(B, 2, as.numeric(B0), "-")
  se  <- sqrt(rowSums((d %*% V) * d))
  z   <- qnorm(1 - (1 - conf) / 2)
  data.frame(delta_na = fit$grid,
             OR   = exp(lp - lp0),
             下限 = exp(lp - lp0 - z * se),
             上限 = exp(lp - lp0 + z * se))
}

## ---- 预设亚组（方案 5.5：仅两个，结论一律标注探索性）--------------------------
SUBGROUPS <- list(
  "首次血钠 <135 vs ≥135" = function(d) d$na_first < 135,
  "有无慢性肾病"           = function(d) d$ckd == 1
)

## 亚组各层的**显示名**。层 = "是" 一律表示上面那个条件成立。
##
## 这张对照表必须与 SUBGROUPS 紧挨着放：图与表若从亚组名里切字符串推标签，
## 一定会错，而且错得看不出来——
##   ① "首次血钠 <135 vs ≥135" 的 "是" 层其实是 **<135 组**，按 "vs" 切出来
##      得到 "≥135"，正好**标反**；
##   ② "有无慢性肾病" 里根本没有 "vs"，切出来还是原名，于是图例成了
##      "有无慢性肾病：有无慢性肾病（是）"。
## ①是方向性错误（读者会把 <135 组的估计安到 ≥135 组头上），实测已发生过。
SUBGROUP_LEVELS <- list(
  "首次血钠 <135 vs ≥135" = c("是" = "首次血钠 <135", "否" = "首次血钠 ≥135"),
  "有无慢性肾病"          = c("是" = "有慢性肾病",     "否" = "无慢性肾病")
)

#' 亚组成员行的显示标签（图与表共用，保证两处口径一致）
#' @param 亚组,层  fit_subgroups() 输出的两列
subgroup_labels <- function(亚组, 层) {
  mapply(function(a, l) {
    lv <- SUBGROUP_LEVELS[[a]]
    if (is.null(lv) || !l %in% names(lv)) paste0(a, "：", l) else unname(lv[[l]])
  }, as.character(亚组), as.character(层), USE.NAMES = FALSE)
}

#' 亚组内跑调整模型，并与其余人群做交互检验
fit_subgroups <- function(d, adj = .pick(BASIC_X, d)) {
  res <- list()
  for (nm in names(SUBGROUPS)) {
    g <- SUBGROUPS[[nm]](d)
    if (all(is.na(g))) next
    for (lv in c(TRUE, FALSE)) {
      sub <- d[!is.na(g) & g == lv, , drop = FALSE]
      if (nrow(sub) < 30 || sum(sub$aki) < 10) next
      m  <- glm(.fml("aki", c("delta_na", adj)), data = sub, family = binomial())
      co <- summary(m)$coefficients["delta_na", ]
      res[[length(res) + 1L]] <- data.frame(
        亚组 = nm, 层 = if (lv) "是" else "否",
        n = nrow(sub), 事件数 = sum(sub$aki),
        OR = exp(co[1]), 下限 = exp(co[1] - 1.96 * co[2]),
        上限 = exp(co[1] + 1.96 * co[2]), P = co[4], row.names = NULL)
    }
    d$..g <- as.integer(g)
    mi <- glm(.fml("aki", c("delta_na", "..g", "delta_na:..g", adj)),
              data = d, family = binomial())
    co <- summary(mi)$coefficients
    if ("delta_na:..g" %in% rownames(co))
      res[[length(res) + 1L]] <- data.frame(
        亚组 = nm, 层 = "交互 P", n = nrow(d), 事件数 = sum(d$aki),
        OR = NA, 下限 = NA, 上限 = NA, P = co["delta_na:..g", 4], row.names = NULL)
  }
  do.call(rbind, res)
}
