## =============================================================================
## cluster_eicu.R —— eICU 医院聚类稳健标准误（外审 P0 补充分析）
##
## 动机：eICU-CRD 是多中心数据库，同院患者的检验习惯、诊疗流程、AKI 判定
## 倾向相关，普通 logistic 的独立同分布假设偏乐观（标准误偏小、CI 偏窄）。
## 本文件用 sandwich 聚类稳健协方差（按 hospitalid 聚类）重估基础/扩展模型
## 的 ΔNa 系数，与普通标准误版并排输出；不改变点估计（同一 glm）。
##
## 前置条件（2026-10-31 备妥、待原始数据到货后执行）：
##   30_cohort_eicu.sql 已加 `p.hospitalid,` 一列并重跑 SQL 提取层，
##   使 eicu_C_full.rds$analysis 含 hospitalid。旧 RDS（v1.0.1 归档版）无此列，
##   本脚本会明确报错而不是静默跳过。
##
## GitHub/Zenodo v1.0.1 仓库不含本文件与该 SQL 改动；聚类分析属服务器/
## 本地复现扩展层，结果进稿件前须先回填正文与 eTable。
## =============================================================================

suppressMessages({ library(sandwich); library(lmtest) })

#' 医院聚类稳健 logistic 结果（ΔNa 行）
#'
#' @param d 已按 model_cases() 切过的分析人群，须含 hospitalid
#' @param xs 协变量名向量（暴露 delta_na 由调用方放入或在此追加）
#' @param exposure 暴露列名，默认 "delta_na"
#' @return data.frame(模型, 项, SE类型, OR, 下限, 上限, P, n, 事件数, 医院数)
fit_cluster_or <- function(d, xs, exposure = "delta_na", label = "") {
  if (!"hospitalid" %in% names(d))
    stop("analysis 缺 hospitalid 列：请先在 30_cohort_eicu.sql 加 `p.hospitalid,` 并重跑提取层（见文件头注）")
  if (any(is.na(d$hospitalid)))
    stop("hospitalid 存在缺失：检查 patient 表关联")
  xs <- unique(c(exposure, xs))
  m <- glm(stats::as.formula(paste("aki ~", paste(xs, collapse = " + "))),
           data = d, family = binomial())
  v_plain <- stats::vcov(m)
  v_cl    <- sandwich::vcovCL(m, cluster = ~ hospitalid)
  row <- function(v, tag) {
    s <- lmtest::coeftest(m, vcov. = v)[exposure, ]
    data.frame(模型 = label, 项 = exposure, SE类型 = tag,
               OR = unname(exp(s[1])),
               下限 = unname(exp(s[1] - 1.96 * s[2])),
               上限 = unname(exp(s[1] + 1.96 * s[2])),
               P = unname(s[4]),
               n = nobs(m), 事件数 = sum(d$aki),
               医院数 = length(unique(d$hospitalid)),
               row.names = NULL)
  }
  rbind(row(v_plain, "普通"), row(v_cl, "聚类稳健"))
}

#' 医院层面概览（每院 n / 事件率的离散度，供审稿回复引用）
hospital_overview <- function(d) {
  tab <- aggregate(cbind(n = rep(1, nrow(d)), aki = d$aki),
                   by = list(hospitalid = d$hospitalid), FUN = sum)
  tab$rate <- tab$aki / tab$n
  list(医院数 = nrow(tab),
       中位每院n = stats::median(tab$n),
       每院n四分位 = stats::quantile(tab$n, c(.25, .75)),
       院事件率范围 = range(tab$rate),
       头部5院占比 = sum(head(sort(tab$n, decreasing = TRUE), 5)) / sum(tab$n))
}
