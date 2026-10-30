## =============================================================================
## cohort.R —— 队列组装与逐层留存（方案表 6 的 N0–N8）
##
## 职责：调用提取层拿到长表 → 交给 common.R 的共享判定层 → 输出分析人群。
## 本文件不含判定逻辑（判定全在 common.R），也不含库特定 SQL（提取全在 sql/）。
##
## 纳入（方案 3.2）与排除（方案 3.4）逐条对应：
##   纳入① 成人 ≥18 岁                      → sql/30（N1）
##   纳入② 首次符合条件的 ICU 住院            → 本文件末（N3，见下方说明）
##   纳入③ 急性心衰表型                      → sql/30（N2）
##   纳入④ 界标时仍在 ICU 且存活              → 本文件（N4）
##   纳入⑤ 0–24h ≥2 次血钠（时刻不同）        → common.R::delta_na（N5）
##   纳入⑥ 界标前可判定 + 界标后有比较值       → 本文件（N6/N7）
##   排除① 透析依赖 / 终末期肾病              → sql/30 + sql/42
##   排除② 界标前已满足 KDIGO                → common.R::kdigo_first_event（N6）
##   排除③ 界标前已开始 RRT                   → sql/42 → 本文件
##   排除④ 缺乏肌酐比较基线                   → 与 N6 的"无法判定"同源
##   排除⑤⑥ 时间戳/单位异常、血钠超范围        → sql/20 的值域与单位断言
##
## N3 的位置：方案 3.5 的原文是"取第一个同时满足纳入标准全部条目、且不触发任何
## 排除标准者"，即首次住院的判定要作用在**完全合格集合**上，故本文件把它放在最后
## 一步。这与流程图里 N3 的编号位置不同，但含义更严格；"因非首次而排除"的人数
## 单独列出（方案 3.5 明文要求报告）。
## =============================================================================

## 脚本所在目录（不能依赖 extract.R 里的同名变量，会先于它求值）
.here <- local({
  a <- commandArgs(trailingOnly = FALSE)
  s <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(s) && nzchar(s)) dirname(normalizePath(s)) else "."
})
source(file.path(.here, "extract.R"))
source(file.path(.here, "common.R"))

#' 跑一遍 MIMIC-IV 主队列
#' @param cohort 表型口径 "A" / "B" / "C"（主分析用 "B"）
#' @param restrict 是否把长表预先裁到队列内的 id（省算力，不改变结果）
#' @param pre_from 界标前判定窗起点，默认 0（方案 4.2.6 原文"界标前（0–24h）"）；
#'                 设 -PRE_H 可复现含院前实测值的版本，作 S 系列敏感性分析
#' @param ex 已提取的列表（extract_mimic() / extract_mimic_all() 的返回值）。
#'           为 NULL 时自行提取；跑批时传入可避免重复扫 chartevents。
#' @param na_fun 暴露构造函数，默认 delta_na（0–24h 末次 − 首次）。S1/S2 传变体。
#' @param cr_fun 结局构造函数，默认按 pre_from/rrt_post 调 kdigo_first_event。S7 传变体。
#' @param first_only 是否施加"首次符合条件的住院"限制（纳入②）。S10 传 FALSE。
#' @param uo 仅 S8 使用：list(urine = 长表, weight = 具名向量 id→kg)。给了就按
#'   "肌酐标准 ∪ 尿量标准"的复合口径判定结局与界标前状态；不给则纯粹按肌酐（主分析）。
run_mimic_cohort <- function(cohort = "B", restrict = TRUE, pre_from = 0, ex = NULL,
                             na_fun = NULL, cr_fun = NULL, first_only = TRUE,
                             uo = NULL) {
  if (is.null(ex)) ex <- extract_mimic(cohort = cohort)
  ## 口径切换：优先用 ex 里带的对应表型表。
  ## **不能只看 ex$cohort**——传入 ex 时 cohort 参数原本被完全忽略，
  ## 于是 cohort = "A" 拿到的是 B 的人群（参见 extract.R::.all_cohorts 的说明）。
  coh <- if (!is.null(ex$cohorts) && cohort %in% names(ex$cohorts))
           ex$cohorts[[cohort]]
         else ex$cohort
  if (is.null(coh))
    stop(sprintf("提取结果中没有表型 %s 的人群表", cohort))
  ids <- coh$id

  sod <- ex$sodium;     if (restrict) sod <- sod[sod$id %in% ids, , drop = FALSE]
  cr  <- ex$creatinine; if (restrict) cr  <- cr[cr$id  %in% ids, , drop = FALSE]
  rrt <- ex$rrt;        if (restrict) rrt <- rrt[rrt$id %in% ids, , drop = FALSE]

  na_df <- if (is.null(na_fun)) delta_na(sod) else na_fun(sod)

  ## 界标前最低实测肌酐（方案表 9"基线肌酐 = 界标前最低实测值"）。
  ## 含院前实测值：方案 4.2.2 明确把入院前 7 天内的实测肌酐纳入比较。
  base_cr <- do.call(rbind, lapply(split(cr, cr$id), function(x) {
    x <- x[x$t_hours <= LANDMARK_H & !is.na(x$value), , drop = FALSE]
    if (!nrow(x)) return(NULL)
    data.frame(id = x$id[1], creat_baseline = min(x$value), stringsAsFactors = FALSE)
  }))

  ## RRT：界标前已在做（排除③）与界标后启动（次要结局按 3 期计，方案表 8）
  ## 两次调用 kdigo_first_event：先不带 rrt 拿到 aki_pre（RRT 不影响"界标前已有 AKI"
  ## 的肌酐判定，而排除③已单独处理），再把 rrt_post 传进去拿分期。
  rrt_pre  <- unique(rrt$id[rrt$t_hours <= LANDMARK_H])
  rrt_post <- unique(rrt$id[rrt$t_hours >  LANDMARK_H & rrt$t_hours <= HORIZON_H])
  cr_df <- if (is.null(cr_fun))
    kdigo_first_event(cr, pre_from = pre_from, rrt = rrt_post) else cr_fun(cr)

  d <- coh
  d <- merge(d, na_df,   by = "id", all.x = TRUE, sort = FALSE)
  d <- merge(d, cr_df,   by = "id", all.x = TRUE, sort = FALSE)
  d <- merge(d, base_cr, by = "id", all.x = TRUE, sort = FALSE)

  ## ---- S8 复合口径：肌酐标准 ∪ 尿量标准 ----------------------------------------
  ## 只在传入 uo 时启用。肌酐与尿量是**并列的两套 KDIGO 判据**，任一阳性即为 AKI，
  ## 界标前的"已满足 KDIGO"与界标后的结局都必须同样处理，否则会出现"界标前已少尿
  ## 但被留在队列里"这种口径错位。
  d$aki_creat <- d$aki; d$aki_pre_creat <- d$aki_pre
  d$comp <- FALSE
  if (!is.null(uo)) {
    uo_w <- uo$weight; names(uo_w) <- as.character(names(uo_w))
    uo_df <- kdigo_uo_first(uo$urine, uo_w, pre_from = pre_from)
    if (is.null(uo_df)) stop("S8：尿量长表为空，无法按尿量标准判定")
    uo_df$id <- as.character(uo_df$id)
    d <- merge(d, uo_df, by = "id", all.x = TRUE, sort = FALSE)
    for (v in c("uo_pre", "uo_aki")) d[[v]][is.na(d[[v]])] <- FALSE
    for (v in c("n_uo_pre", "n_uo_post", "n_uo_pre_eval", "n_uo_post_eval"))
      d[[v]][is.na(d[[v]])] <- 0L
    ## 取材时点取两套标准中较早者；两者都未发生时自然都为 NA
    t1 <- d$t_event; t2 <- d$uo_t_event
    d$t_event <- ifelse(is.na(t1), t2, ifelse(is.na(t2), t1, pmin(t1, t2)))
    d$aki_pre <- d$aki_pre_creat | d$uo_pre
    d$aki     <- d$aki_creat     | d$uo_aki
    d$comp    <- TRUE
  }

  ## ---- 各项标志 ----
  d$has_delta_na   <- !is.na(d$delta_na)
  d$esrd_icd       <- d$esrd_excluded == 1
  d$esrd_chart     <- d$id %in% ex$esrd_chart
  d$esrd           <- d$esrd_icd | d$esrd_chart
  d$rrt_pre        <- d$id %in% rrt_pre
  d$rrt_post       <- d$id %in% rrt_post
  d$landmark_ok    <- d$still_in_icu == 1 & d$died_within_24h == 0
  ## 界标前（含院前）肌酐 < 2 次 → 构造不出任何比较值，判不了界标前有无 AKI。
  ## 复合口径下，尿量标准不依赖"比较基线"（用的是绝对阈值），只要界标前有可评估的
  ## 尿量速率就能判，故两条路径任一可行即视为可判定。
  ## ⚠️ 不能写成 `d$comp & d$n_uo_pre_eval > 0L`：单用肌酐时 uo 列为 NULL，
  ## `NULL > 0L` 得到 logical(0)，而 `logical(n) & logical(0)` 在 R 里**不报错**，
  ## 直接塌成 logical(0)，一路带到下面 `$<-` 才炸成"replacement has 0 rows"。
  ## 用取值函数把"列不存在"显式落成全 FALSE，错误就不会藏在几行之外。
  uo_col <- function(v) if (is.null(d[[v]])) rep(0L, nrow(d)) else d[[v]]
  uo_ok_pre  <- d$comp & uo_col("n_uo_pre_eval")  > 0L
  uo_ok_post <- d$comp & uo_col("n_uo_post_eval") > 0L
  d$not_judgeable  <- (is.na(d$n_creat_before) | d$n_creat_before < 2L) & !uo_ok_pre
  ## 界标后是否有可判定的观察（肌酐复测，或复合口径下的可评估尿量速率）。
  ## 单独存一列，供留存的 N7 关口与 evaluable 共用，避免两处各写一遍条件后走样。
  d$has_post_obs   <- (!is.na(d$n_creat_post) & d$n_creat_post > 0L) | uo_ok_post
  ## ⚠️ 用 `%in% FALSE` 而不是 `!`：aki_pre 在"界标前/界标时一条肌酐都没有"时为 NA，
  ## 而 `!NA` 是 NA。此时若同行的 has_delta_na 恰为 TRUE、not_judgeable 又为 FALSE，
  ## 与运算就不再被别的 FALSE 短路吸收，evaluable 变成 NA——接着
  ## `d[d$evaluable, ]` 用含 NA 的逻辑下标取子集，会**静默插入整行 NA**，
  ## 分析人群凭空多出若干行且无人察觉。要求"已知为否"可彻底堵住这条路径。
  ## （现有 MIMIC 数据里 8 个 NA 恰好都被 has_delta_na=FALSE 短路掉，故为等价改写。）
  d$evaluable      <- d$has_delta_na %in% TRUE & d$aki_pre %in% FALSE &
                      d$not_judgeable %in% FALSE & d$has_post_obs %in% TRUE

  ## ---- 逐层留存（按方案表 6 的节点名，按实际判定顺序筛） ----
  ## 每一行同时给出"本节点排除人数"与"排除后留存人数"两列。
  ## **不能只用一列**：前三个节点（N0/N1/N2）拿到的是**留存数**，其后各节点
  ## 算出的是**排除数**，混在一列里会让 留存 = N0 − cumsum(人数) 变成负数
  ## （早先版本实测末行为 −15,124），流程图与正文都无法自洽。
  keep  <- rep(TRUE, nrow(d))
  n0    <- ex$attrition[["N0_icu_stays_all"]]
  n1    <- ex$attrition[["N1_adult"]]
  n2    <- nrow(d)
  n_run <- n2
  at <- list(
    data.frame(节点 = "N0 数据库内全部 ICU 住院", 排除原因 = "—",
               排除 = 0L, 留存 = n0, stringsAsFactors = FALSE),
    data.frame(节点 = "N1 成人（≥18 岁）", 排除原因 = "未满 18 岁",
               排除 = n0 - n1, 留存 = n1, stringsAsFactors = FALSE),
    data.frame(节点 = sprintf("N2 急性心衰表型 %s", cohort), 排除原因 = "不符合表型",
               排除 = n1 - n2, 留存 = n2, stringsAsFactors = FALSE))

  ## 每调用一次：记下本节点排除数与排除后的留存数，并把这一步的排除并进 keep
  step <- function(cond, label, reason) {
    hit <- keep & cond
    n_run <<- n_run - sum(hit)
    data.frame(节点 = label, 排除原因 = reason, 排除 = sum(hit),
               留存 = n_run, stringsAsFactors = FALSE)
  }

  ## 排除① 透析依赖 / 终末期肾病（两条路径分列，便于核查边际贡献）
  ## 两条路径的**来源表**两库不同（MIMIC：diagnoses_icd + chartevents 问诊项；
  ## eICU：diagnosis 码 + treatment 的慢性肾衰透析记录），故节点名随库切换。
  ## 判定逻辑本身完全一致，差别只在"这条路径取自哪张表"。
  lab_a <- if (identical(ex$db, "eicu"))
    "　└ 排除①a 诊断码（N18.6/585.6/Z99.2/V45.1x）"
  else "　└ 排除①a 诊断码（N186/5856/Z992/V4511/V4512）"
  lab_b <- if (identical(ex$db, "eicu"))
    "　└ 排除①b 慢性肾衰透析治疗（treatment）"
  else "　└ 排除①b 入科问诊'透析患者'（chartevents 225126）"
  at[[4]] <- step(d$esrd_icd,   lab_a, "终末期肾病/透析依赖")
  keep <- keep & !d$esrd_icd
  at[[5]] <- step(d$esrd_chart, lab_b, "透析依赖")
  keep <- keep & !d$esrd_chart

  ## N4 界标时仍在 ICU 且存活（排除原因按方案表 6 拆分两类）
  at[[6]] <- step(d$died_within_24h == 1, "N4 界标时仍在 ICU 且存活 — 24 h 内死亡", "24 h 内死亡")
  keep <- keep & d$died_within_24h == 0
  at[[7]] <- step(d$still_in_icu == 0, "　└ 24 h 内转出 ICU", "24 h 内转出")
  keep <- keep & d$still_in_icu == 1

  ## 排除③ 界标前已开始 RRT
  at[[8]] <- step(d$rrt_pre, "　└ 排除③ 界标前已开始 RRT", "界标前已行 RRT")
  keep <- keep & !d$rrt_pre

  ## N5 0–24h 至少两次血钠
  at[[9]] <- step(!d$has_delta_na, "N5 0–24 h 内 ≥2 次血钠", "血钠不足 2 次")
  keep <- keep & d$has_delta_na

  ## N6 界标前可判定为无 AKI（排除原因拆分："已有 AKI" 与 "无法判定"）
  at[[10]] <- step(d$aki_pre %in% TRUE, "N6 界标前可判定为无 AKI — 已有 AKI", "界标前已满足 KDIGO")
  keep <- keep & !(d$aki_pre %in% TRUE)
  at[[11]] <- step(d$not_judgeable, "　└ 界标前无法判定（肌酐 < 2 次）", "缺乏比较基线")
  keep <- keep & !d$not_judgeable

  ## N7 界标后有可判定的肌酐比较值
  at[[12]] <- step(!d$has_post_obs,
                   "N7 界标后有可判定的肌酐比较值", "界标后无肌酐复测")
  keep <- keep & d$has_post_obs

  ## N3 首次符合条件的 ICU 住院（作用在完全合格集合上，见文件头说明）
  if (first_only) {
    ord <- order(d$subject_id, d$intime)
    first_ok <- rep(FALSE, nrow(d))
    first_ok[ord][!duplicated(d$subject_id[ord])] <- TRUE
    at[[13]] <- step(keep & !first_ok, "　└ 因非首次符合条件的住院而排除", "非首次住院")
    keep <- keep & first_ok
  } else {
    ## S10：不施加首次限制，该节点排除 0 人
    at[[13]] <- data.frame(节点 = "　└ （S10：不做首次限制）", 排除原因 = "—",
                           排除 = 0L, 留存 = n_run, stringsAsFactors = FALSE)
  }

  att <- do.call(rbind, at)

  d$evaluable <- d$evaluable & keep          # N3 之后的最终分析人群
  ## 逻辑下标里只要混进一个 NA，`d[..., ]` 就会**静默插入一行全 NA**，
  ## 分析人群凭空变大而无人报错。宁可在这里硬停。
  if (anyNA(d$evaluable))
    stop(sprintf("evaluable 含 %d 个 NA，无法安全取子集", sum(is.na(d$evaluable))))
  analysis <- d[d$evaluable, , drop = FALSE]

  ## 末两行不是筛选关口，而是结果计数：事件数与分析人群规模
  at_final <- data.frame(
    节点 = c("N8 结局窗内检出 AKI", "最终分析人群"),
    排除原因 = c("—", "—"),
    排除 = c(NA_integer_, NA_integer_),
    留存 = c(sum(analysis$aki), nrow(analysis)),
    stringsAsFactors = FALSE)

  ## 留存列必须一路走到底与分析人群规模吻合，且与 keep 自洽
  stopifnot(sum(keep) == nrow(analysis),
            att$留存[nrow(att)] == nrow(analysis))

  list(data = d, analysis = analysis,
       attrition = rbind(att, at_final),
       sodium = sod, creatinine = cr, na_df = na_df, cr_df = cr_df,
       rrt = rrt, rrt_pre = rrt_pre, rrt_post = rrt_post)
}

#' 把协变量表并进分析人群
#'
#' 协变量在 SQL 侧按 cohort_C 备好（A ⊆ B ⊆ C，一份表覆盖三个口径），此处只按 id
#' 取子集。合并症与治疗类二分类变量缺失一律记 0（"无该诊断/未用该药"），连续变量
#' 保持 NA 不填补——缺失机制与填补策略留到 5.6 一节统一交代。
attach_covariates <- function(analysis, cov) {
  stopifnot("id" %in% names(cov))
  ## **幂等**：已经被并过协变量的人群再次调用时，必须原样返回。
  ## 否则 merge 遇到同名列会生成 ckd.x / ckd.y，而 BASELINE_VARS 与 BASIC_X 找的是
  ## ckd——变量整列"消失"，基线表会从 30 行悄悄缩到 9 行（实测如此），
  ## 模型则直接报"对象未找到"。已落盘的 mimic_B_full.rds 里 analysis 就是并过
  ## 协变量的版本，出表脚本若再并一次必踩此坑。
  dup <- setdiff(intersect(names(cov), names(analysis)), "id")
  if (length(dup)) {
    warning(sprintf("协变量 %s 已在分析人群中，本次合并跳过这些列",
                    paste(dup, collapse = ", ")), call. = FALSE)
    cov <- cov[, setdiff(names(cov), dup), drop = FALSE]
  }
  keep <- setdiff(names(cov), c("stay_id", "hadm_id", "subject_id"))
  ## 全部协变量都已在人群里 → 没有可并的列，直接返回（否则 cov[, "id"] 会退化成
  ## 向量，merge 报错）
  if (!length(keep)) return(analysis)
  d <- merge(analysis, cov[, c("id", keep), drop = FALSE],
             by = "id", all.x = TRUE, sort = FALSE)
  for (v in c("ckd", "diabetes", "hypertension", "chd", "atrial_fib", "valvular",
              "vasoactive_any", "vasopressor_any", "inotrope_any",
              "diuretic_any", "hypertonic_saline",
              "mech_vent_invasive", "mech_vent_noninvasive", "mech_vent_any"))
    if (v %in% names(d)) d[[v]][is.na(d[[v]])] <- 0L
  for (v in c("fluid_in_24h_ml", "urine_out_24h_ml", "fluid_balance_24h_ml",
              "diuretic_furosemide_eq_mg"))
    if (v %in% names(d)) d[[v]][is.na(d[[v]])] <- 0

  ## 血管活性药剂量 nee_max（去甲肾上腺素当量，mcg/kg/min）已在 SQL 侧 COALESCE 为 0，
  ## 未使用者即为 0，不留 NA。**不能留 NA**：glm 默认完整案例删除，若使用者的
  ## 剂量为 NA，扩展模型会悄悄退化成"只对部分用药者"做的分析（此前 vasoactive_max_rate
  ## 版本实测把有效样本从 2,683 压到 725）。
  ##
  ## 一致性断言：有当量剂量者必须确实用过血管加压药。反之不成立且**不应当断言**——
  ## 用了加压药但所有记录都超出临床可信区间（单位/录入错误）者当量恰为 0，
  ## 这是设计使然，由 vasopressor_any 二分类项承担"是否使用"的信息。
  if (all(c("vasopressor_any", "nee_max") %in% names(d))) {
    stopifnot(!any(is.na(d$nee_max)),
              !any(d$nee_max > 0 & d$vasopressor_any != 1))
  }
  d
}

#' 与库无关的别名
#'
#' 本函数**只依赖 ex 的形状，不含任何 MIMIC 特有逻辑**（库差异全部留在提取适配层：
#' SQL 不同、根目录不同，cohort_* 里的 patientunitstayid 已别名为 stay_id）。
#' eICU 侧按这个名字调用，免得读代码的人以为 eICU 的队列走的是另一套判定。
run_cohort <- run_mimic_cohort

#' 打印留存表
print_attrition <- function(x) {
  cat("\n=== 逐层留存（方案表 6 节点）===\n")
  cat(sprintf("  %-46s %10s %10s\n", "节点", "本层排除", "留存"))
  for (i in seq_len(nrow(x$attrition)))
    cat(sprintf("  %-46s %10s %10d\n", x$attrition$节点[i],
                ifelse(is.na(x$attrition$排除[i]), "—", x$attrition$排除[i]),
                x$attrition$留存[i]))
  a <- x$analysis
  cat(sprintf("\n  事件数 %d / 分析人群 %d = %.1f%%\n",
              sum(a$aki), nrow(a), 100 * mean(a$aki)))
  cat(sprintf("  EPV（10 协变量）= %.1f\n", sum(a$aki) / 10))
  invisible(x$attrition)
}
