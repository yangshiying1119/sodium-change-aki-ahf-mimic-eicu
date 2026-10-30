## =============================================================================
## tables.R —— 描述性统计与基线表（方案 5.2 / 4.1.4）
##
## 分组口径按方案 4.1.4 预先固定的三档：ΔNa ≤ −3（下降）、−3 < ΔNa < 3（稳定）、
## ΔNa ≥ 3 mmol/L（升高）。该界限只用于展示与临床解释，**不作为主分析的暴露定义**
## （主分析用连续 delta_na）。
##
## 连续变量的分布型别判定规则在提取前固定：偏度绝对值 > 1 记中位数（四分位数），
## 否则记均值（标准差）。不按组间 P 值挑选呈现方式。
## =============================================================================

NA_GROUP_LEVELS <- c("ΔNa ≤ −3（下降）", "−3 < ΔNa < 3（稳定）", "ΔNa ≥ 3（升高）")

na_group <- function(x) {
  factor(ifelse(is.na(x), NA_character_,
                ifelse(x <= -3, NA_GROUP_LEVELS[1],
                       ifelse(x < 3, NA_GROUP_LEVELS[2], NA_GROUP_LEVELS[3]))),
         levels = NA_GROUP_LEVELS)
}

.skew <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x)
  if (n < 3L) return(0)
  m <- mean(x); s <- sd(x)
  if (!is.finite(s) || s == 0) return(0)
  n / ((n - 1) * (n - 2)) * sum(((x - m) / s)^3)
}

#' 连续变量：偏度 >1 记中位数（P25–P75），否则记均值（标准差）
#'
#' @param style "auto"（按本向量偏度自判）、"median"、"mean"。**分组列必须传
#'   固定的 style**：分布型别是提取前按整体分布定死的呈现规则，若每个亚组各自
#'   判一次，同一行会出现"140.7（4.9）"与"139.0（136.0–141.0）"并排的情况，
#'   读者无从比较。
fmt_cont <- function(x, digits = 1L, style = c("auto", "median", "mean")) {
  style <- match.arg(style)
  x <- suppressWarnings(as.numeric(x)); x <- x[!is.na(x)]
  if (!length(x)) return("—")
  if (style == "auto") style <- if (abs(.skew(x)) > 1) "median" else "mean"
  if (style == "median") {
    q <- quantile(x, c(.25, .5, .75), names = FALSE)
    sprintf(paste0("%.", digits, "f（%.", digits, "f–%.", digits, "f）"),
            q[2], q[1], q[3])
  } else {
    sprintf(paste0("%.", digits, "f（%.", digits, "f）"), mean(x), sd(x))
  }
}

#' 分类变量：n（%）；分母为改组非缺失人数
#'
#' denom 用于**表 8 那一类按列给出组内人数的表**：同一列里每一行的百分比都必须
#' 以该列人数为分母，否则会出现"人数 112 例、AKI 26 例（25.0%）"这种同一列内
#' 分母自相矛盾的写法（26/112 = 23.2%）。组内该变量存在缺失时，缺省的非缺失
#' 分母会静默地换掉读者手上那个分母，故由调用方显式传入。
fmt_cat <- function(x, level = NULL, denom = NULL) {
  x <- x[!is.na(x)]
  n <- length(x); if (!n) return("—")
  if (!is.null(denom) && denom > 0) n <- denom
  if (is.null(level)) {
    if (is.logical(x) || all(x %in% c(0, 1))) {
      k <- sum(as.numeric(x) == 1)
      return(sprintf("%d（%.1f%%）", k, 100 * k / n))
    }
    return(sprintf("%d", n))
  }
  k <- sum(as.character(x) == level)
  sprintf("%d（%.1f%%）", k, 100 * k / n)
}

#' 组间比较：连续用 Kruskal-Wallis，分类用卡方
.p_cont <- function(x, g) {
  ok <- !is.na(x) & !is.na(g)
  if (length(unique(g[ok])) < 2L) return(NA_real_)
  tryCatch(kruskal.test(x[ok], factor(g[ok]))$p.value, error = function(e) NA_real_)
}
.p_cat <- function(x, g) {
  ok <- !is.na(x) & !is.na(g)
  if (length(unique(g[ok])) < 2L) return(NA_real_)
  tryCatch(suppressWarnings(chisq.test(table(x[ok], g[ok]))$p.value),
           error = function(e) NA_real_)
}

## 多水平分类变量的中文标签。
## 不直接印数据库里的英文原串（"Medical/Surgical Intensive Care Unit (MICU/SICU)"
## 这类长串会把表格撑爆且与全文中文体例不符）；未在表内的水平原样保留，
## 保证新出现的水平不会被静默吞掉。
LEVEL_LABELS <- list(
  first_careunit = c(
    "Coronary Care Unit (CCU)"                       = "冠脉监护病房（CCU）",
    "Medical Intensive Care Unit (MICU)"             = "内科监护病房（MICU）",
    "Cardiac Vascular Intensive Care Unit (CVICU)"   = "心血管外科监护病房（CVICU）",
    "Medical/Surgical Intensive Care Unit (MICU/SICU)" = "内外科监护病房（MICU/SICU）",
    "Surgical Intensive Care Unit (SICU)"            = "外科监护病房（SICU）",
    "Trauma SICU (TSICU)"                            = "创伤外科监护病房（TSICU）",
    "Neuro Intermediate"                             = "神经科中间监护",
    "Neuro Surgical Intensive Care Unit (Neuro SICU)" = "神经外科监护病房",
    "Neuro Stepdown"                                 = "神经科降阶梯病房",
    "PACU"                                           = "麻醉后恢复室（PACU）"),
  admission_type = c(
    "EW EMER."                      = "急诊",
    "URGENT"                        = "紧急",
    "OBSERVATION ADMIT"             = "观察入院",
    "DIRECT EMER."                  = "直接急诊入院",
    "SURGICAL SAME DAY ADMISSION"   = "手术当日入院",
    "ELECTIVE"                      = "择期",
    "EU OBSERVATION"                = "急诊观察"))

#' 把一个变量的水平按频数降序排列，并译成中文标签
#'
#' 返回 raw/label 两列：数值统计必须用**原始水平**去匹配（fmt_cat 比的是数据库里的
#' 原值），中文标签只用于显示。直接拿标签去匹配会得到全 0 的计数。
.multi_levels <- function(x, nm) {
  x <- x[!is.na(x)]
  if (!length(x)) return(data.frame(raw = character(0), label = character(0)))
  raw <- names(sort(table(x), decreasing = TRUE))
  lab <- LEVEL_LABELS[[nm]]
  label <- if (is.null(lab)) raw else {
    hit <- lab[raw]
    ifelse(is.na(hit), raw, unname(hit))
  }
  data.frame(raw = raw, label = label, stringsAsFactors = FALSE)
}

#' 基线表的变量清单（名称 → 显示名、类型、分类水平）
#' 类型 "multi"：多水平分类，逐水平各占一行（前面加全角缩进），P 值记在首行。
BASELINE_VARS <- list(
  list("age",             "年龄（岁）",                 "cont"),
  list("gender",          "男性",                       "cat", "M"),
  list("first_careunit",  "ICU 类型",                   "multi"),
  list("admission_type",  "入院类型",                   "multi"),
  list("na_first",        "首次血钠（mmol/L）",          "cont", 1),
  list("na_last",         "末次血钠（mmol/L）",          "cont", 1),
  list("delta_na",        "ΔNa（mmol/L）",              "cont", 1),
  list("n_na",            "血钠检测次数",                "cont", 0),
  list("creat_baseline",  "界标前最低肌酐（mg/dL）",      "cont", 2),
  list("ckd",             "慢性肾病",                    "cat"),
  list("diabetes",        "糖尿病",                      "cat"),
  list("hypertension",    "高血压",                      "cat"),
  list("chd",             "冠心病",                      "cat"),
  list("atrial_fib",      "心房颤动",                    "cat"),
  list("valvular",        "心脏瓣膜病",                  "cat"),
  list("sbp_min_24h",     "0–24h 最低收缩压（mmHg）",     "cont", 0),
  list("hr_max_24h",      "0–24h 最高心率（次/分）",      "cont", 0),
  list("rr_max_24h",      "0–24h 最高呼吸频率（次/分）",   "cont", 0),
  list("gcs_min_24h",     "0–24h 最低 GCS（分）",        "cont", 0),
  list("glucose_max_24h", "0–24h 最高血糖（mg/dL）",      "cont", 0),
  list("hgb_min_24h",     "0–24h 最低血红蛋白（g/dL）",   "cont", 1),
  list("mech_vent_any",   "0–24h 机械通气",              "cat"),
  list("vasoactive_any",  "0–24h 血管活性药",            "cat"),
  list("vasopressor_any", "　└ 其中血管加压药",           "cat"),
  list("inotrope_any",    "　└ 其中强心药",               "cat"),
  ## 血管活性药剂量只在**使用者**中报告：全人群中位数与四分位数恒为 0（未用者记 0），
  ## 整行印成"0.000（0.000–0.000）"没有任何信息量。使用者占比已由上一行的
  ## vasopressor_any 给出，此处只描述"用了的人用了多大剂量"。
  list("nee_max",         "0–24h 去甲肾上腺素当量最大速率，仅使用者（mcg/kg/min）",
       "cont", 3, function(d) !is.na(d$vasopressor_any) & d$vasopressor_any == 1),
  list("diuretic_any",    "0–24h 利尿剂",                "cat"),
  list("diuretic_furosemide_eq_mg", "0–24h 利尿剂（呋塞米当量 mg）", "cont", 0),
  list("fluid_in_24h_ml", "0–24h 液体入量（mL）",        "cont", 0),
  list("urine_out_24h_ml","0–24h 尿量（mL）",            "cont", 0)
)

## 系数名 → 显示名。模型的系数名就是代码里的变量名（delta_na / hgb_min_24h），
## 直接印进论文表读者看不懂，也不该指望审稿人去猜；这里按 BASELINE_VARS 的显示名
## 回译。性别在系数里只出现 M 这一个水平（genderM），单独记一条。
.VAR_LABEL <- local({
  m <- c("(Intercept)" = "常数项", genderM = "男性")
  for (v in BASELINE_VARS) m[[v[[1]]]] <- v[[2]]
  m
})

#' 把系数名换成论文里能读的显示名；查不到的原样返回（宁可露出原名，也不要悄悄改成空）
var_label <- function(nm) {
  hit <- .VAR_LABEL[nm]
  ifelse(is.na(hit), nm, unname(hit))
}

#' 生成基线表（整体 + 三档 ΔNa 分组）
#'
#' @param extra 可选的附加人群列：具名 list(列名 = 逻辑向量)，长度须等于 nrow(d)。
#'   方案 5.2 要求基线表**单列**"结局不可判定者"与"未复测者"的基线特征，
#'   这是判断选择偏倚方向的关键，故不是一个可有可无的装饰列。
#' @return data.frame(变量, 整体, 组1, 组2, 组3, [附加列...], P, 缺失)
baseline_table <- function(d, extra = NULL) {
  g <- na_group(d$delta_na)
  n <- nrow(d)
  for (nm in names(extra))
    stopifnot(is.logical(extra[[nm]]), length(extra[[nm]]) == n)
  rows <- list()
  spec <- list()   # 与 rows 一一对应，记录该行由哪个变量、哪个水平产生
  for (v in BASELINE_VARS) {
    nm <- v[[1]]
    if (!nm %in% names(d)) next
    x <- d[[nm]]
    if (v[[3]] == "cont") {
      dg <- if (length(v) >= 4) v[[4]] else 1L
      ## 可选的子集函数（第 5 个元素）：只在某一人群内描述该变量
      sub <- if (length(v) >= 5 && is.function(v[[5]])) v[[5]](d) else rep(TRUE, n)
      sub[is.na(sub)] <- FALSE
      x <- x[sub]; gs <- g[sub]
      ## 呈现型别按**整体**分布定一次，四个列共用，保证同一行可比
      xv <- suppressWarnings(as.numeric(x)); xv <- xv[!is.na(xv)]
      st <- if (length(xv) && abs(.skew(xv)) > 1) "median" else "mean"
      rows[[length(rows) + 1L]] <- data.frame(
        变量 = v[[2]],
        整体 = fmt_cont(x, dg, st),
        a = fmt_cont(x[gs == NA_GROUP_LEVELS[1]], dg, st),
        b = fmt_cont(x[gs == NA_GROUP_LEVELS[2]], dg, st),
        c = fmt_cont(x[gs == NA_GROUP_LEVELS[3]], dg, st),
        P = .p_cont(as.numeric(x), gs),
        缺失 = sum(is.na(x)),
        row.names = NULL)
      spec[[length(spec) + 1L]] <- list(nm = nm, lv = NULL, type = "cont",
                                        digits = dg, style = st, sub = sub)
    } else if (v[[3]] == "multi") {
      ## 多水平分类：变量名单独占一行（不重复计数），下面逐水平列 n（%）。
      ## 早先按 "cat" 处理时 fmt_cat(level = NULL) 对多水平字符变量只返回总人数，
      ## "ICU 类型"整行印成 "2700"，等于把这一整块信息丢了。
      lvs <- .multi_levels(x, nm)
      pv  <- .p_cat(x, g)
      rows[[length(rows) + 1L]] <- data.frame(
        变量 = paste0(v[[2]], "，n（%）"), 整体 = "", a = "", b = "", c = "",
        P = pv, 缺失 = sum(is.na(x)), row.names = NULL)
      spec[[length(spec) + 1L]] <- list(nm = nm, lv = NULL, type = "head")
      for (k in seq_len(nrow(lvs))) {
        rows[[length(rows) + 1L]] <- data.frame(
          变量 = paste0("　", lvs$label[k]),
          整体 = fmt_cat(x, lvs$raw[k]),
          a = fmt_cat(x[g == NA_GROUP_LEVELS[1]], lvs$raw[k]),
          b = fmt_cat(x[g == NA_GROUP_LEVELS[2]], lvs$raw[k]),
          c = fmt_cat(x[g == NA_GROUP_LEVELS[3]], lvs$raw[k]),
          P = NA_real_, 缺失 = NA_integer_, row.names = NULL)
        spec[[length(spec) + 1L]] <- list(nm = nm, lv = lvs$raw[k], type = "cat")
      }
    } else {
      lv <- if (length(v) >= 4) v[[4]] else NULL
      ## 二分类一律标注 "n（%）"：早先只有指定了水平的（男性）才标，
      ## 慢性肾病/糖尿病等整行只有数没有单位，读者无从判断是人数还是别的
      rows[[length(rows) + 1L]] <- data.frame(
        变量 = paste0(v[[2]], "，n（%）"),
        整体 = fmt_cat(x, lv),
        a = fmt_cat(x[g == NA_GROUP_LEVELS[1]], lv),
        b = fmt_cat(x[g == NA_GROUP_LEVELS[2]], lv),
        c = fmt_cat(x[g == NA_GROUP_LEVELS[3]], lv),
        P = .p_cat(x, g),
        缺失 = sum(is.na(x)),
        row.names = NULL)
      spec[[length(spec) + 1L]] <- list(nm = nm, lv = lv, type = "cat")
    }
  }
  out <- do.call(rbind, rows)
  names(out)[3:5] <- NA_GROUP_LEVELS

  ## ---- 附加人群列（方案 5.2：结局不可判定者 / 未复测者单列） ----
  if (!is.null(extra)) for (k in names(extra)) {
    idx <- extra[[k]]; idx[is.na(idx)] <- FALSE
    vals <- character(length(spec))
    for (j in seq_along(spec)) {
      sp <- spec[[j]]
      v <- Filter(function(z) z[[1]] == sp$nm, BASELINE_VARS)[[1]]
      if (sp$type == "head") { vals[j] <- ""; next }
      ## 连续变量的子集函数（如"仅血管加压药使用者"）在附加列上同样要生效
      if (sp$type == "cont" && !is.null(sp$sub)) idx_j <- idx & sp$sub else idx_j <- idx
      if (sp$type == "cat" && !is.null(sp$lv)) {
        vals[j] <- fmt_cat(d[[sp$nm]][idx_j], sp$lv)
      } else {
        dg <- if (!is.null(sp$digits)) sp$digits else 1L
        vals[j] <- fmt_cont(d[[sp$nm]][idx_j], dg,
                            if (is.null(sp$style)) "auto" else sp$style)
      }
    }
    out[[k]] <- vals
  }
  out
}

#' 打印基线表
print_baseline <- function(tab) {
  cat(sprintf("\n%-34s %-22s %-22s %-22s %-22s %s\n",
              "变量", "整体", NA_GROUP_LEVELS[1], NA_GROUP_LEVELS[2], NA_GROUP_LEVELS[3], "P"))
  for (i in seq_len(nrow(tab)))
    cat(sprintf("%-34s %-22s %-22s %-22s %-22s %s\n",
                tab$变量[i], tab$整体[i], tab[[3]][i], tab[[4]][i], tab[[5]][i],
                ifelse(is.na(tab$P[i]), "—", format.pval(tab$P[i], digits = 2, eps = 1e-4))))
  invisible(tab)
}
