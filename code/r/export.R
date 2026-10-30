## =============================================================================
## export.R —— 把各分析模块的结果对象整理成论文表并落盘
##
## 落盘两份：out/tables/*.csv（UTF-8 带 BOM，Excel 直接双击不乱码）与
## out/tables/表.xlsx（一个工作簿多张表，供撰写时取材与二次排版）。
##
## 设计取舍：**表格内容全部来自已算好的结果对象，本文件不做任何统计**。
## 唯一允许的计算是把数值格式化成"0.933（0.903–0.964）"这类字符串——
## 一旦在这里重算，表与图、正文就可能悄悄对不上。
## =============================================================================

suppressMessages({ library(openxlsx) })

## .fmt_or() / .fmt_p() 定义在 common.R 的"数值呈现"一节：图与表要报同一套格式，
## 只能有一处实现。本文件按惯例在 common.R 之后 source，直接取用即可。

#' 写 CSV（UTF-8 + BOM，Excel 与 R 读入都不乱码）
.write_csv <- function(d, name, dir) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  con <- file(file.path(dir, paste0(name, ".csv")), open = "wb")
  on.exit(close(con), add = TRUE)
  ## 用 con 直接写，fwrite 的 BOM 参数在老版本不可用，这里显式写 BOM 字节
  writeBin(charToRaw("\xEF\xBB\xBF"), con)
  utils::write.csv(d, con, row.names = FALSE, fileEncoding = "UTF-8")
  invisible(file.path(dir, paste0(name, ".csv")))
}

## 工作簿页签名不得含 : \ / ? * [ ]，且 ≤31 字符
.sheet_name <- function(nm) {
  nm <- gsub("[:\\\\/\\?\\*\\[\\]]", "_", nm)
  substr(nm, 1, 31)
}

#' 把一组表写成一个 xlsx 工作簿，每表一页
.write_xlsx <- function(tabs, path) {
  wb <- createWorkbook()
  hs <- createStyle(textDecoration = "bold", halign = "center",
                    border = "bottom", fgFill = "#EFF3F7")
  for (nm in names(tabs)) {
    d <- tabs[[nm]]
    addWorksheet(wb, .sheet_name(nm))
    writeData(wb, .sheet_name(nm), d, headerStyle = hs)
    freezePane(wb, .sheet_name(nm), firstRow = TRUE)
    setColWidths(wb, .sheet_name(nm), cols = seq_len(ncol(d)),
                 widths = "auto")
  }
  saveWorkbook(wb, path, overwrite = TRUE)
  invisible(path)
}

## ---- 各表的构造 --------------------------------------------------------------

#' 表 1 基线特征（整体 + ΔNa 三档）
tbl_baseline <- function(bl) {
  out <- bl[, c("变量", "整体", names(bl)[3:5], "P", "缺失")]
  out$P <- .fmt_p(out$P)
  names(out) <- c("变量", "整体", "ΔNa ≤ −3", "−3 < ΔNa < 3", "ΔNa ≥ 3",
                  "P", "缺失数")
  out
}

#' 表 2 主分析：三档模型下 ΔNa 的效应（其余系数另表）
#'
#' 只保留"OR（95% CI）"这一个呈现串：早先还并列 OR / 下限 / 上限 三个四位
#' 小数原始列，同一行里同一个估计值印了两遍，读者会误以为后面三列是别的模型。
#' 全精度数值仍在 results.json 的 models$tiers 里，核算不受影响。
tbl_main <- function(tiers) {
  d <- tiers[tiers$项 == "delta_na", , drop = FALSE]
  data.frame(
    模型 = d$模型,
    `OR（95% CI）` = .fmt_or(d$OR, d$下限, d$上限),
    P = .fmt_p(d$P), n = d$n, 事件数 = d$事件数,
    check.names = FALSE, row.names = NULL)
}

#' 表 2b 扩展模型完整系数（放附录，供审稿人核算）
tbl_full_coef <- function(tiers) {
  d <- tiers[tiers$模型 == "扩展", , drop = FALSE]
  data.frame(
    变量 = var_label(d$项),
    `OR（95% CI）` = .fmt_or(d$OR, d$下限, d$上限),
    P = .fmt_p(d$P), row.names = NULL,
    ## 全角括号与 % 都不是合法的 R 名称字符，data.frame() 默认的
    ## check.names=TRUE 会把列名悄悄改写成 "OR.95..CI."（写进 CSV 后
    ## 表头就成了这个样子）。必须显式关掉。
    check.names = FALSE)
}

#' 表 3 剂量反应：整体/非线性检验 + 关键节点上的 OR
#'
#' 关键点取 −6 / −4 / −2 / 0 / +2 / +4 / +6，均在曲线定义域内按最近网格点取值；
#' 定义域外的点直接不出现，不做外推。
tbl_rcs <- function(fit, curve, points = c(-6, -4, -2, 0, 2, 4, 6)) {
  head <- data.frame(
    项目 = c("整体关联检验 P", "非线性检验 P", "节点位置（mmol/L）",
             "参考点（mmol/L）"),
    值 = c(.fmt_p(fit$p_overall), .fmt_p(fit$p_nonlinear),
           paste(sprintf("%g", fit$knots), collapse = " / "), "0"),
    row.names = NULL)
  xr <- range(curve$delta_na)
  pts <- points[points >= xr[1] & points <= xr[2]]
  rows <- do.call(rbind, lapply(pts, function(p) {
    ## 参考点单独处理：样条曲线在 ΔNa = 0 处的 OR 按构造恰为 1，
    ## 但网格是固定的 200 点序列，未必正好落在 0 上，取最近点会得到
    ## 0.998（0.996–1.000）这类**看似不是 1** 的值。那不是估计结果，
    ## 只是网格取点造成的显示误差，照抄进表会被读成"参考点 OR 不等于 1"。
    if (abs(p - 0) < 1e-9) {
      return(data.frame(ΔNa = "0（参考）", `OR（95% CI）` = "1.00（参考）",
                        check.names = FALSE, row.names = NULL))
    }
    i <- which.min(abs(curve$delta_na - p))
    data.frame(ΔNa = sprintf("%+g", p),
               `OR（95% CI）` = .fmt_or(curve$OR[i], curve$下限[i], curve$上限[i]),
               check.names = FALSE, row.names = NULL)
  }))
  list(检验 = head, 节点OR = rows)
}

#' 表 4 预设亚组
tbl_subgroup <- function(sg) {
  d <- sg
  ## 分层行的"层"列换成显示名，不要直接把"是/否"抛给读者：光看
  ## "首次血钠 <135 vs ≥135 / 是" 得自己翻译一遍才知道"是"指哪一边，而图 4
  ## 用的是显示名——表与图标签不一致，读者就得在两处之间来回对照。
  ## 映射表 SUBGROUP_LEVELS 与 SUBGROUPS 紧挨着放在 models.R，两处同源。
  d$层 <- ifelse(d$层 == "交互 P", "交互 P", subgroup_labels(d$亚组, d$层))
  d$`OR（95% CI）` <- ifelse(is.na(d$OR) | d$层 == "交互 P", "—",
                             .fmt_or(d$OR, d$下限, d$上限))
  out <- data.frame(
    亚组 = d$亚组, 层 = d$层, n = d$n, 事件数 = d$事件数,
    `OR（95% CI）` = d$`OR（95% CI）`, P = .fmt_p(d$P),
    check.names = FALSE, row.names = NULL)
  ## 交互 P 行的 P 列在源表里就是交互 P，标出来免得读者误当层内 P
  out$P[out$层 == "交互 P"] <- paste0("交互 ", out$P[out$层 == "交互 P"])
  out
}

#' 表 5 敏感性分析
#'
#' @param extra 额外追加的结果行（S8 复合口径、复测概率加权等各自模块算出的行），
#'   列名与主表一致。放在这里统一格式化，避免在各自模块里再写一遍数字格式化。
tbl_sensitivity <- function(s, extra = NULL) {
  if (!is.null(extra) && nrow(extra))
    s <- rbind(s, extra[, names(s), drop = FALSE])
  data.frame(
    分析 = s$分析, n = s$n, 事件数 = s$事件数,
    `OR（95% CI）` = ifelse(is.na(s$OR), s$备注,
                            .fmt_or(s$OR, s$下限, s$上限)),
    P = .fmt_p(s$P), 备注 = s$备注,
    check.names = FALSE, row.names = NULL)
}

#' 表 6 竞争风险：Fine-Gray 次分布风险 + Aalen-Johansen 累积发生率
tbl_competing <- function(cmp) {
  fg <- cmp$finegray
  fg_out <- data.frame(
    变量 = var_label(fg$项),
    `sHR（95% CI）` = .fmt_or(fg$sHR, fg$下限, fg$上限),
    P = .fmt_p(fg$P), row.names = NULL,
    ## 同 tbl_full_coef：不关掉 check.names，表头会被改写成 "sHR.95..CI." 写进 CSV
    check.names = FALSE)
  aj <- cmp$aj_summary
  ## 三行里有一行是"仍在观察"的补集，论文只需 AKI 与竞争事件两条曲线
  aj <- aj[aj$状态 %in% c("AKI", "转出或死亡"), , drop = FALSE]
  ## 累积发生率是 0–1 的比例，而表头与图 5 的纵轴都是百分数；照原样印 0.2154
  ## 会让读者在两处之间换算一次，换算错了就是 21.5 与 0.2 的量级之差。
  names(aj) <- c("状态", "累积发生率（72 h）")
  aj[[2]] <- sprintf("%.1f%%", 100 * aj[[2]])
  list(finegray = fg_out, aj = aj)
}

#' 表 7 逐层留存
tbl_attrition <- function(att) {
  out <- att
  out$排除[is.na(out$排除)] <- "—"
  names(out) <- c("节点", "排除原因", "本层排除", "留存")
  out
}

## ---- 偏倚分析各表（方案 5.6 第四步 / 5.7）--------------------------------------

#' 表 8 可判定结局者与各类未纳入者的基线对比
#'
#' 列由 bias.R 的 .eligibility_class 给定，各类互斥且合起来正好是未纳入者全体，
#' 故不设"合计"列——加了会诱导读者把互斥类别相加。
tbl_bias_compare <- function(cmp) {
  cmp
}

#' 表 9 血钠监测强度分层的复测比例与 AKI 检出率
tbl_bias_monitoring <- function(rm) {
  out <- rm
  names(out) <- c("0–24h 血钠检测次数", "n", "界标后复测比例",
                  "纳入分析比例", "AKI 检出率")
  out
}

#' 表 10 界标后复测肌酐的预测模型（仅界标前协变量）
tbl_retest_model <- function(w) {
  co <- summary(w$model)$coefficients
  data.frame(
    变量 = var_label(rownames(co)),
    `OR（95% CI）` = .fmt_or(exp(co[, 1]),
                             exp(co[, 1] - 1.96 * co[, 2]),
                             exp(co[, 1] + 1.96 * co[, 2])),
    P = .fmt_p(co[, 4]),
    check.names = FALSE, row.names = NULL)
}

#' 一次性导出全部表
#'
#' @param sens_extra 追加到表 5 末尾的结果行（S8 复合口径、复测概率加权）。来自
#'   各自模块（sensitivity.R / bias.R），此处只做拼接与格式化。
#' @param bias 偏倚分析结果（bias.R 的产出）。为 NULL 时表 8–10 不出现，
#'   以便在偏倚分析尚未跑完时也能先出主表。
#' @param extra_tabs 具名列表，追加到表集末尾（如 S8 的肌酐/尿量标准交叉表）。
#'   放在这里而不是导出后再补，是为了让 xlsx 工作簿一次写全。
#' @return 表名 → data.frame 的列表（同时落盘）
export_db <- function(full, models, sens, cmp, sens_extra = NULL, bias = NULL,
                      extra_tabs = NULL,
                      dir = file.path(OUT_DIR, "tables"),
                      xlsx = "表_急性心衰血钠变化与AKI_MIMIC.xlsx") {
  bl   <- baseline_table(attach_covariates(full$analysis, full$covariates))
  rcs  <- tbl_rcs(models$rcs, models$curve, points = c(-6, -4, -2, 0, 2, 4, 6))
  cp   <- tbl_competing(cmp)

  tabs <- list(
    表1_基线特征      = tbl_baseline(bl),
    表2_主分析        = tbl_main(models$tiers),
    表2b_扩展模型系数 = tbl_full_coef(models$tiers),
    表3_剂量反应检验  = rcs$检验,
    表3b_剂量反应节点 = rcs$节点OR,
    表4_预设亚组      = tbl_subgroup(models$subgroups),
    表5_敏感性分析    = tbl_sensitivity(sens, extra = sens_extra),
    表6_竞争风险_FineGray = cp$finegray,
    表6b_累积发生率   = cp$aj,
    表7_逐层留存      = tbl_attrition(full$attrition)
  )

  if (!is.null(bias)) {
    tabs$表8_可判定与未纳入对比 <- tbl_bias_compare(bias$对比)
    tabs$表9_监测强度分层       <- tbl_bias_monitoring(bias$监测)
    tabs$表10_复测概率模型      <- tbl_retest_model(bias$ipw)
  }

  for (nm in names(extra_tabs)) tabs[[nm]] <- extra_tabs[[nm]]

  ## 语言层：OUT_LANG=en 时把表名与全部字符内容换成英文。
  ## 放在最后一步是有意的——上面所有构造仍按中文写，翻译只发生在落盘前，
  ## 于是中英文两版的数值不可能对不上（同一批对象，只换了标签）。
  if (is_en()) {
    tabs <- lapply(tabs, tr_df)
    names(tabs) <- tr_tabname(names(tabs))
  }

  for (nm in names(tabs)) .write_csv(tabs[[nm]], nm, dir)
  .write_xlsx(tabs, file.path(dir, xlsx))
  tabs
}

#' MIMIC-IV 主队列的表（历史入口，中文落盘位置与文件名保持不变，零回归）
#'
#' 英文版写到 tables_en/，与中文目录并列而不是覆盖——中文稿与英文投稿件
#' 常常要同时存在，覆盖式输出会让人分不清手上这份是哪个版本。
export_mimic <- function(...)
  export_db(..., dir = file.path(OUT_DIR, out_subdir("tables", "tables_en")),
            xlsx = out_subdir("表_急性心衰血钠变化与AKI_MIMIC.xlsx",
                              "Table_AHF_Sodium_AKI_MIMIC.xlsx"))

#' eICU-CRD 重复队列的表
#'
#' 与主队列**同一套构造函数**，只有落盘目录与工作簿名不同——两库的"口径一致"
#' 由共用代码保证，不靠人工对齐。本库没有血气钠，S4/S4b 不在敏感性表内出现。
export_eicu <- function(...)
  export_db(..., dir = file.path(OUT_DIR, out_subdir("tables_eicu", "tables_eicu_en")),
            xlsx = out_subdir("表_急性心衰血钠变化与AKI_eICU.xlsx",
                              "Table_AHF_Sodium_AKI_eICU.xlsx"))
