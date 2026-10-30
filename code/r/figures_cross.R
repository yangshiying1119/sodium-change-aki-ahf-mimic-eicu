## =============================================================================
## figures_cross.R —— 双库并排图（方案 6.4 需要"一致性判定"的直观依据）
##
##   图 6  ΔNa 与 AKI 的剂量反应曲线：MIMIC-IV 主队列 与 eICU-CRD 重复队列叠加
##   图 7  两库三档模型的效应量并排森林图
##
## 设计取舍：两库的曲线**不做任何重新标定**，各画各的 OR 与区间——重复队列的
## 意义正在于"同一套口径在另一套数据库上得到什么估计"，叠加前先做平滑或缩放
## 就等于把待检验的东西先抹掉。配色沿用 figures.R：主队列蓝、重复队列赭红。
## =============================================================================

.here_x <- local({
  a <- commandArgs(trailingOnly = FALSE)
  s <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(s) && nzchar(s)) dirname(normalizePath(s)) else "."
})
source(file.path(.here_x, "figures.R"))

COL_MAIN <- "#2166AC"   # MIMIC-IV 主队列
COL_REP  <- "#B2182B"   # eICU-CRD 重复队列

## ---- 图 6 双库剂量反应叠加 ----------------------------------------------------
#' @param cm,ce 两库的 RCS 曲线（delta_na / OR / 下限 / 上限）
#' @param fm,fe 两库的 fit_rcs() 结果（取 P 值标注）
#' @param lab_m,lab_e 图例中的队列名
fig_rcs_cross <- function(cm, ce, fm = NULL, fe = NULL,
                          lab_m = "MIMIC-IV 主队列", lab_e = "eICU-CRD 重复队列",
                          dir = "../out/figures", file = "fig6_rcs_cross.png") {
  path <- .figpath(dir, file)

  draw <- function() {
    par(mar = c(4.8, 4.8, 3.4, 1.4), las = 1)
    yl <- range(c(cm$下限, cm$上限, ce$下限, ce$上限, 1), na.rm = TRUE)
    xl <- range(c(cm$delta_na, ce$delta_na), na.rm = TRUE)
    plot(NA, xlim = xl, ylim = yl, log = "y",
         xlab = if (is_en()) expression(Delta * "Na (0-24 h, mmol/L)")
                else expression(Delta * "Na（0–24 h，mmol/L）"),
         ylab = tr("AKI 的比值比（95% CI）"), main = "")
    .main_text(tr("图 6  早期血钠变化与 AKI 的剂量反应关系：两库对比"))

    abline(h = 1, lty = 2, col = COL_REF, lwd = 1.2)
    abline(v = 0, lty = 3, col = COL_REF, lwd = 1)

    .band <- function(cv, col) {
      polygon(c(cv$delta_na, rev(cv$delta_na)), c(cv$下限, rev(cv$上限)),
              col = adjustcolor(col, 0.15), border = NA)
      lines(cv$delta_na, cv$OR, lwd = 2.4, col = col)
    }
    .band(ce, COL_REP)
    .band(cm, COL_MAIN)

    ## 队列名 + 整体 P 的条目按语言各写一套模板：先翻模板再 sprintf 的话，
    ## 队列名本身也要翻，两处都得管；这里一次写清，避免漏翻括号里的"整体 P"。
    .plab <- function(lab, fit) {
      p <- if (is.null(fit)) "—" else .fmt_p(fit$p_overall)
      if (is_en()) sprintf("%s (overall P = %s)", tr(lab), p)
      else sprintf("%s（整体 P = %s）", lab, p)
    }
    legend("topright", bty = "n", cex = 0.82,
           legend = c(.plab(lab_m, fm), .plab(lab_e, fe), tr("参考 ΔNa = 0")),
           lty = c(1, 1, 2), lwd = c(2.4, 2.4, 1.2),
           col = c(COL_MAIN, COL_REP, COL_REF))
  }

  .fig_render(path, 8.2, 5.8, draw)
  invisible(path)
}

## ---- 图 7 双库并排森林图 ------------------------------------------------------
#' @param tm,te 两库的 models$tiers 数据框（需含 模型 / 项 / OR / 下限 / 上限 / P）
#' @param lab_m,lab_e 两库的名称前缀
fig_tiers_cross <- function(tm, te, lab_m = "MIMIC-IV", lab_e = "eICU-CRD",
                            dir = "../out/figures", file = "fig7_tiers_cross.png") {
  ## 只取 ΔNa 那一行，三档模型各一条
  .sel <- function(ti, lab) {
    ti <- ti[ti$项 == "delta_na", , drop = FALSE]
    ## 保持方案里的档位顺序，而不是数据框里的出现顺序
    ord <- match(c("未调整", "基础", "扩展"), ti$模型)
    ti  <- ti[ord[!is.na(ord)], , drop = FALSE]
    ti$标签 <- if (is_en()) sprintf("%s: %s", tr(lab), tr(ti$模型))
               else sprintf("%s　%s", lab, ti$模型)
    ti
  }
  d <- rbind(.sel(tm, lab_m), .sel(te, lab_e))
  d <- d[!is.na(d$OR), , drop = FALSE]
  n <- nrow(d)
  if (!n) { cat("  没有可画的模型行，跳过图 7\n"); return(invisible(NULL)) }

  ## 与图 3/4 走同一套森林图引擎（forestploter），只是逐行配色区分两库。
  ## 表头写成"模型（ΔNa 每 +1 mmol/L）"而不是只写"模型"：forestploter 的脚注是
  ## 从表格左边缘起排的**单行**文字，左列越窄，脚注越容易横着伸进刻度标签里
  ## （实测"均为 ΔNa 每升高…"那版正好压在 0.90 上）。把暴露单位放进表头，
  ## 既把左列撑宽给脚注腾出位置，也顺势交代了效应量的口径。
  n_m <- sum(grepl(lab_m, d$标签, fixed = TRUE))
  row_col <- c(rep(COL_MAIN, n_m), rep(COL_REP, n - n_m))
  fig_forest(d, "标签", "图 7  两库三档模型的效应量并排",
             note = "两库共用同一套判定与建模代码。",
             dir = dir, file = file,
             header = "模型（ΔNa 每 +1 mmol/L）", row_col = row_col)
}
