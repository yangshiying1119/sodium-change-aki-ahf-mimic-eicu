## =============================================================================
## figures.R —— 论文图（方案 5.7）
##
##   图 1  逐层留存流程图（N0–N8）
##   图 2  ΔNa 与 AKI 的限制性立方样条剂量反应曲线
##   图 3  敏感性分析森林图
##   图 4  预设亚组森林图
##   图 5  Aalen-Johansen 累积发生率曲线
##
## 出图：同一份绘制代码同时出 PNG（400 dpi，供排版与投稿系统预览）与 PDF
## （矢量线图）。多数期刊对线条图要求矢量格式，只给位图会被退稿要求重传。
##
## 中文用 Noto Sans SC。PNG 走 ragg、PDF 走 cairo，两者都经 fontconfig 取字，
## 不需要 showtext/sysfonts 那套手动注册。
##
## 不用 ggplot2 画流程图——boxes-and-arrows 用基础图形更好控制，且流程图必须
## 严格对应 cohort.R 输出的留存表，不接受手工誊抄的数字。
## =============================================================================

suppressMessages({ library(ragg); library(grid); library(forestploter) })

## 界标/结局窗等常数只在 common.R 定义一次，此处沿用，不另立一份
.here_fig <- local({
  a <- commandArgs(trailingOnly = FALSE)
  s <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(s) && nzchar(s)) dirname(normalizePath(s)) else "."
})
source(file.path(.here_fig, "common.R"))
## 亚组的显示标签由模型层定义（SUBGROUP_LEVELS 必须紧挨着 SUBGROUPS 才不会被改漂），
## 图这边只是取用。models.R 本身不依赖提取层，单独 source 是安全的。
source(file.path(.here_fig, "models.R"))

## 中文用 Noto Sans SC；英文版换 Arial 等宽度替代字 Liberation Sans。
## 字体随语言走，不由调用方传——否则总有一处忘了传，出图会是中文字体画英文。
FONT    <- if (is_en()) FONT_EN else "Noto Sans SC"
FIG_DPI <- 400
COL_AKI <- "#B2182B"   # 结局色
COL_CMP <- "#2166AC"   # 竞争事件色
COL_REF <- "grey45"

.figpath <- function(dir, name) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  file.path(dir, name)
}

#' 一份绘制代码，同时出 PNG 与 PDF
#'
#' **不要在每个图函数里各自 dev.off()**：早先的写法是 `.fig_open()` 里
#' `par(xpd = NA)`，为的是让森林图左右留白里的 mtext 不被裁掉。但 xpd=NA 是
#' **整设备**生效的，于是所有 `abline()` / `grid()` 也跟着画到绘图区之外——
#' 图 2 的 OR=1 参考线、图 3/图 7 的竖网格线、图 5 的网格全都贯穿整张画布，
#' 是最容易被审稿人一眼看出的低级瑕疵。
#' 实测（`par(xpd=FALSE)` + 宽 mai）确认 mtext 画在边距里**不依赖 xpd**，
#' 因此这里统一 xpd=FALSE，只有流程图那一个需要越界的标题自行开 xpd。
#'
#' @param draw 无参闭包，里面写该图的绘制代码；PNG 与 PDF 共用它，
#'   所以两种格式不可能出现内容不一致。
.fig_render <- function(path, w, h, draw) {
  agg_png(path, width = w, height = h, units = "in", res = FIG_DPI,
          background = "white")
  par(family = FONT, xpd = FALSE, bg = "white")
  draw()
  dev.off()

  pdf_path <- sub("\\.png$", ".pdf", path)
  if (isTRUE(capabilities("cairo"))) {
    cairo_pdf(pdf_path, width = w, height = h, bg = "white", family = FONT)
  } else {
    ## 没有 cairo 时 pdf() 对 CJK 支持不佳，但至少保证矢量文件存在
    pdf(pdf_path, width = w, height = h, bg = "white")
  }
  par(family = FONT, xpd = FALSE, bg = "white")
  draw()
  dev.off()
  invisible(path)
}

#' 中日韩字符按 2 个西文宽度计，用于估框宽
#'
#' 必须**逐元素**处理并返回等长向量：早先写成 `sum(... strsplit(s,"")[[1]] ...)`，
#' 传向量时 `[[1]]` 只取第一个元素，于是 `max(.cw(标签))` 实际返回的是**第一行**
#' 标签的宽度，森林图左边距据此算出 0.98 in，长标签全被切掉。
.cw <- function(s) {
  vapply(as.character(s), function(x) {
    ch <- strsplit(x, "", fixed = TRUE)[[1]]
    sum(ifelse(grepl("[一-鿿　-〿＀-￯]", ch), 2L, 1L))
  }, numeric(1), USE.NAMES = FALSE)
}

#' 图题按绘图区实际宽度自动缩号
#'
#' 图题是居中排在绘图区上方的，宽过绘图区就被设备直接裁掉——不报错、不留痕。
#' 中文题按"一个汉字占两个西文宽"估算即可，英文题却是同一条信息的两倍长
#' （实测图 2 的题尾 "and AKI" 整段消失）。所以不用估，直接用 strwidth 量：
#' 超了按比例缩 cex。缩过号的题仍比刻度醒目，被裁掉半句则等于没有题。
#' 必须在 plot() 之后调用（strwidth 要有活动绘图区）。
.main_text <- function(txt, cex = 1.0, frac = 0.96) {
  w <- tryCatch(strwidth(txt, units = "inches", cex = cex), error = function(e) NA_real_)
  avail <- par("pin")[1] * frac
  if (is.finite(w) && w > avail) cex <- cex * avail / w
  title(main = txt, cex.main = cex)
  invisible(cex)
}

## ---- 图 1 流程图 --------------------------------------------------------------
#' @param att      cohort.R 的留存表（含 节点 / 排除原因 / 排除 / 留存）
#' @param n_event  结局事件数
#' @param cohort   表型标签
fig_flow <- function(att, n_event, cohort = "B", dir = "../out/figures",
                     file = "fig1_flow.png") {
  path <- .figpath(dir, file)

  draw <- function() {
    ## 唯一需要 xpd 的图：图题放在绘图窗之上，否则被裁掉
    ## mar 必须压到最小：默认边距按"行"算，在 9 英寸高的画布上光上下就吃掉
    ## 1.8 英寸，流程图会缩在中间一条，四周全是白边。
    par(xpd = NA, mar = c(0.5, 0.5, 0.5, 0.5))
    plot.new(); plot.window(c(0, 100), c(0, 100))

    hb  <- 7.0                  # 主链方框高
    xm  <- 26; wm <- 46         # 主链：居中 x 与宽
    xe  <- 76; we <- 46         # 右侧排除框

    box <- function(x, y, w, h, label, sub = NULL, fill = "#EFF3F7",
                    lwd = 1.2, cex = 0.86, border = "grey30") {
      rect(x - w/2, y - h/2, x + w/2, y + h/2, col = fill, border = border, lwd = lwd)
      text(x, y + if (is.null(sub)) 0 else h * 0.17, label, cex = cex)
      if (!is.null(sub)) text(x, y - h * 0.20, sub, cex = cex * 1.10, font = 2)
    }
    ## 排除框的文字用**排除原因**而不是节点名：节点名带"排除①a 诊断码
    ## （N186/5856/…）"这类长串，塞进框里必然溢出；排除原因既短又是流程图惯例。
    ## 第二行同时给"排除 n 例"与"余 m 例"：读者不用翻表 7 就能把 N0 一路减到
    ## 最终人群，逐层留存表与流程图之间的算术闭环在图上就是自明的。
    excl_box <- function(y, reason, n, left, h = 5.2) {
      rect(xe - we/2, y - h/2, xe + we/2, y + h/2,
           col = "#FBEFEF", border = "grey45", lwd = 1)
      text(xe, y + h * 0.17, tr(reason), cex = 0.78, col = "grey15")
      text(xe, y - h * 0.20,
           tr(sprintf("排除 %s 例（余 %s 例）",
                      format(n, big.mark = ","), format(left, big.mark = ","))),
           cex = 0.82, font = 2)
      segments(xm + wm/2, y, xe - we/2, y, lwd = 1, col = "grey45")
    }

    ys <- c(92, 83, 74, 8)      # 三级纳入 + 最终人群的框中心

    ## 先翻译**模板**再 sprintf：模板里的 %s 在词表中原样保留，
    ## 若先 sprintf 再翻译，数字进了串里就匹配不上了。
    labs <- tr(c("N0 数据库内全部 ICU 住院", "N1 成人（≥18 岁）",
                 sprintf(tr("N2 急性心衰表型（定义 %s）"), cohort),
                 "最终分析人群"))
    n_analysis <- att$留存[nrow(att)]
    subs <- tr(c(sprintf("%s 例", format(att$留存[1], big.mark = ",")),
                 sprintf("%s 例", format(att$留存[2], big.mark = ",")),
                 sprintf("%s 例", format(att$留存[3], big.mark = ",")),
                 sprintf("%s 例；结局窗内检出 AKI\n%s 例（%.1f%%）",
                         format(n_analysis, big.mark = ","),
                         format(n_event, big.mark = ","), 100 * n_event / n_analysis)))
    for (i in 1:3) {
      box(xm, ys[i], wm, hb, labs[i], subs[i])
      arrows(xm, ys[i] - hb/2, xm, ys[i + 1] + hb/2,
             length = 0.085, angle = 22, lwd = 1.15, col = "grey25")
    }
    ## 最终人群框与主链同宽：早先写成 wm + 6，左边就顶到画布边上，与上面三个
    ## 框左右不齐，看着像没对齐的错版
    box(xm, ys[4], wm, hb + 3.0, labs[4], subs[4], fill = "#DCE9F5", lwd = 1.8)

    ## N1、N2 两道关口的排除框，画在对应箭头的右侧
    excl_box(ys[1] - hb/2 - (ys[1] - ys[2] - hb)/2, att$排除原因[2], att$排除[2], att$留存[2])
    excl_box(ys[2] - hb/2 - (ys[2] - ys[3] - hb)/2, att$排除原因[3], att$排除[3], att$留存[3])

    ## N2 → 最终分析人群之间的各道排除关口，右侧排成一列。
    ## 末两行（N8 事件数、最终分析人群）不是筛选关口，排除列为 NA，故排除在列。
    tail_rows <- 4:(nrow(att) - 2)
    tail_rows <- tail_rows[!is.na(att$排除[tail_rows]) & att$排除[tail_rows] > 0]

    y_top <- ys[3] - hb/2 - 2.0            # 主链竖线的上端
    y_end <- ys[4] + (hb + 3.0)/2          # 最终人群框的上沿：竖线必须落在这里
    segments(xm, y_top, xm, y_end + 1.6, lwd = 1.15, col = "grey25")
    ## 箭头**向下**指进最终人群框。早先写成 arrows(y_end - 0.5 → y_end)，
    ## 起点在终点之下，箭头反着朝上，正好背离流程方向。
    arrows(xm, y_end + 1.6, xm, y_end, length = 0.085, angle = 22,
           lwd = 1.15, col = "grey25")
    if (length(tail_rows)) {
      yy <- seq(y_top - 1.2, y_end + 1.6, length.out = length(tail_rows))
      h  <- min(5.2, (y_top - y_end) / length(tail_rows) * 0.92)
      for (k in seq_along(tail_rows)) {
        i <- tail_rows[k]
        excl_box(yy[k], att$排除原因[i], att$排除[i], att$留存[i], h = h)
      }
    }
    text(50, 98.5, tr("图 1  研究人群逐层留存"), cex = 1.12, font = 2)
  }

  ## 画布按内容反推：绘图窗已铺满画布（mar≈0），故按方框的实际比例定尺寸即可
  .fig_render(path, 7.3, 9.2, draw)
  invisible(path)
}

## ---- 图 2 RCS 剂量反应 --------------------------------------------------------
#' @param curve  rcs_curve() 的结果（delta_na / OR / 下限 / 上限）
#' @param fit    fit_rcs() 的结果（取 data 画暴露分布、取 P 值标注）
fig_rcs <- function(curve, fit, dir = "../out/figures", file = "fig2_rcs.png") {
  path <- .figpath(dir, file)

  draw <- function() {
    par(mar = c(4.8, 4.8, 3.4, 1.4), las = 1)
    yl <- range(c(curve$下限, curve$上限, 1), na.rm = TRUE)
    plot(curve$delta_na, curve$OR, type = "n", log = "y", ylim = yl,
         xlab = if (is_en()) expression(Delta * "Na (0-24 h, mmol/L)")
                else expression(Delta * "Na（0–24 h，mmol/L）"),
         ylab = tr("AKI 的比值比（95% CI）"), main = "")
    .main_text(tr("图 2  早期血钠变化与 AKI 的剂量反应关系"))
    ## 先画参考线，再压曲线与色带，免得虚线盖在估计上
    abline(h = 1, lty = 2, col = COL_REF, lwd = 1.2)
    abline(v = 0, lty = 3, col = COL_REF, lwd = 1)
    ## 暴露分布（地毯图，示意各处的数据密度）。只画落在横轴范围内者，
    ## 否则 rug 会因超出坐标范围被裁掉并报警告。
    xr <- range(curve$delta_na)
    dx <- fit$data$delta_na[fit$data$delta_na >= xr[1] & fit$data$delta_na <= xr[2]]
    rug(jitter(dx, amount = 0.12), col = "grey72", ticksize = 0.022)
    polygon(c(curve$delta_na, rev(curve$delta_na)),
            c(curve$下限, rev(curve$上限)),
            col = adjustcolor("#2166AC", 0.16), border = NA)
    lines(curve$delta_na, curve$OR, lwd = 2.4, col = "#2166AC")
    legend("topright", bty = "n", cex = 0.84,
           legend = tr(c("OR（95% CI）", "参考 ΔNa = 0",
                         sprintf(tr("整体 P = %s"), .fmt_p(fit$p_overall)),
                         sprintf(tr("非线性 P = %s"), .fmt_p(fit$p_nonlinear)))),
           lty = c(1, 2, NA, NA), lwd = c(2.4, 1.2, NA, NA),
           col = c("#2166AC", COL_REF, NA, NA))
  }

  .fig_render(path, 7.4, 5.6, draw)
  invisible(path)
}

## ---- 谱系森林图（图 3 / 图 4 / 图 7 共用）------------------------------------
#' 用 forestploter 排森林图（https://github.com/adayim/forestploter）
#'
#' 它是 grid/gtable 布局：表头、标签列、数值列、CI 列自然对齐，是现在 R 里做
#' 投稿级森林图的通行做法。早先那版是基础图形手工排的——可靠，但要靠 strwidth
#' 逐次试算左右边距，标签一长就得重调，而且数值列与 CI 列的对齐靠自己拿捏。
#'
#' 两处刻意的取舍：
#'  * **参考行不用 `is_summary`**：那会在该行画一个"合并估计"的菱形，而"主分析"
#'    是基准值不是 pooled 值，菱形会被读成 meta 分析的合并 OR。改为普通方块，
#'    仅把标签与方块标红加粗。
#'  * CI 列的宽度由**空白列的字符数**控制（forestploter 的惯例），不是像素——
#'    所以要按图幅反推，见 `CI_PAD` 与 `.thin_ticks()`。
#'

## CI 列（也就是 x 轴）的物理宽度 = 这一列里空格的数量。刻度能放下几个标签，
## 完全由它决定，所以这里集中定义一次，别在各处再写魔数。
CI_PAD <- 56

#' 抽稀挤在一起的刻度标签
#'
#' CI 列只有几十个空格宽（56 空格 ≈ 1.96 英寸），两位小数的刻度稍密就贴上了。
#' 实测 xlim 只跨 0.12 时，0.94 与 0.96 会连成 "0.940.96"，比没有刻度还难读。
#' 这里按"刻度标签占轴长的比例"估算所需最小间距，再从左往右贪心保留。
#'
#' @param at     候选刻度（pretty() 的结果）
#' @param xlim   x 轴范围（对数刻度）
#' @param digits 标签小数位
.thin_ticks <- function(at, xlim, digits) {
  if (length(at) < 2) return(at)
  ## 轴宽 ≈ CI_PAD 个空格 × 0.28 em；一个数字 ≈ 0.56 em（无衬线字体的常见取值）。
  ## 数字个数取 digits + 2（小数点 + 整数位的 0）。
  axis_pt <- CI_PAD * 0.28 * 9
  lab_pt  <- (digits + 2) * 0.56 * 9 * 0.85
  need    <- 1.3 * lab_pt / axis_pt * diff(log(xlim))   # 对数轴，间距按比值算
  keep <- at[1]
  for (v in at[-1])
    if (abs(log(v) - log(keep[length(keep)])) >= need) keep <- c(keep, v)
  keep
}

#' @param tab     含 标签列、OR、下限、上限、P 的数据框
#' @param lbl     标签列名
#' @param title   图题
#' @param note    图下脚注
#' @param header  第一列的表头文字
#' @param ref_row 需标红加粗的标签值（通常是"主分析"）
#' @param row_col 可选，逐行的 CI 方块与标签配色；给 NULL 则统一蓝色、
#'   参考行标红。两库并排（图 7）靠它区分主队列与重复队列。
fig_forest <- function(tab, lbl, title, note = NULL, xlim = NULL,
                       dir = "../out/figures", file = "fig_forest.png",
                       ref_row = NULL, header = "分析", row_col = NULL) {
  path <- .figpath(dir, file)
  d <- tab[!is.na(tab$OR), , drop = FALSE]
  n <- nrow(d)
  if (!n) { cat("  无可画的模型行，跳过", file, "\n"); return(invisible(NULL)) }

  is_ref <- if (is.null(ref_row)) rep(FALSE, n) else as.character(d[[lbl]]) %in% ref_row
  ## 配色：默认全部蓝、参考行红；给了 row_col 就按行用（图 7 用它分两库）
  base_col <- if (is.null(row_col)) rep("#2166AC", n) else rep_len(row_col, n)
  if (is.null(row_col) && any(is_ref)) base_col[is_ref] <- "#B2182B"

  if (is.null(xlim)) {
    r <- range(c(d$下限, d$上限, 1), na.rm = TRUE)
    xlim <- c(max(0.3, r[1] - 0.04 * diff(r)), min(2.2, r[2] + 0.04 * diff(r)))
  }
  ## 刻度按数据范围现算，且**最多 4 个**：CI 列宽度由下面空白列的字符数定死，
  ## 刻度一密标签就互相压在一起——实测 5 个两位小数刻度（0.90 0.95 1.00 1.05
  ## 1.10）在同一行糊成一团，比没有刻度还糟。
  at <- pretty(xlim, 5)
  at <- at[at >= xlim[1] & at <= xlim[2]]
  if (!length(at)) at <- mean(xlim)
  if (length(at) > 4) at <- at[round(seq(1, length(at), length.out = 4))]
  digits <- if (diff(range(at)) < 0.3) 2 else 1
  at <- .thin_ticks(at, xlim, digits)

  ## 标签先翻译再算宽度：英文串的字符数比中文多（"S1 最早两次血钠之差" 9 字 vs
  ## "S1 ΔNa from the first two sodium measurements" 46 字符），拿中文原串估宽会
  ## 把英文标签裁掉半截——实测图 3 的左列只剩 "he first two sodium measurements"。
  labs_disp <- tr(as.character(d[[lbl]]))
  hdr_disp  <- tr(header)

  ## 画布尺寸按内容反推：forestploter 把整张表**居中**放，画布给大了就是四周
  ## 一圈白边，投出去显得图很虚。宽度由左列（标签或表头，取宽者）决定，其余各列
  ## 的宽度基本固定，故用常数补足。中文按 2 个西文宽度计。
  w <- max(5.6, 0.055 * max(.cw(labs_disp), .cw(hdr_disp)) + 5.0)
  h <- max(2.4, 0.30 * n + 1.60)

  show <- data.frame(
    a = labs_disp,
    b = paste(rep(" ", CI_PAD), collapse = ""),
    ## 走 .fmt_or() 而不是就地 sprintf：括号的全/半角要跟表 5 一致，
    ## 图的 CI 列与表的"OR（95% CI）"列必须是同一串字符
    c = .fmt_or(d$OR, d$下限, d$上限),
    ## 走 .fmt_p() 而不是 format.pval()：后者整列一起格式化时，会按列里最小的
    ## 那个 P 对齐小数位，一个 0.00067 就把整列撑成 5 位小数（"0.24281"），
    ## 与表 2/表 5 的 4 位小数对不上。
    e = .fmt_p(d$P),
    check.names = FALSE, stringsAsFactors = FALSE)
  names(show) <- c(hdr_disp, " ", tr("OR（95% CI）"), "P")

  draw <- function() {
    sty <- forest_style(
      base_size = 9, base_family = FONT,
      ci = gpar(col = "#2166AC", fill = "#2166AC", lwd = 1.4),
      ci_pch = 15, ci_t_height = 0.18,
      ref_line = gpar(col = COL_REF, lty = 2, lwd = 1.2),
      header = gpar(fontfamily = FONT, cex = 0.88, fontface = "bold"),
      xlab = gpar(fontfamily = FONT, cex = 0.88, col = "grey20"),
      title = gpar(fontfamily = FONT, cex = 1.05, fontface = "bold"),
      title_just = "center",
      footnote = gpar(fontfamily = FONT, cex = 0.70, col = "grey35"))

    p <- forest(show, est = d$OR, lower = d$下限, upper = d$上限,
                ci_column = 2, ref_line = 1, sizes = 0.72, style = sty)
    p <- set_xaxis(p, xlim = xlim, ticks_at = at, ticks_digits = digits)
    p <- set_labs(p, xlab = tr("比值比（OR，对数刻度）"),
                  title = tr(title), footnote = if (is.null(note)) NULL else tr(note))

    ## CI 方块按行配色（图 7 靠它区分两库；其余图统一蓝、参考行红）。
    ## **必须显式 which = "ci"**：edit_plot 的 which 缺省是 "text"，那样会去改
    ## 正文文字格（而且 col 会被内部 +1、row 会被加上表头偏移），CI 图元本身
    ## 一动不动——实测表现为"整张图全是蓝的"，且不报任何错。
    ## which = "ci" 时 row/col 均按原值使用，不再偏移，故 row 就是正文行号。
    for (cc in unique(base_col)) {
      p <- edit_plot(p, row = which(base_col == cc), col = 2, part = "body",
                     which = "ci", gp = gpar(col = cc, fill = cc))
    }
    if (any(is_ref)) {
      p <- edit_plot(p, row = which(is_ref), col = 1, which = "text",
                     gp = gpar(col = "#B2182B", fontface = "bold", fontfamily = FONT))
      p <- edit_plot(p, row = which(is_ref), col = 3:4, which = "text",
                     gp = gpar(fontface = "bold", fontfamily = FONT))
    }
    grid.draw(p)
  }

  .fig_render(path, w, h, draw)
  invisible(path)
}

## ---- 图 5 Aalen-Johansen 累积发生率 -------------------------------------------
#' @param d  competing.R::build_competing() 的输出
fig_aj <- function(d, dir = "../out/figures", file = "fig5_aj.png") {
  path <- .figpath(dir, file)
  suppressMessages(library(survival))

  ## 水平名同时是 survfit 的 pstate 列名（下面按 "转出或死亡" 取列），
  ## 所以**不能在因子层面翻译**——只在图例等展示处 tr()。
  st <- factor(d$status, levels = c(0, 1, 2),
               labels = c("censored", "AKI", "转出或死亡"))
  fit <- survfit(Surv(d$time, st) ~ 1, data = d)
  ts <- seq(LANDMARK_H, HORIZON_H, by = 1)
  s <- summary(fit, times = ts, extend = TRUE)

  draw <- function() {
    ## 下边距要同时容下刻度、轴标题、脚注三行：
    ## 早先 mar 底部 4.6 且轴标题用 plot() 默认 line（=3），脚注又画在 line=3.3，
    ## 两者**重叠**打印在同一行——"ICU 入科后小时数"与"风险集 n = …"糊成一团。
    ## 现在轴标题改显式 mtext(line=2.0)，脚注退到 line=3.9，各自独占一行。
    par(mar = c(5.6, 4.8, 3.4, 1.4), las = 1)
    ymax <- max(s$pstate, na.rm = TRUE) * 1.12
    plot(ts, s$pstate[, "AKI"], type = "n", ylim = c(0, ymax),
         xlab = "", ylab = tr("累积发生率"), main = "")
    .main_text(tr("图 5  界标后 AKI 与竞争事件的累积发生率"))
    grid(col = "grey90", lty = 1)
    ## 三条曲线：仍在观察（补到 1）、AKI、转出或死亡
    for (nm in c("AKI", "转出或死亡")) {
      col <- if (nm == "AKI") COL_AKI else COL_CMP
      lines(ts, s$pstate[, nm], lwd = 2.4, col = col)
    }
    abline(v = LANDMARK_H, lty = 3, col = COL_REF)
    mtext(tr("ICU 入科后小时数"), side = 1, line = 2.2, cex = 0.92)
    legend("topleft", bty = "n", lwd = 2.4, cex = 0.88,
           col = c(COL_AKI, COL_CMP),
           legend = tr(c("AKI（结局）", "转出或死亡（竞争事件）")))
    ## 人数带千分位：图 1 的方框写 "94,458"，这里若写 "2656" 就是同一篇图里两种写法
    mtext(tr(sprintf("风险集 n = %s；界标后 72 h 内 AKI %s 例，转出或死亡 %s 例",
                     format(nrow(d), big.mark = ","),
                     format(sum(d$status == 1L), big.mark = ","),
                     format(sum(d$status == 2L), big.mark = ","))),
          side = 1, line = 3.9, cex = 0.74, col = "grey30")
  }

  .fig_render(path, 7.4, 5.4, draw)
  invisible(path)
}
