## 把全部已落盘的结果汇总成一份 results.json，供成稿脚本（build_manuscript.py）取数。
##
## 为什么要多这一步：论文正文里的每一个数字都应当**自动来自算好的结果对象**，
## 而不是从终端输出里手抄。手抄一次就可能抄错一次，而且改了分析之后正文不会跟着变。
## 本文件只做"读数 → 落 JSON"，不做任何统计。
##
## 双库结构：pack() 按库前缀打包，两库走的是**同一个函数**。方案 6.4 把"跨库口径
## 漂移"列为主要风险，对策是两库共用同一套判定与建模代码；汇总层同样共用，
## 差异只体现在读入的文件名前缀上。MIMIC 的键位保持在顶层（成稿脚本据此取数，
## 不动即为零回归），eICU 的结果整体放在 o$eICU 之下。
source("extract.R")
source("common.R")
source("cohort.R")
source("models.R")

rd <- function(f) if (file.exists(file.path(OUT_DIR, f)))
  readRDS(file.path(OUT_DIR, f)) else NULL

## eICU 侧用哪个表型口径作为重复队列，在此一处指定。
## eICU 的编码习惯以 I50.9 / 428.0（未特指）为主，急性子码全库仅 818 个单元，
## 故方案原定的表型 B 在本库只剩 210 例 / 44 事件（EPV 4.4），撑不起任何调整模型；
## 表型 C 得 2,812 例 / 569 事件（EPV 56.9），与 MIMIC 主队列（2,700 / 718）规模相当。
EICU_CC <- "C"

## 文件名先按 "<库>_<口径>_x.rds" 找，找不到再退回 "<库>_x.rds"——
## MIMIC 侧历史文件名不带口径字母，靠这条回退保持零回归。
rd_any <- function(...) {
  for (f in c(...)) {
    p <- file.path(OUT_DIR, f)
    if (file.exists(p)) return(readRDS(p))
  }
  NULL
}

## 数字一律四舍五入到论文会写的位数；NA 落成 null 由 jsonlite 处理
r4 <- function(x, k = 4) if (is.null(x) || !length(x) || is.na(x)) NULL else round(as.numeric(x), k)
r2 <- function(x) r4(x, 2)
r1 <- function(x) r4(x, 1)

## ---------------------------------------------------------------------------
## pack()：把一个库的全部结果对象打成一层 list。前缀 p 取 "mimic" / "eicu"。
## ---------------------------------------------------------------------------
pack <- function(p, cc = "B") {
  full <- rd_any(sprintf("%s_%s_full.rds", p, cc), sprintf("%s_full.rds", p))
  mods <- rd_any(sprintf("%s_%s_models.rds", p, cc), sprintf("%s_models.rds", p))
  sens <- rd_any(sprintf("%s_%s_sensitivity.rds", p, cc), sprintf("%s_sensitivity.rds", p))
  cmp  <- rd_any(sprintf("%s_%s_competing.rds", p, cc), sprintf("%s_competing.rds", p))
  s8   <- rd_any(sprintf("%s_%s_s8_uo.rds", p, cc), sprintf("%s_s8_uo.rds", p))
  bias <- rd_any(sprintf("%s_%s_bias.rds", p, cc), sprintf("%s_bias.rds", p))
  bcm  <- rd_any(sprintf("%s_%s_bias_compare.rds", p, cc), sprintf("%s_bias_compare.rds", p))
  disc <- rd_any(sprintf("%s_%s_discrim.rds", p, cc), sprintf("%s_discrim.rds", p))
  if (is.null(full) && is.null(mods) && is.null(sens)) return(NULL)

  o <- list()

  ## ---- 队列 ----
  if (!is.null(full)) {
    a <- full$analysis
    keep <- model_cases(a)
    o$队列 <- list(
      逐层留存 = lapply(seq_len(nrow(full$attrition)), function(i) list(
        节点 = full$attrition$节点[i],
        排除原因 = if (is.na(full$attrition$排除原因[i])) NULL else full$attrition$排除原因[i],
        本层排除 = if (is.na(full$attrition$排除[i])) NULL else full$attrition$排除[i],
        留存 = full$attrition$留存[i])),
      分析人群 = nrow(a),
      事件数 = sum(a$aki),
      事件率 = r4(100 * mean(a$aki), 2),
      模型人群 = sum(keep),
      ## 完整病例集里的事件数：摘要要写"其中多少例进入调整模型"，只给人数不够，
      ## 事件数必须同时在手，否则只能由全样本事件数反推，那是猜的。
      事件数_模型 = sum(a$aki[keep]),
      因协变量缺失剔除 = sum(!keep),
      EPV = r2(sum(a$aki) / length(.pick(BASIC_X, a))),
      EPV_扩展 = r2(sum(a$aki) / length(.pick(c(BASIC_X, EXTEND_X), a))))
  }

  ## ---- 主模型三档 ----
  if (!is.null(mods)) {
    ti <- mods$tiers
    o$主模型 <- lapply(unique(ti$模型), function(nm) {
      r <- ti[ti$模型 == nm & ti$项 == "delta_na", ]
      list(模型 = nm, n = r$n, 事件数 = r$事件数, OR = r4(r$OR, 3),
           下限 = r4(r$下限, 3), 上限 = r4(r$上限, 3), P = r4(r$P, 6))
    })
    o$扩展模型系数 <- lapply(seq_len(nrow(ti)), function(i)
      if (ti$模型[i] == "扩展") list(变量 = ti$项[i], OR = r4(ti$OR[i], 3),
                                    下限 = r4(ti$下限[i], 3), 上限 = r4(ti$上限[i], 3),
                                    P = r4(ti$P[i], 4)) else NULL)
    o$扩展模型系数 <- Filter(Negate(is.null), o$扩展模型系数)

    cv <- mods$curve
    o$剂量反应 <- list(
      整体P = r4(mods$rcs$p_overall, 5), 非线性P = r4(mods$rcs$p_nonlinear, 5),
      节点 = round(mods$rcs$knots, 2),
      曲线 = lapply(seq(1, nrow(cv), length.out = 0), identity))
    ## 关键节点上的 OR（取曲线上的最近点，不外推）
    pts <- c(-6, -3, 0, 3, 6)
    pts <- pts[pts >= min(cv$delta_na) & pts <= max(cv$delta_na)]
    o$剂量反应$节点OR <- lapply(pts, function(x) {
      i <- which.min(abs(cv$delta_na - x))
      list(ΔNa = r1(cv$delta_na[i]), OR = r4(cv$OR[i], 3),
           下限 = r4(cv$下限[i], 3), 上限 = r4(cv$上限[i], 3))
    })

    sg <- mods$subgroups
    o$亚组 <- lapply(seq_len(nrow(sg)), function(i)
      list(亚组 = sg$亚组[i], 层 = sg$层[i], n = sg$n[i], 事件数 = sg$事件数[i],
           OR = r4(sg$OR[i], 3), 下限 = r4(sg$下限[i], 3),
           上限 = r4(sg$上限[i], 3), P = r4(sg$P[i], 4)))
  }

  ## ---- 敏感性分析 ----
  if (!is.null(sens)) {
    o$敏感性 <- lapply(seq_len(nrow(sens)), function(i)
      list(分析 = sens$分析[i], n = sens$n[i], 事件数 = sens$事件数[i],
           OR = r4(sens$OR[i], 3), 下限 = r4(sens$下限[i], 3),
           上限 = r4(sens$上限[i], 3), P = r4(sens$P[i], 6),
           备注 = if (is.na(sens$备注[i])) NULL else sens$备注[i]))
  }

  ## ---- 竞争风险 ----
  if (!is.null(cmp)) {
    fg <- cmp$finegray
    r  <- fg[fg$项 == "delta_na", ]
    o$竞争风险 <- list(
      sHR = r4(r$sHR, 3), 下限 = r4(r$下限, 3), 上限 = r4(r$上限, 3), P = r4(r$P, 6),
      累积发生率 = lapply(seq_len(nrow(cmp$aj_summary)), function(i)
        list(状态 = cmp$aj_summary$状态[i],
             值 = r4(cmp$aj_summary$累积发生率[i], 4))))
  }

  ## ---- S8 尿量复合口径 ----
  if (!is.null(s8)) {
    r <- s8$行
    ## 复合口径**同时**改动了"界标前已满足 KDIGO"的排除与结局判定，故与主分析
    ## 不是同一人群上的比较。把该口径下的队列规模与 N6 排除数一并带出，
    ## 成稿时才不必写死数字（写死就会随流水线变动而与表 7 打架）。
    .n6 <- function(rr) {
      a <- rr$attrition
      i <- grep("N6", a$节点)
      if (!length(i)) return(NA_integer_)
      as.integer(a$排除[i[1]])
    }
    o$S8 <- list(n = r$n, 事件数 = r$事件数, OR = r4(r$OR, 3),
                 下限 = r4(r$下限, 3), 上限 = r4(r$上限, 3), P = r4(r$P, 6),
                 队列人数 = nrow(s8$队列$analysis),
                 界标前已满足KDIGO = .n6(s8$队列),
                 交叉 = lapply(seq_len(nrow(s8$交叉)), function(i)
                   as.list(s8$交叉[i, , drop = FALSE])),
                 仅尿量检出 = nrow(s8$仅尿量))
  }

  ## ---- 偏倚分析 ----
  if (!is.null(bcm)) {
    tab <- table(bcm$分组)
    o$可判定分组 <- lapply(names(tab), function(k) list(类别 = k, n = as.integer(tab[[k]])))
  }
  if (!is.null(bias)) {
    o$复测加权 <- list(OR = r4(bias$模型$OR, 3), 下限 = r4(bias$模型$下限, 3),
                       上限 = r4(bias$模型$上限, 3), P = r4(bias$模型$P, 6),
                       n = bias$模型$n, 事件数 = bias$模型$事件数,
                       备注 = bias$模型$备注,
                       截断下 = r4(bias$ipw$截断点[1], 3),
                       截断上 = r4(bias$ipw$截断点[2], 3))
    o$复测模型 <- lapply(seq_len(nrow(bias$监测)), function(i)
      as.list(bias$监测[i, , drop = FALSE]))
  }

  ## ---- 判别力、校准与临床量级（方案 5.3）----
  if (!is.null(disc)) {
    o$判别与校准 <- list(
      n = disc$n, 事件数 = disc$事件数,
      每3mmol = list(OR = r4(disc$每3mmol$OR, 3),
                     下限 = r4(disc$每3mmol$下限, 3),
                     上限 = r4(disc$每3mmol$上限, 3),
                     P = r4(disc$每3mmol$P, 6)),
      LRT = list(卡方 = r2(disc$LRT$卡方), P = r4(disc$LRT$P, 6)),
      HL = list(卡方 = r2(disc$HL$卡方), 自由度 = disc$HL$自由度,
                P = r4(disc$HL$P, 4)),
      判别 = list(AUC = r4(disc$判别$AUC, 3),
                  AUC下限 = r4(disc$判别$AUC下限, 3),
                  AUC上限 = r4(disc$判别$AUC上限, 3),
                  AUC无暴露 = r4(disc$判别$AUC无暴露, 3),
                  差值 = r4(disc$判别$差值, 4),
                  差值P = r4(disc$判别$差值P, 4)),
      校正 = list(c校正 = r4(disc$校正$c校正, 3),
                  c校正下限 = r4(disc$校正$c校正下限, 3),
                  c校正上限 = r4(disc$校正$c校正上限, 3),
                  乐观量 = r4(disc$校正$乐观量, 4),
                  斜率校正 = r4(disc$校正$斜率校正, 3),
                  B = disc$校正$B),
      校准曲线 = lapply(seq_len(nrow(disc$校准曲线)), function(i) list(
        组 = as.character(disc$校准曲线$组[i]),
        n = disc$校准曲线$n[i],
        平均预测 = r4(disc$校准曲线$平均预测[i], 4),
        实测发生率 = r4(disc$校准曲线$实测发生率[i], 4))),
      标准化概率 = lapply(seq_len(nrow(disc$标准化概率)), function(i)
        as.list(disc$标准化概率[i, , drop = FALSE])),
      层内人数 = as.list(disc$层内人数))
  }

  o
}

o <- list(生成时间 = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))

## MIMIC-IV 主队列：键位直接铺在顶层，成稿脚本按既有路径取数
mim <- pack("mimic")
if (!is.null(mim)) o <- c(o, mim)

## eICU-CRD 重复队列：整体置于 o$eICU 之下
eic <- pack("eicu", EICU_CC)
if (!is.null(eic)) {
  eic$表型口径 <- EICU_CC
  o$eICU <- eic
}

## eICU 侧按方案**原定**表型的对照结果。该口径在该库只剩 210 例 / 44 事件
## （EPV 4.4），不足以支撑调整模型，故不作重复队列主分析，但必须如实报告——
## 只报改用宽口径之后的数字、不交代原口径的结果，等于替读者做了取舍。
eicB <- pack("eicu", "B")
if (!is.null(eicB)) {
  eicB$表型口径 <- "B"
  o$eICU严格表型 <- eicB
}

## ---------------------------------------------------------------------------
## 补充分析：严重度梯度、绝对风险差、E 值（_run11.R）与阴性对照暴露（_run13.R）。
##
## 这两块不由 pack() 生成，因为它们的 rds 来自独立驱动，且不随表型字母以外的
## 口径变化。缺文件时整块缺席——成稿脚本按 g(..., default=None) 取数，相应段落
## 自动跳过，不会写半句话。
## ---------------------------------------------------------------------------
补充 <- function(prefix, cc) {
  adv <- rd_any(sprintf("%s_%s_advance.rds", prefix, cc),
                sprintf("%s_advance.rds", prefix))
  neg <- rd_any(sprintf("%s_%s_negctl.rds", prefix, cc),
                sprintf("%s_negctl.rds", prefix))
  out <- list()
  if (!is.null(adv)) {
    out$严重度 <- lapply(seq_len(nrow(adv$严重度)), function(i)
      as.list(adv$严重度[i, , drop = FALSE]))
    ## 保留 6 位小数：正文把比例乘 100 后取 1 位小数，若这里只留 4 位，
    ## 边界值（如 −0.041492）会被提前截成 −0.0415，正文写成 −4.2 而表里是 −4.1。
    out$风险 <- lapply(adv$风险, function(r) list(
      从     = r4(r$from, 6),  到     = r4(r$to, 6),
      从下限 = r4(r$from_ci[1], 6), 从上限 = r4(r$from_ci[2], 6),
      到下限 = r4(r$to_ci[1], 6),   到上限 = r4(r$to_ci[2], 6),
      差值   = r4(r$diff, 6),
      差值下限 = r4(r$diff_ci[1], 6), 差值上限 = r4(r$diff_ci[2], 6)))
    out$基础效应 <- list(OR = r4(adv$基础效应$OR, 3),
                         下限 = r4(adv$基础效应$下限, 3),
                         上限 = r4(adv$基础效应$上限, 3),
                         P = r4(adv$基础效应$P, 6))
    out$E值 <- lapply(seq_len(nrow(adv$E值)), function(i)
      as.list(adv$E值[i, , drop = FALSE]))
    out$基线风险 <- r4(adv$基线风险, 6)
    out$n <- adv$n; out$事件数 <- adv$事件数
  }
  if (!is.null(neg)) {
    ## 结果对象里的暴露是列名（delta_na 等）；成稿脚本按 "ΔNa"/"ΔK"/"ΔCl"
    ## 取行，这里必须换掉，否则 3.10 节的段落会静默整段消失。
    ne <- neg$效应
    disp <- c(delta_na = "ΔNa", delta_k = "ΔK", delta_cl = "ΔCl",
              "delta_na（同时调整 delta_k）" = "ΔNa（同时调整 ΔK）",
              "delta_na（同时调整 delta_cl）" = "ΔNa（同时调整 ΔCl）",
              "delta_cl（同时调整 delta_na）" = "ΔCl（同时调整 ΔNa）",
              "delta_na（每 1 SD，调整 delta_cl）" = "ΔNa（每 1 SD，同时调整 ΔCl）",
              "delta_cl（每 1 SD，调整 delta_na）" = "ΔCl（每 1 SD，同时调整 ΔNa）")
    ne$暴露 <- unname(disp[ne$暴露])
    out$阴性对照 <- lapply(seq_len(nrow(ne)), function(i)
      as.list(ne[i, , drop = FALSE]))
    ## 只保留上三角（ΔNa–ΔK、ΔNa–ΔCl、ΔK–ΔCl），避免正文引用时还要看对称位
    rho <- neg$相关
    out$相关 <- list(Na_K = r4(rho["delta_na", "delta_k"], 3),
                     Na_Cl = r4(rho["delta_na", "delta_cl"], 3),
                     K_Cl = r4(rho["delta_k", "delta_cl"], 3))
    out$非缺失 <- as.list(neg$非缺失)
    ## 互调模型的共线性与标准化尺度：正文用它解释"每 1 SD"一行怎么来的
    if (!is.null(neg$共线)) out$共线 <- list(
      r = r4(neg$共线$r, 3), VIF = r4(neg$共线$VIF, 2),
      SD_Na = r4(neg$共线$SD_Na, 2), SD_Cl = r4(neg$共线$SD_Cl, 2))
  }
  if (length(out)) out else NULL
}

.supp <- 补充("mimic", "B")
if (!is.null(.supp)) o$补充分析 <- .supp
if (!is.null(o$eICU)) {
  .supp_e <- 补充("eicu", EICU_CC)
  if (!is.null(.supp_e)) o$eICU$补充分析 <- .supp_e
}

## ---------------------------------------------------------------------------
## 双库对照：方案 6.4 的一致性判定需要把两库的主估计并排看，这里先把并排所需
## 的字段抽出来。方向一致性与区间重叠与否的**判定**留给成稿脚本按 6.4 的措辞写，
## 本文件只给原始估计，不在数据层替结论定性。
## ---------------------------------------------------------------------------
if (!is.null(o$主模型) && !is.null(eic$主模型)) {
  .主 <- function(blk, nm) {
    for (x in blk) if (identical(x$模型, nm)) return(x)
    NULL
  }
  o$双库对照 <- lapply(c("未调整", "基础", "扩展"), function(nm) {
    m <- .主(o$主模型, nm); e <- .主(eic$主模型, nm)
    list(模型 = nm,
         MIMIC_n = if (is.null(m)) NULL else m$n,
         MIMIC_OR = if (is.null(m)) NULL else m$OR,
         MIMIC_下限 = if (is.null(m)) NULL else m$下限,
         MIMIC_上限 = if (is.null(m)) NULL else m$上限,
         eICU_n = if (is.null(e)) NULL else e$n,
         eICU_OR = if (is.null(e)) NULL else e$OR,
         eICU_下限 = if (is.null(e)) NULL else e$下限,
         eICU_上限 = if (is.null(e)) NULL else e$上限)
  })
}

## ---------------------------------------------------------------------------
## 跨库汇总（方案 6.4）：两库**不做个体数据层面的合并**。若须给出总体估计，只能
## 对两库已得的效应量作随机效应汇总，并在文中标明该汇总为二次计算。
## 这里按 DerSimonian–Laird 法在 log OR 尺度上汇总基础模型，只用 base R，
## 不引入 metafor 依赖（R 4.1.2 装新版包有版本约束，能不加依赖就不加）。
## ---------------------------------------------------------------------------
if (!is.null(o$主模型) && !is.null(eic$主模型)) {
  m0 <- .主(o$主模型, "基础"); e0 <- .主(eic$主模型, "基础")
  if (!is.null(m0) && !is.null(e0) &&
      !any(is.na(c(m0$OR, m0$下限, m0$上限,
                   e0$OR, e0$下限, e0$上限)))) {
    ors <- c(m0$OR, e0$OR)
    los <- c(m0$下限, e0$下限)
    his <- c(m0$上限, e0$上限)
    yi  <- log(ors)
    ## 95% CI 的两端间距对应 2 × 1.96 个标准误
    sei <- (log(his) - log(los)) / (2 * 1.96)
    wi  <- 1 / sei^2
    ybar <- sum(wi * yi) / sum(wi)
    Q    <- sum(wi * (yi - ybar)^2)
    k    <- length(yi)
    C    <- sum(wi) - sum(wi^2) / sum(wi)
    tau2 <- max(0, (Q - (k - 1)) / C)
    ws   <- 1 / (sei^2 + tau2)
    mu   <- sum(ws * yi) / sum(ws)
    sem  <- sqrt(1 / sum(ws))
    o$双库汇总 <- list(
      模型 = "基础模型",
      方法 = "随机效应（DerSimonian–Laird），log OR 尺度",
      说明 = paste0("二次计算：两库各自独立筛选、独立建模，",
                    "此处仅汇总两库已得的效应量，原始记录未做任何合并"),
      库数 = k,
      OR = r4(exp(mu), 3),
      下限 = r4(exp(mu - 1.96 * sem), 3),
      上限 = r4(exp(mu + 1.96 * sem), 3),
      tau2 = r4(tau2, 4),
      Q = r4(Q, 3),
      Q的P = r4(pchisq(Q, df = k - 1, lower.tail = FALSE), 3),
      I2 = r4(if (Q > 0) 100 * max(0, (Q - (k - 1)) / Q) else 0, 1))
  }
}

if (!requireNamespace("jsonlite", quietly = TRUE))
  stop("缺少 jsonlite 包，无法写 results.json")
jsonlite::write_json(o, file.path(OUT_DIR, "results.json"),
                     auto_unbox = TRUE, pretty = TRUE, null = "null",
                     digits = NA)
cat(sprintf("已写 %s（%d 个顶层字段）\n",
            file.path(OUT_DIR, "results.json"), length(o)))
