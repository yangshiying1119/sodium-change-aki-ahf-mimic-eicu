## =============================================================================
## _diag_s4.R —— S4（血气钠替换生化钠）出现零结果的归因诊断
##
## S4 结果与主分析明显不同（OR 1.016 vs 0.933），须区分两种可能：
##   (甲) **选择效应**：血气只在较重的患者身上查（血气钠 ≥2 次者仅占分析人群
##        约两成），这个亚群本身的暴露-结局关系就不同于全人群；
##   (乙) **测量误差**：血气钠与生化钠存在系统偏倚与随机误差，ΔNa 被稀释，
##        效应量向 OR=1 衰减（回归稀释）。
## 两者含义完全不同：甲说明 S4 的零结果不代表"标本类型影响"，乙才是。
##
## 判别方法：
##   1. 配对比较——同一患者、相近时刻两种标本的钠值之差（系统偏倚与离散度）；
##   2. **把主分析原样限制到 S4 的子人群**（仍用生化钠）。若该子集里 OR 已经
##      接近 1，则 S4 的零结果主要来自选择；若仍接近 0.933，则来自测量误差。
## =============================================================================

source("extract.R"); source("cohort.R"); source("models.R")

ex   <- readRDS(file.path(OUT_DIR, "mimic_extract_B.rds"))
full <- readRDS(file.path(OUT_DIR, "mimic_B_full.rds"))
a    <- full$analysis
adj  <- .pick(BASIC_X, a)

.or <- function(d, lbl) {
  m  <- glm(.fml("aki", c("delta_na", adj)), data = d, family = binomial())
  co <- summary(m)$coefficients["delta_na", ]
  cat(sprintf("  %-34s n=%5d 事件=%4d  OR %.3f（%.3f–%.3f）  P=%s\n",
              lbl, nobs(m), sum(d$aki), exp(co[1]),
              exp(co[1] - 1.96 * co[2]), exp(co[1] + 1.96 * co[2]),
              format.pval(co[4], digits = 2, eps = 1e-4)))
}

## ---- 1. 配对比较：同一患者、同一时刻（±0.25h）两种标本 ----
cat("=== 1. 生化钠 vs 血气钠 配对比较（同患者、时刻差 ≤0.25h）===\n")
ch <- ex$sodium;       names(ch)[names(ch) == "value"] <- "na_chem"
bg <- ex$sodium_blood; names(bg)[names(bg) == "value"] <- "na_bg"
ch <- ch[!is.na(ch$na_chem), ]; bg <- bg[!is.na(bg$na_bg), ]
m <- merge(ch, bg, by = "id", suffixes = c("", ".bg"))
m <- m[abs(m$t_hours - m$t_hours.bg) <= SAME_TIME_TOL_H, ]
if (nrow(m)) {
  m <- m[order(m$id, abs(m$t_hours - m$t_hours.bg)), ]
  m <- m[!duplicated(paste(m$id, m$t_hours)), ]
  d <- m$na_chem - m$na_bg
  cat(sprintf("  配对点数：%d（覆盖 %d 名患者）\n", nrow(m), length(unique(m$id))))
  cat(sprintf("  生化钠均值 %.2f，血气钠均值 %.2f\n", mean(m$na_chem), mean(m$na_bg)))
  cat(sprintf("  差值（生化−血气）：均值 %+.2f，SD %.2f，中位 %+.2f，四分位距 %.2f–%.2f\n",
              mean(d), sd(d), median(d), quantile(d, .25), quantile(d, .75)))
  cat(sprintf("  相关系数 r = %.3f\n", cor(m$na_chem, m$na_bg)))
  cat(sprintf("  差值绝对值 >2 mmol/L 的比例：%.1f%%；>4：%.1f%%\n",
              100 * mean(abs(d) > 2), 100 * mean(abs(d) > 4)))
  ## 两种标本各自测得的 ΔNa 在配对子集里的一致性
  cat("\n  同患者两种标本算出的 ΔNa（配对子集内，各自 0–24h 末次−首次）：\n")
  f <- function(dt, v) {
    dt <- dt[!is.na(dt[[v]]) & dt$t_hours >= 0 & dt$t_hours <= LANDMARK_H, ]
    r <- do.call(rbind, lapply(split(dt, dt$id), function(x) {
      x <- x[order(x$t_hours), ]
      if (nrow(x) < 2) return(NULL)
      data.frame(id = x$id[1], d = x[[v]][nrow(x)] - x[[v]][1])
    }))
    r
  }
  d1 <- f(ex$sodium, "value"); d2 <- f(ex$sodium_blood, "value")
  dd <- merge(d1, d2, by = "id", suffixes = c("_chem", "_bg"))
  if (nrow(dd) > 5) {
    cat(sprintf("    可同时算出者 %d 例；ΔNa 相关 r = %.3f\n",
                nrow(dd), cor(dd$d_chem, dd$d_bg)))
    cat(sprintf("    ΔNa 差值（生化−血气）：均值 %+.2f，SD %.2f\n",
                mean(dd$d_chem - dd$d_bg), sd(dd$d_chem - dd$d_bg)))
    cat(sprintf("    两法 ΔNa 符号不一致者：%d 例（%.1f%%）\n",
                sum(sign(dd$d_chem) != sign(dd$d_bg)),
                100 * mean(sign(dd$d_chem) != sign(dd$d_bg))))
  }
} else cat("  无配对点\n")

## ---- 2. 把主分析限制到 S4 子人群（仍用生化钠） ----
cat("\n=== 2. 选择效应检验：主分析（生化钠）限制到 S4 子人群 ===\n")
n_bg <- table(ex$sodium_blood$id[ex$sodium_blood$t_hours >= 0 &
                                ex$sodium_blood$t_hours <= LANDMARK_H])
have_bg2 <- names(n_bg)[n_bg >= 2]
a$in_s4 <- a$id %in% have_bg2
.or(a, "全分析人群（主分析）")
.or(a[a$in_s4, , drop = FALSE], "其中：有 ≥2 次血气钠者")
.or(a[!a$in_s4, , drop = FALSE], "其中：无 ≥2 次血气钠者")

cat("\n=== 3. 两组的基线差异（提示选择方向）===\n")
for (v in c("age", "na_first", "creat_baseline", "sbp_min_24h", "gcs_min_24h",
            "mech_vent_any", "vasoactive_any", "aki")) {
  if (!v %in% names(a)) next
  g1 <- a[[v]][a$in_s4]; g0 <- a[[v]][!a$in_s4]
  cat(sprintf("  %-18s 有血气 %7.2f ／ 无血气 %7.2f\n", v, mean(g1, na.rm = TRUE),
              mean(g0, na.rm = TRUE)))
}
