## 随机抽 20 例手工核对 ΔNa 与 KDIGO 判定（计划中的验证步骤）
##
## 目的不是再算一遍，而是**把判定所依据的原始序列摆出来**，让复核者能逐行手工验算：
##   - 血钠：0–24h 内逐次实测值（同时间已按方案合并）、首末取值、ΔNa
##   - 肌酐：界标前（含院前负值）到 72h 的完整序列，每点附该点的 48h/7d 比较低值
##     与判出的分期；事件点单独标注
## 输出为纯文本，便于逐条比对与留档（属内部产物，不进客户交付物）。
source("cohort.R")
source("tables.R")

full <- readRDS(file.path(OUT_DIR, "mimic_B_full.rds"))
res  <- full$data           # 全表型人群，含未纳入者，便于看排除过程
a    <- full$analysis

set.seed(20261001)
samp <- a[sample(nrow(a), 20), , drop = FALSE]
samp <- samp[order(samp$id), , drop = FALSE]

sod_all <- merge_same_time(full$sodium)
cr_all  <- merge_same_time(full$creatinine)

out <- character(0)
say <- function(...) out <<- c(out, sprintf(...))

say("随机抽 20 例判定复核（seed = 20261001）")
say("分组列：ΔNa = 0–24h 末次 − 首次血钠；AKI = 界标后 (24, 72] 首次满足 KDIGO 肌酐标准")
say("")

for (i in seq_len(nrow(samp))) {
  id <- samp$id[i]
  say("================================================================")
  say("第 %2d 例   stay_id = %s", i, id)
  say("  ΔNa = %+.1f mmol/L（首次 %.1f → 末次 %.1f，检测 %d 次，首末间隔 %.1f h）",
      samp$delta_na[i], samp$na_first[i], samp$na_last[i],
      samp$n_na[i], samp$t_last[i] - samp$t_first[i])
  say("  界标前最低肌酐 %.2f mg/dL；界标前肌酐 %d 次、界标后 %d 次",
      samp$creat_baseline[i], samp$n_creat_before[i], samp$n_creat_post[i])
  say("  判定：AKI = %s；事件时刻 = %s；首次分期 = %s；最高分期 = %s",
      ifelse(samp$aki[i], "是", "否"),
      ifelse(is.na(samp$t_event[i]), "—", sprintf("%.1f h", samp$t_event[i])),
      samp$stage_first[i], samp$stage_max[i])
  if (is.na(samp$t_event[i]))
    say("        （事件点比较值：48h 最低 %s，7d 最低 %s）",
        ifelse(is.na(samp$low48_at_event[i]), "—", sprintf("%.2f", samp$low48_at_event[i])),
        ifelse(is.na(samp$low7_at_event[i]),  "—", sprintf("%.2f", samp$low7_at_event[i])))
  else
    say("        （事件点肌酐 %.2f；48h 最低 %.2f；7d 最低 %.2f）",
        samp$creat_first[i], samp$low48_at_event[i], samp$low7_at_event[i])

  s <- sod_all[sod_all$id == id, , drop = FALSE]
  s <- s[s$t_hours >= 0 & s$t_hours <= LANDMARK_H, , drop = FALSE]
  say("  --- 0–24h 血钠（%d 次）---", nrow(s))
  for (j in seq_len(nrow(s)))
    say("      t=%6.2f h   Na=%5.1f%s", s$t_hours[j], s$value[j],
        ifelse(j == 1L, "   <== 首次（暴露起点）",
        ifelse(j == nrow(s), "   <== 末次（暴露终点）", "")))

  x <- cr_all[cr_all$id == id, , drop = FALSE]
  x <- x[x$t_hours <= HORIZON_H, , drop = FALSE]
  x <- x[order(x$t_hours), , drop = FALSE]
  k <- .kdigo_series(x$t_hours, x$value)
  say("  --- 界标前至 72h 肌酐（%d 次）---", nrow(k))
  for (j in seq_len(nrow(k)))
    say("      t=%7.2f h   Cr=%5.2f   48h最低=%5s   7d最低=%5s   分期=%s%s",
        k$t_hours[j], k$creat[j],
        ifelse(is.na(k$low48[j]), "    -", sprintf("%.2f", k$low48[j])),
        ifelse(is.na(k$low7[j]),  "    -", sprintf("%.2f", k$low7[j])),
        k$stage[j],
        ifelse(!is.na(samp$t_event[i]) && abs(k$t_hours[j] - samp$t_event[i]) < 1e-9,
               "   <== 事件", ""))
  say("")
}

say("================================================================")
say("复核要点：")
say("  1) ΔNa 是否等于末次减首次，且两次都在 0–24h 内")
say("  2) 每一点的 48h/7d 比较低值是否**只取自该点之前**的实测值（不得包含之后的值）")
say("  3) 事件点是否确实落在 (24, 72]，且其分期确由 48h +0.3 或 1.5/2/3 倍基线之一触发")
say("  4) 界标前已满足 KDIGO 者是否已被排除（不在本表中）")

f <- file.path(OUT_DIR, "logs", "review20.txt")
writeLines(out, f, useBytes = TRUE)
cat(paste(out, collapse = "\n"), "\n")
cat(sprintf("\n已写入 %s\n", f))
