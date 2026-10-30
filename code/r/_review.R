## 随机抽 20 例手工核对 ΔNa 与 KDIGO 判定（计划中的验证步骤）
source("cohort.R")
res <- readRDS(file.path(OUT_DIR, "mimic_cohort_B.rds"))
a <- res$analysis

set.seed(20260930)
samp <- a[sample(nrow(a), 20), c("id", "delta_na", "na_first", "na_last",
                                 "aki", "t_event", "stage_first", "stage_max",
                                 "low48_at_event", "low7_at_event", "creat_first")]
samp <- samp[order(samp$id), ]

cr <- res$creatinine   # 从已落盘的提取结果读取，不再重扫原始表

for (i in seq_len(nrow(samp))) {
  id <- samp$id[i]
  x <- cr[cr$id == id, ]; x <- x[order(x$t_hours), ]
  s <- .kdigo_series(x$t_hours, x$value)
  cat(sprintf("\n--- %s | ΔNa=%+.1f (%.0f→%.0f) | aki=%s t_event=%s stage=%s ---\n",
              id, samp$delta_na[i], samp$na_first[i], samp$na_last[i],
              samp$aki[i], ifelse(is.na(samp$t_event[i]), "-", samp$t_event[i]), samp$stage_first[i]))
  show <- s[order(abs(s$t_hours))[1:min(8, nrow(s))], ]
  show <- show[order(show$t_hours), ]
  for (j in seq_len(nrow(show)))
    cat(sprintf("   t=%8.1fh  Cr=%5.2f  low48=%5s  low7=%5s  stage=%s%s\n",
                show$t_hours[j], show$creat[j],
                ifelse(is.na(show$low48[j]), "  -", sprintf("%.2f", show$low48[j])),
                ifelse(is.na(show$low7[j]),  "  -", sprintf("%.2f", show$low7[j])),
                show$stage[j],
                ifelse(!is.na(samp$t_event[i]) && show$t_hours[j] == samp$t_event[i], "   <== 事件", "")))
}
