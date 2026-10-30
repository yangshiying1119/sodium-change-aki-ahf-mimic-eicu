source("cohort.R")
ex <- extract_mimic("B"); ids <- ex$cohort$id
cr  <- ex$creatinine[ex$creatinine$id %in% ids, ]
sod <- ex$sodium[ex$sodium$id %in% ids, ]
nd  <- delta_na(sod)                     # 可算 ΔNa 的人群
nd  <- nd[!is.na(nd$delta_na), c("id","delta_na")]

seq_by_id <- split(cr, cr$id)
cat("=== 界标前窗口敏感性（分析人群 = 可算ΔNa 且 界标前未达标 且 界标后有复测）===\n")
cat(sprintf("  可算 ΔNa: %d 例；其中有肌酐序列 %d 例\n",
            nrow(nd), sum(nd$id %in% names(seq_by_id))))

for (lo in c(0, -168, -Inf)) {
  rows <- lapply(names(seq_by_id), function(id) {
    x <- seq_by_id[[id]]; s <- .kdigo_series(x$t_hours, x$value)
    pre  <- s$t_hours >= lo & s$t_hours <= 24
    post <- s$t_hours >  24 & s$t_hours <= 72
    data.frame(id = id,
               aki_pre = any(pre & !is.na(s$stage) & s$stage >= 1),
               aki     = any(post & !is.na(s$stage) & s$stage >= 1),
               n_post  = sum(post), stringsAsFactors = FALSE)
  })
  d <- merge(nd, do.call(rbind, rows), by = "id")
  ev <- d[!d$aki_pre & d$n_post > 0, ]
  cat(sprintf("  界标前窗 [%-5s, 24]: aki_pre 排除 %5d → 分析 %5d, 事件 %4d (%.1f%%)\n",
              if (is.infinite(lo)) "-Inf" else as.character(lo),
              sum(d$aki_pre), nrow(ev), sum(ev$aki), 100*mean(ev$aki)))
}
