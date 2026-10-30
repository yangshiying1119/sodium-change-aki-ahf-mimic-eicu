## A 口径在**方案完整版**流水线下的队列与事件数
## 用途：修订一（主分析 A→B）的依据此前用的是作废版数字，须重算后重述。
## 提取层与 B 共用同一份缓存，故此处只重跑组装，不重扫库。
source("cohort.R")
source("models.R")

ex <- readRDS(file.path(OUT_DIR, "mimic_extract_B.rds"))

for (coh in c("A", "B", "C")) {
  r <- run_mimic_cohort(coh, ex = ex)
  a <- r$analysis
  cat(sprintf("\n===== 表型 %s =====\n", coh))
  for (i in seq_len(nrow(r$attrition)))
    cat(sprintf("  %-46s %10s %10d\n", r$attrition$节点[i],
                ifelse(is.na(r$attrition$排除[i]), "—", r$attrition$排除[i]),
                r$attrition$留存[i]))
  cat(sprintf("  分析人群 %d，事件 %d（%.1f%%），EPV(10) = %.1f\n",
              nrow(a), sum(a$aki), 100 * mean(a$aki), sum(a$aki) / 10))
  ## 主模型（未调整 + 基础），口径与主分析一致
  if (nrow(a) > 100 && sum(a$aki) > 30) {
    m0 <- glm(aki ~ delta_na, data = a, family = binomial())
    c0 <- summary(m0)$coefficients["delta_na", ]
    cat(sprintf("  未调整 OR %.3f（%.3f–%.3f）P = %s\n",
                exp(c0[1]), exp(c0[1] - 1.96 * c0[2]), exp(c0[1] + 1.96 * c0[2]),
                format.pval(c0[4], digits = 2, eps = 1e-4)))
  }
  saveRDS(r, file.path(OUT_DIR, sprintf("mimic_cohort_%s.rds", coh)))
}
