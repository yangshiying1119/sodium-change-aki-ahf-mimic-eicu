## 跑 MIMIC-IV 主队列并落盘（提取一次，后续复核/建模直接读 rds）
source("cohort.R")
res <- run_mimic_cohort("B")
print_attrition(res)
cat("\n=== 化验长表规模 ===\n")
cat(sprintf("  血钠 %d 行 / %d stay\n", nrow(res$sodium), length(unique(res$sodium$id))))
cat(sprintf("  肌酐 %d 行 / %d stay\n", nrow(res$creatinine), length(unique(res$creatinine$id))))

saveRDS(list(analysis = res$analysis, data = res$data, attrition = res$attrition,
             sodium = res$sodium, creatinine = res$creatinine, na_df = res$na_df),
        file.path(OUT_DIR, "mimic_cohort_B.rds"))
cat(sprintf("\n已存 %s\n", file.path(OUT_DIR, "mimic_cohort_B.rds")))
