## 敏感性分析（S1–S11）+ 竞争风险 + 偏倚分析
source("extract.R")
source("sensitivity.R")
source("competing.R")

ex <- readRDS(file.path(OUT_DIR, "mimic_extract_B.rds"))
full <- readRDS(file.path(OUT_DIR, "mimic_B_full.rds"))
a <- full$analysis

## ---- 主分析对照行（基础模型，与 run_sensitivity 同一口径）----
## 与 _run3.R 一样先取完整病例，否则这一行会和表 2 的"基础"行对不上
## （一个 2683、一个 2656），审稿人一眼就看出来。
am <- a[model_cases(a), , drop = FALSE]
xs <- c("delta_na", .pick(BASIC_X, am))
m0 <- glm(.fml("aki", xs), data = am, family = binomial())
c0 <- summary(m0)$coefficients["delta_na", ]
main_row <- data.frame(分析 = "主分析（表型 B）", n = nobs(m0),
                       事件数 = sum(am$aki),
                       OR = exp(c0[1]),
                       下限 = exp(c0[1] - 1.96 * c0[2]),
                       上限 = exp(c0[1] + 1.96 * c0[2]),
                       P = c0[4], 备注 = "", row.names = NULL)

cat("=== 敏感性分析（S1–S10 + 附加）===\n")
## d_main 传的是**未按 model_cases 剔除**的 a：多重插补的意义就在于把那批人补回来
## cohort 显式传 B：与主队列主分析的表型一致。默认值虽也是 B，
## 但显式写出后，将来若改动默认值不会静默改变本表的人群。
st <- run_sensitivity(ex, main_row, d_main = a, cohort = "B")
print_sensitivity(st)
saveRDS(st, file.path(OUT_DIR, "mimic_sensitivity.rds"))

## ---- S11 竞争风险 ----
cat("\n=== S11 竞争风险 ===\n")
d <- build_competing(am)
cat(sprintf("  进入风险集 %d 例；AKI %d，转出/死亡 %d，删失 %d\n",
            nrow(d), sum(d$status == 1L), sum(d$status == 2L), sum(d$status == 0L)))

fg <- fit_finegray(d)
r <- fg[fg$项 == "delta_na", ]
cat(sprintf("  Fine-Gray 次分布风险：sHR %.3f（%.3f–%.3f）  P = %s\n",
            r$sHR, r$下限, r$上限, format.pval(r$P, digits = 2, eps = 1e-4)))
cat("  （对照：普通 logistic OR 见上表主分析行）\n")

aj <- fit_aalen_johansen(d)
cat("  72h 累积发生率（Aalen-Johansen，不做模型假设）：\n")
for (i in seq_len(nrow(aj$汇总)))
  cat(sprintf("    %-12s %.4f\n", aj$汇总$状态[i], aj$汇总$累积发生率[i]))

saveRDS(list(finegray = fg, aj_summary = aj$汇总, competing_data = d),
        file.path(OUT_DIR, "mimic_competing.rds"))
cat(sprintf("\n已存 %s\n", file.path(OUT_DIR, "mimic_sensitivity.rds")))
