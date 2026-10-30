## eICU 侧敏感性分析（S1–S11）+ 竞争风险
##
## 与 MIMIC 侧 _run4.R 同构。**唯一实质差别在 S4/S4b**：本库 lab 表里 sodium 只有
## 血清/血浆一种，没有血气钠这一路，故 extract_eicu_all() 返回零行的 sodium_blood，
## run_sensitivity() 的 `nrow() > 0` 判据会让 S4/S4b 自然跳过——是"本库不适用"，
## 不是"跑失败"。成稿时两库并排的敏感性表要据此说明，不能留空白行。
source("extract.R")
source("sensitivity.R")
source("competing.R")

args <- commandArgs(trailingOnly = TRUE)
CC <- if (length(args) >= 1 && nzchar(args[1])) toupper(args[1]) else "B"

ex   <- readRDS(file.path(OUT_DIR, "eicu_B_extract.rds"))
full <- readRDS(file.path(OUT_DIR, sprintf("eicu_%s_full.rds", CC)))
a    <- full$analysis

## ---- 主分析对照行（基础模型，与 run_sensitivity 同一口径）----
am <- a[model_cases(a), , drop = FALSE]
xs <- c("delta_na", .pick(BASIC_X, am))
m0 <- glm(.fml("aki", xs), data = am, family = binomial())
c0 <- summary(m0)$coefficients["delta_na", ]
main_row <- data.frame(分析 = sprintf("主分析（表型 %s）", CC), n = nobs(m0),
                       事件数 = sum(am$aki),
                       OR = exp(c0[1]),
                       下限 = exp(c0[1] - 1.96 * c0[2]),
                       上限 = exp(c0[1] + 1.96 * c0[2]),
                       P = c0[4], 备注 = "", row.names = NULL)

cat("=== 敏感性分析（S1–S10 + 附加）===\n")
## cohort = CC 必须传：每一项变体都会重跑队列组装，默认口径是 B，
## 漏传会让整张敏感性表退回另一个表型的人群（见 sensitivity.R 的说明）
st <- run_sensitivity(ex, main_row, d_main = a, cohort = CC)
print_sensitivity(st)
saveRDS(st, file.path(OUT_DIR, sprintf("eicu_%s_sensitivity.rds", CC)))

## ---- S11 竞争风险 ----
cat("\n=== S11 竞争风险 ===\n")
d <- build_competing(am)
cat(sprintf("  进入风险集 %d 例；AKI %d，转出/死亡 %d，删失 %d\n",
            nrow(d), sum(d$status == 1L), sum(d$status == 2L), sum(d$status == 0L)))

fg <- fit_finegray(d)
r <- fg[fg$项 == "delta_na", ]
cat(sprintf("  Fine-Gray 次分布风险：sHR %.3f（%.3f–%.3f）  P = %s\n",
            r$sHR, r$下限, r$上限, format.pval(r$P, digits = 2, eps = 1e-4)))

aj <- fit_aalen_johansen(d)
cat("  72h 累积发生率（Aalen-Johansen，不做模型假设）：\n")
for (i in seq_len(nrow(aj$汇总)))
  cat(sprintf("    %-12s %.4f\n", aj$汇总$状态[i], aj$汇总$累积发生率[i]))

saveRDS(list(finegray = fg, aj_summary = aj$汇总, competing_data = d),
        file.path(OUT_DIR, sprintf("eicu_%s_competing.rds", CC)))
cat(sprintf("\n已存 eicu_%s_sensitivity.rds / eicu_%s_competing.rds\n", CC, CC))
