## eICU 侧的 S8（尿量标准的复合结局）+ 选择与测量偏倚分析。
## 与 MIMIC 侧 _run7.R 同构，两库共用 sensitivity.R / bias.R 里的同一批函数；
## 唯一的差别是**每一处会重跑队列组装的调用都必须显式带上本库的表型口径**——
## run_mimic_cohort 的默认口径是 "B"，漏传会让整张表退回另一个表型的人群。
##
## 用法：Rscript _runE7.R [口径字母]，缺省 C
source("extract.R")
source("sensitivity.R")   # 内含 cohort.R / models.R
source("bias.R")

args <- commandArgs(trailingOnly = TRUE)
CC   <- if (length(args) >= 1 && nzchar(args[1])) toupper(args[1]) else "C"

ex   <- readRDS(file.path(OUT_DIR, "eicu_B_extract.rds"))   # 缓存与口径无关
full <- readRDS(file.path(OUT_DIR, sprintf("eicu_%s_full.rds", CC)))
d    <- attach_covariates(full$data, ex$covariates)

## ---- S8 复合口径 ----
cat("=== S8 结局改用复合口径（肌酐标准 ∪ 尿量标准）===\n")
s8 <- run_s8_uo(ex, cohort = CC)
cat(sprintf("  尿量长表 %d 行；有体重者 %d / %d\n",
            nrow(ex$urine), sum(!is.na(ex$weight)), length(ex$weight)))
for (i in seq_len(nrow(s8$队列$attrition)))
  cat(sprintf("  %-46s %10s %10d\n", s8$队列$attrition$节点[i],
              ifelse(is.na(s8$队列$attrition$排除[i]), "—",
                     s8$队列$attrition$排除[i]), s8$队列$attrition$留存[i]))
print(s8$交叉)
r <- s8$行
cat(sprintf("  基础模型：OR %.3f（%.3f–%.3f）  P = %s   n = %d，事件 = %d\n",
            r$OR, r$下限, r$上限, format.pval(r$P, digits = 2, eps = 1e-4),
            r$n, r$事件数))
cat(sprintf("  %s\n", r$备注))
## 表型口径写进分析名，两库并排时才看得出是哪一版
s8$行$分析 <- sprintf("S8 复合口径（肌酐 ∪ 尿量，表型 %s）", CC)
saveRDS(s8, file.path(OUT_DIR, sprintf("eicu_%s_s8_uo.rds", CC)))

## ---- 偏倚分析 ①：可判定 vs 不可判定 ----
cat("\n=== 偏倚分析① 可判定结局者与各类未纳入者的比较 ===\n")
cls <- .eligibility_class(d)
print(table(cls))
cmp <- compare_evaluable(d, cls)
for (i in seq_len(nrow(cmp)))
  cat(sprintf("  %-24s %s\n", cmp$变量[i],
              paste(sprintf("%-16s", unlist(cmp[i, -1])), collapse = " ")))
saveRDS(list(分组 = cls, 对比 = cmp),
        file.path(OUT_DIR, sprintf("eicu_%s_bias_compare.rds", CC)))

cat("\n--- 血钠监测强度与复测/AKI 检出（方案 5.7 检出偏倚）---\n")
rm_ <- retest_by_monitoring(d)
for (i in seq_len(nrow(rm_)))
  cat(sprintf("  %-8s n=%5d  复测 %-8s  入组 %-8s  AKI %s\n",
              rm_$血钠检测次数[i], rm_$n[i], rm_$界标后复测比例[i],
              rm_$分析人群占比[i], rm_$AKI检出率[i]))

## ---- 偏倚分析 ②：复测概率权重 ----
cat("\n=== 偏倚分析② 复测概率权重（IPW，探索性）===\n")
w <- fit_retest_ipw(d)
print_retest_ipw(w)
ip <- ips_retest_model(d, w)
cat(sprintf("\n  %s：OR %.3f（%.3f–%.3f）  P = %s   n = %d，事件 = %d\n",
            ip$分析, ip$OR, ip$下限, ip$上限,
            format.pval(ip$P, digits = 2, eps = 1e-4), ip$n, ip$事件数))
## 对照值现算，不写死——主分析的 OR 会随口径与完整病例集变动
a_ref <- d[d$evaluable %in% TRUE, , drop = FALSE]
a_ref <- a_ref[model_cases(a_ref), , drop = FALSE]
m_ref <- glm(.fml("aki", c("delta_na", .pick(BASIC_X, a_ref))),
             data = a_ref, family = binomial())
co_ref <- summary(m_ref)$coefficients["delta_na", ]
cat(sprintf("  对照：主分析（未加权，同一批 %d 例）OR %.3f（%.3f–%.3f）\n",
            nrow(a_ref), exp(co_ref[1]),
            exp(co_ref[1] - 1.96 * co_ref[2]), exp(co_ref[1] + 1.96 * co_ref[2])))
saveRDS(list(ipw = w, 模型 = ip, 对比 = cmp, 监测 = rm_),
        file.path(OUT_DIR, sprintf("eicu_%s_bias.rds", CC)))

cat(sprintf("\n已存 %s / %s\n",
            file.path(OUT_DIR, sprintf("eicu_%s_s8_uo.rds", CC)),
            file.path(OUT_DIR, sprintf("eicu_%s_bias.rds", CC))))
