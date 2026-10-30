## 主分析：三档模型 + RCS 剂量反应 + 预设亚组（读 _run2 落盘的 rds）
source("extract.R")   # 只为拿 OUT_DIR 等路径常量，不会连库
source("models.R")
source("tables.R")

f <- file.path(OUT_DIR, "mimic_B_full.rds")
stopifnot(file.exists(f))
a0 <- readRDS(f)$analysis
cat(sprintf("分析人群 %d，事件 %d（%.1f%%）\n", nrow(a0), sum(a0$aki),
            100 * mean(a0$aki)))

## ---- 三档模型统一分母 ----
## 见 models.R 的 model_cases()：glm 的按行删除是静默的，不先切一次，
## 未调整/基础/扩展三个模型的分母就会各不相同。
keep <- model_cases(a0)
a <- a0[keep, , drop = FALSE]
cat(sprintf("模型人群 %d，事件 %d（%.1f%%）；因协变量缺失剔除 %d 例\n",
            nrow(a), sum(a$aki), 100 * mean(a$aki), sum(!keep)))
if (any(!keep)) {
  ## 剔除由哪些变量的缺失造成——论文里要写清楚，否则"少了 44 人"无法交代
  miss <- vapply(.pick(MODEL_VARS, a0),
                 function(v) sum(is.na(a0[[v]]) & keep == FALSE), integer(1))
  miss <- sort(miss[miss > 0], decreasing = TRUE)
  if (length(miss))
    cat(sprintf("  剔除原因（按变量计的缺失人次，一人可占多项）：%s\n",
                paste(sprintf("%s %d", names(miss), miss), collapse = "；")))
}
cat("\n")

cat("=== 基线表（按 ΔNa 三档）===\n")
print_baseline(baseline_table(a))
cat("\n")

## ---- 主模型三档 ----
cat("=== 主模型（暴露：ΔNa 每升高 1 mmol/L）===\n")
ti <- fit_tiers(a)
for (nm in unique(ti$模型)) {
  r <- ti[ti$模型 == nm & ti$项 == "delta_na", ]
  cat(sprintf("  %-6s OR %.3f（%.3f–%.3f）  P = %s   n = %d\n",
              nm, r$OR, r$下限, r$上限,
              format.pval(r$P, digits = 2, eps = 1e-4), r$n))
}
cat("\n--- 基础模型全部系数 ---\n")
b <- ti[ti$模型 == "基础", ]
for (i in seq_len(nrow(b)))
  cat(sprintf("  %-26s OR %.3f（%.3f–%.3f）  P = %s\n",
              b$项[i], b$OR[i], b$下限[i], b$上限[i],
              format.pval(b$P[i], digits = 2, eps = 1e-4)))

## ---- RCS ----
cat("\n=== RCS 剂量反应（4 节点，参考 ΔNa = 0）===\n")
rc <- fit_rcs(a)
cat(sprintf("  节点（分位数 %.0f/%.0f/%.0f/%.0f%%）：%s\n",
            100 * RCS_KNOTS_Q[1], 100 * RCS_KNOTS_Q[2],
            100 * RCS_KNOTS_Q[3], 100 * RCS_KNOTS_Q[4],
            paste(sprintf("%.2f", rc$knots), collapse = ", ")))
cat(sprintf("  整体 P = %s   非线性 P = %s\n",
            format.pval(rc$p_overall, digits = 3, eps = 1e-5),
            format.pval(rc$p_nonlinear, digits = 3, eps = 1e-5)))
cv <- rcs_curve(rc)
for (x in c(-6, -3, 0, 3, 6)) {
  i <- which.min(abs(cv$delta_na - x))
  cat(sprintf("  ΔNa = %+5.1f  OR %.3f（%.3f–%.3f）\n",
              cv$delta_na[i], cv$OR[i], cv$下限[i], cv$上限[i]))
}

## ---- 预设亚组 ----
cat("\n=== 预设亚组（结论标注探索性）===\n")
sg <- fit_subgroups(a)
for (i in seq_len(nrow(sg)))
  cat(sprintf("  %-22s %-7s n=%5d 事件=%4d  %s\n", sg$亚组[i], sg$层[i],
              sg$n[i], sg$事件数[i],
              if (is.na(sg$OR[i])) sprintf("交互 P = %s", format.pval(sg$P[i], digits = 2, eps = 1e-4))
              else sprintf("OR %.3f（%.3f–%.3f） P = %s", sg$OR[i], sg$下限[i], sg$上限[i],
                           format.pval(sg$P[i], digits = 2, eps = 1e-4))))

saveRDS(list(tiers = ti, rcs = rc, curve = cv, subgroups = sg),
        file.path(OUT_DIR, "mimic_B_models.rds"))
cat(sprintf("\n已存 %s\n", file.path(OUT_DIR, "mimic_B_models.rds")))
