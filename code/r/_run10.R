## 判别力与校准（方案 5.3）。两库各跑一次，各落一份 RDS 供 _run9.R 汇总。
##
## 用法：Rscript _run10.R [mimic|eicu] [口径字母]
##   MIMIC 侧口径固定 B，eICU 侧缺省 C（与 _run9.R 的 EICU_CC 一致）。
source("extract.R")
source("common.R")
source("cohort.R")
source("models.R")
source("discrim.R")

args <- commandArgs(trailingOnly = TRUE)
db   <- if (length(args) >= 1 && nzchar(args[1])) tolower(args[1]) else "mimic"
cc   <- if (length(args) >= 2 && nzchar(args[2])) toupper(args[2]) else
        if (db == "mimic") "B" else "C"

full <- readRDS(file.path(OUT_DIR, sprintf("%s_%s_full.rds", db, cc)))

## 直接用分析人群取完整病例：full$analysis 已并过协变量，
## 与三档模型是同一批人，AUC 与主分析的 OR 才在同一个分母上。
d <- full$analysis
d <- d[model_cases(d), , drop = FALSE]
cat(sprintf("%s 口径 %s：完整病例 %d 例，事件 %d\n", db, cc, nrow(d), sum(d$aki)))

r <- discrim_all(d)

cat(sprintf("每 3 mmol/L：OR %.3f（%.3f–%.3f），P = %.4g\n",
            r$每3mmol$OR, r$每3mmol$下限, r$每3mmol$上限, r$每3mmol$P))
cat(sprintf("似然比检验：χ² = %.2f，P = %.4g\n", r$LRT$卡方, r$LRT$P))
cat(sprintf("AUC %.3f（%.3f–%.3f）；去掉 ΔNa 后 %.3f；ΔAUC %.4f，P = %.4g\n",
            r$判别$AUC, r$判别$AUC下限, r$判别$AUC上限, r$判别$AUC无暴露,
            r$判别$差值, r$判别$差值P))
cat(sprintf("Hosmer–Lemeshow：χ² = %.2f，df = %d，P = %.4g（最小期望 %.1f）\n",
            r$HL$卡方, r$HL$自由度, r$HL$P, r$HL$格子最小期望))
cat(sprintf("表观校准斜率（恒为 1，仅作对照）：%.3f，表观截距 %.3f\n",
            r$校准$表观斜率, r$校准$表观截距))
if (!is.null(r$校正))
  cat(sprintf("bootstrap 校正（B=%d）：c 统计量 %.3f（%.3f–%.3f），乐观量 %.4f，校准斜率 %.3f\n",
              r$校正$B, r$校正$c校正, r$校正$c校正下限, r$校正$c校正上限,
              r$校正$乐观量, r$校正$斜率校正))
cat("校准曲线（十分位）：\n"); print(round(r$校准曲线[, -1], 3))
cat("标准化预测概率（%）：\n"); print(round(100 * r$标准化概率[, -1, drop = FALSE], 1))
cat("层内人数：\n"); print(r$层内人数)

## 模型对象不入 RDS（体积大且含 data 副本），只留出图要用的预测值
out <- r[c("n", "事件数", "每3mmol", "LRT", "HL", "校准", "校准曲线",
           "判别", "校正", "标准化概率", "层内人数")]
out$预测 = data.frame(结局 = as.integer(d$aki),
                      预测概率 = stats::fitted(r$模型),
                      线性预测子 = stats::predict(r$模型, type = "link"))
saveRDS(out, file.path(OUT_DIR, sprintf("%s_%s_discrim.rds", db, cc)))
cat(sprintf("已写 %s\n", file.path(OUT_DIR, sprintf("%s_%s_discrim.rds", db, cc))))
