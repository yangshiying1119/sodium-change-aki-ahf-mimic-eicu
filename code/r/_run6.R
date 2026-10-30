## 出表：读已落盘的结果，导出论文表（不重跑模型）
source("extract.R")
source("common.R")      # 界标/结局窗等常数（competing.R 的默认参数要用）
source("cohort.R")      # attach_covariates（基线表要并协变量）
source("models.R")      # subgroup_labels()：亚组表要用与图 4 同一套显示名
source("tables.R")
source("competing.R")
source("export.R")

full <- readRDS(file.path(OUT_DIR, "mimic_B_full.rds"))
mods <- readRDS(file.path(OUT_DIR, "mimic_B_models.rds"))
sens <- readRDS(file.path(OUT_DIR, "mimic_sensitivity.rds"))
cmp  <- readRDS(file.path(OUT_DIR, "mimic_competing.rds"))

## AJ 汇总在建表前重算一次：早先 fit_aalen_johansen 的汇总行是把 1×3 概率矩阵
## 按列回收出来的，三个状态记成了同一组数（详见 competing.R 的注释）。
## 风险集数据已在 rds 里，重算只是读取，不需要重跑模型。
cmp$aj_summary <- fit_aalen_johansen(cmp$competing_data)$汇总
cat("AJ 累积发生率（72 h）：\n"); print(cmp$aj_summary)

## 偏倚分析与 S8 各自单独跑（_run7.R），本步只负责"在的话就一起出表"。
## 这样即便偏倚部分还没跑，主表也能先出来；顺序不构成依赖。
.pick_rds <- function(f) if (file.exists(f)) readRDS(f) else NULL
s8   <- .pick_rds(file.path(OUT_DIR, "mimic_s8_uo.rds"))
bias <- .pick_rds(file.path(OUT_DIR, "mimic_bias.rds"))

## 追加到表 5 的行：S8 复合口径 + 复测概率加权。两者都是"换口径看结论稳不稳"，
## 与 S1–S11 同族，故并入同一张表而非另起。
sens_extra <- NULL
if (!is.null(s8))   sens_extra <- rbind(sens_extra, s8$行)
if (!is.null(bias)) sens_extra <- rbind(sens_extra, bias$模型)
if (is.null(sens_extra))
  cat("提示：尚未找到 S8 / 偏倚分析结果，表 5 只含 S1–S5，表 8–10 不导出\n")

## S8 的交叉表也要留档：仅肌酐检出 vs 复合口径检出的重叠，
## 是"尿量标准补回了多少事件"的直接证据。
extra_tabs <- NULL
if (!is.null(s8)) extra_tabs <- list(表5b_肌酐与尿量标准交叉 = s8$交叉)

tabs <- export_mimic(full, mods, sens, cmp,
                     sens_extra = sens_extra, bias = bias,
                     extra_tabs = extra_tabs)

cat("\n=== 已导出表 ===\n")
for (nm in names(tabs))
  cat(sprintf("  %-24s %3d 行 × %2d 列\n", nm, nrow(tabs[[nm]]), ncol(tabs[[nm]])))
tab_dir <- file.path(OUT_DIR, out_subdir("tables", "tables_en"))
cat("\n表件目录：", tab_dir, "\n")
for (f in list.files(tab_dir, full.names = TRUE))
  cat(sprintf("  %-56s %8.1f KB\n", basename(f), file.size(f) / 1024))
