## eICU 重复队列的论文表导出。与 MIMIC 侧 _run6.R 同构：只读已落盘的结果对象，
## 不重跑任何模型；两库走的是**同一套** export_db()，只有落盘目录与工作簿名不同。
##
## 用法：Rscript _runEC.R [口径字母]   缺省用 C（见 _run9.R 中 EICU_CC 的说明）
source("extract.R")
source("common.R")      # 界标/结局窗常数（competing.R 的默认参数要用）
source("cohort.R")      # attach_covariates（基线表要并协变量）
source("models.R")      # subgroup_labels()：亚组表要用与图 11 同一套显示名
source("tables.R")
source("competing.R")
source("export.R")

args <- commandArgs(trailingOnly = TRUE)
CC   <- if (length(args) >= 1 && nzchar(args[1])) toupper(args[1]) else "C"

full <- readRDS(file.path(OUT_DIR, sprintf("eicu_%s_full.rds", CC)))
mods <- readRDS(file.path(OUT_DIR, sprintf("eicu_%s_models.rds", CC)))
sens <- readRDS(file.path(OUT_DIR, sprintf("eicu_%s_sensitivity.rds", CC)))
cmp  <- readRDS(file.path(OUT_DIR, sprintf("eicu_%s_competing.rds", CC)))

## AJ 汇总按 MIMIC 侧同样的理由重算：风险集数据已在 rds 内，只是重读，不重跑模型。
cmp$aj_summary <- fit_aalen_johansen(cmp$competing_data)$汇总
cat("AJ 累积发生率（72 h）：\n"); print(cmp$aj_summary)

## S8 复合口径与偏倚分析在 eICU 侧尚未跑（见正文跨库差异一节），
## 在的话就一并出表；不在则表 5b、表 8–10 不出现，不报错。
.pick_rds <- function(f) if (file.exists(f)) readRDS(f) else NULL
s8   <- .pick_rds(file.path(OUT_DIR, sprintf("eicu_%s_s8_uo.rds", CC)))
bias <- .pick_rds(file.path(OUT_DIR, sprintf("eicu_%s_bias.rds", CC)))

sens_extra <- NULL
if (!is.null(s8))   sens_extra <- rbind(sens_extra, s8$行)
if (!is.null(bias)) sens_extra <- rbind(sens_extra, bias$模型)

extra_tabs <- NULL
if (!is.null(s8)) extra_tabs <- list(表5b_肌酐与尿量标准交叉 = s8$交叉)

## 主分析那一行的分析名带上表型口径，避免与 MIMIC 的表对照时看不出是哪一版口径
sens$分析 <- sub("^主分析.*$", sprintf("主分析（表型 %s）", CC), sens$分析)

tabs <- export_eicu(full, mods, sens, cmp,
                    sens_extra = sens_extra, bias = bias,
                    extra_tabs = extra_tabs)

cat("\n=== 已导出表（eICU，表型", CC, "）===\n")
for (nm in names(tabs))
  cat(sprintf("  %-24s %3d 行 × %2d 列\n", nm, nrow(tabs[[nm]]), ncol(tabs[[nm]])))
tab_dir <- file.path(OUT_DIR, out_subdir("tables_eicu", "tables_eicu_en"))
cat("\n表件目录：", tab_dir, "\n")
for (f in list.files(tab_dir, full.names = TRUE))
  cat(sprintf("  %-56s %8.1f KB\n", basename(f), file.size(f) / 1024))
