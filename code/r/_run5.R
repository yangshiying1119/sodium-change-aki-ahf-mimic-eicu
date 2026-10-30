## 出图：读已落盘的结果，画图 1–5（不重跑模型）
source("extract.R")
source("figures.R")

fig_dir <- file.path(OUT_DIR, out_subdir("figures", "figures_en"))
full <- readRDS(file.path(OUT_DIR, "mimic_B_full.rds"))
mods <- readRDS(file.path(OUT_DIR, "mimic_B_models.rds"))
sens <- readRDS(file.path(OUT_DIR, "mimic_sensitivity.rds"))
cmp  <- readRDS(file.path(OUT_DIR, "mimic_competing.rds"))

## ---- 图 1 流程图 ----
fig_flow(full$attrition, sum(full$analysis$aki), cohort = "B", dir = fig_dir)
cat("图 1 已出\n")

## ---- 图 2 RCS ----
fig_rcs(mods$curve, mods$rcs, dir = fig_dir)
cat("图 2 已出\n")

## ---- 图 3 敏感性分析 ----
s <- sens
s$分析 <- sub("^主分析.*$", "主分析", s$分析)
fig_forest(s, "分析", "图 3 敏感性分析：ΔNa 每升高 1 mmol/L 的 OR",
           note = "全部采用基础模型调整；\n红色方块为主分析。",
           ref_row = "主分析", dir = fig_dir, file = "fig3_sensitivity.png")
cat("图 3 已出\n")

## ---- 图 4 预设亚组 ----
sg <- mods$subgroups
sg <- sg[!is.na(sg$OR), , drop = FALSE]
## 标签一律走 subgroup_labels()，不要在驱动脚本里就地拼字符串：
## 早先按 "vs" 切亚组名，把 "<135 组" 标成了 "≥135 组"，方向标反。
sg$标签 <- subgroup_labels(sg$亚组, sg$层)
## 表头不能只写"亚组"：亚组图的左列很短，forestploter 的脚注又从表格左边缘起排，
## 窄左列会把脚注顶到刻度标签上去。表头带一句说明既撑宽左列、也点明了建模方式。
## 脚注同理只能短——"交互 P 见正文表"这类交代放到图题（build_manuscript.py）里，
## 硬塞进图内会把 0.90 那个刻度压住。
fig_forest(sg, "标签", "图 4 预设亚组分析（探索性）",
           note = "各层内均为基础模型调整。",
           dir = fig_dir, file = "fig4_subgroup.png",
           header = "亚组（各层独立建模）")
cat("图 4 已出\n")

## ---- 图 5 累积发生率 ----
fig_aj(cmp$competing_data, dir = fig_dir)
cat("图 5 已出\n")

cat("\n图件目录：", fig_dir, "\n")
for (f in list.files(fig_dir, full.names = TRUE))
  cat(sprintf("  %-28s %8.0f KB\n", basename(f), file.size(f) / 1024))
