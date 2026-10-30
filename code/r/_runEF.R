## 出图（eICU 重复队列 + 双库并排）。只读已落盘的结果，不重跑模型。
##
## 用法：Rscript _runEF.R [eICU 口径字母]，缺省 C
source("extract.R")
source("figures.R")
source("figures_cross.R")

args <- commandArgs(trailingOnly = TRUE)
CC   <- if (length(args) >= 1 && nzchar(args[1])) toupper(args[1]) else "C"

fig_e <- file.path(OUT_DIR, out_subdir("figures_eicu", "figures_eicu_en"))
fig_x <- file.path(OUT_DIR, out_subdir("figures_cross", "figures_cross_en"))

full <- readRDS(file.path(OUT_DIR, sprintf("eicu_%s_full.rds", CC)))
mods <- readRDS(file.path(OUT_DIR, sprintf("eicu_%s_models.rds", CC)))
cmp  <- readRDS(file.path(OUT_DIR, sprintf("eicu_%s_competing.rds", CC)))
sens <- readRDS(file.path(OUT_DIR, sprintf("eicu_%s_sensitivity.rds", CC)))

## ---- eICU 单库图（与主队列同构，文件名加 eicu_ 前缀）----
fig_flow(full$attrition, sum(full$analysis$aki), cohort = CC,
         dir = fig_e, file = "eicu_fig1_flow.png")
cat("eICU 图 1 已出\n")

fig_rcs(mods$curve, mods$rcs, dir = fig_e, file = "eicu_fig2_rcs.png")
cat("eICU 图 2 已出\n")

s <- sens
s$分析 <- sub("^主分析.*$", sprintf("主分析（表型 %s）", CC), s$分析)
fig_forest(s, "分析", sprintf("敏感性分析：ΔNa 每升高 1 mmol/L 的 OR（eICU，表型 %s）", CC),
           note = "全部采用基础模型调整；红色方块为主分析。\n本库无血气钠，故无 S4/S4b。",
           ref_row = s$分析[grepl("^主分析", s$分析)],
           dir = fig_e, file = "eicu_fig3_sensitivity.png")
cat("eICU 图 3 已出\n")

## 与主队列同源：标签走 subgroup_labels()，不要直接把"亚组"列当行标签——
## 该列每个亚组重复两行（"是"/"否"各一行），直接拿来会画成两个一模一样的标签。
sg <- mods$subgroups[!is.na(mods$subgroups$OR), , drop = FALSE]
sg$标签 <- subgroup_labels(sg$亚组, sg$层)
fig_forest(sg, "标签", "预设亚组分析（eICU，探索性）",
           note = "各层内均为基础模型调整。",
           dir = fig_e, file = "eicu_fig4_subgroup.png",
           header = "亚组（各层独立建模）")

cat("eICU 图 4 已出\n")

fig_aj(cmp$competing_data, dir = fig_e, file = "eicu_fig5_aj.png")
cat("eICU 图 5 已出\n")

## ---- 双库并排图 ----
## MIMIC 侧的曲线与模型读主队列对象；两库都取"主分析口径"的那一份结果，
## 不做任何重新标定。
mm   <- readRDS(file.path(OUT_DIR, "mimic_B_models.rds"))
fig_rcs_cross(mm$curve, mods$curve, mm$rcs, mods$rcs,
              lab_m = "MIMIC-IV 主队列（表型 B）",
              lab_e = sprintf("eICU-CRD 重复队列（表型 %s）", CC),
              dir = fig_x, file = "fig6_rcs_cross.png")
cat("图 6 双库剂量反应叠加已出\n")

fig_tiers_cross(mm$tiers, mods$tiers, dir = fig_x, file = "fig7_tiers_cross.png")
cat("图 7 双库并排森林图已出\n")

cat("\n图件目录：\n")
for (d in c(fig_e, fig_x)) {
  cat(" ", d, "\n")
  for (f in list.files(d, full.names = TRUE))
    cat(sprintf("    %-34s %8.0f KB\n", basename(f), file.size(f) / 1024))
}
