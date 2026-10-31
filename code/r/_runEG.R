## _runEG.R —— eICU 医院聚类稳健敏感性（外审 P0；数据到货后运行）
## 用法：Rscript _runEG.R   （须先跑完 SQL 提取层，eicu_C_full.rds 含 hospitalid）
## 输出：out/eicu_cluster.csv（普通 vs 聚类稳健 SE，基础/扩展两档）+ 医院概览
source("extract.R")
source("models.R")
source("cluster_eicu.R")

full <- readRDS(file.path(OUT_DIR, "eicu_C_full.rds"))
a <- full$analysis
ok <- model_cases(a) & !is.na(a$hospitalid)
d <- a[ok, ]
cat(sprintf("建模 n=%d, 事件=%d, 医院=%d\n", nrow(d), sum(d$aki),
            length(unique(d$hospitalid))))

res <- rbind(
  fit_cluster_or(d, .pick(BASIC_X, d), label = "基础"),
  fit_cluster_or(d, c(.pick(BASIC_X, d), .pick(EXTEND_X, d)), label = "扩展")
)
ov <- hospital_overview(d)
cat("\n医院概览：\n"); str(ov)

write.csv(res, file.path(OUT_DIR, "eicu_cluster.csv"), row.names = FALSE)
print(res, digits = 4)
cat("\nSAVED", file.path(OUT_DIR, "eicu_cluster.csv"), "\n")
