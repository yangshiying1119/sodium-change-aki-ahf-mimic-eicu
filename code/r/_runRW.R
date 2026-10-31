## _runRW.R —— 限定 72h 内未转出（结局窗结束时仍在 ICU）子集的敏感性分析
## 动机：主分析结局窗 24-72h 的肌酐取自全院检验表，转出后病房肌酐仍参与判定
## （MIMIC 134/706=19.0%、eICU 54/490=11.0% 的 AKI 系转出后检出）。
## 外审问：主结果是否被这类"转出后检出"撑起来？本脚本把 72h 前出科/死亡者
## 整批剔除，只在结局窗结束时仍在 ICU 的子集上重跑基础/扩展模型。
## 与 competing.R 的 Fine-Gray 互补：那边把转出当竞争事件，这边直接限人群。
.args <- commandArgs(FALSE)
.sf <- sub("^--file=", "", grep("^--file=", .args, value = TRUE))
if (length(.sf)) setwd(dirname(normalizePath(.sf[1])))
suppressMessages({ source("extract.R"); source("common.R"); source("models.R") })
## 复现缓存目录：默认与主管道一致（extract.R 的 OUT_DIR）；换机器可设 SODIUM_OUT 覆盖
OUT <- Sys.getenv("SODIUM_OUT", OUT_DIR)

fmt_or <- function(co, dig = 3) sprintf("%.3f (%.3f-%.3f), P=%s",
  exp(co[1]), exp(co[1]-1.96*co[2]), exp(co[1]+1.96*co[2]),
  format.pval(co[4], digits=3))

rows <- list()
for (pair in list(c("MIMIC-IV", "mimic_B_full.rds"),
                  c("eICU-CRD", "eicu_C_full.rds"))) {
  tag <- pair[1]
  a <- readRDS(file.path(OUT, pair[2]))$analysis
  d <- a[model_cases(a), ]
  stay72 <- as.numeric(difftime(d$outtime, d$intime, units = "hours")) > 72
  cat(sprintf("\n========== %s ==========\n", tag))
  cat(sprintf("完整病例 n=%d, AKI=%d；72h 结局窗结束时仍在 ICU：%d (%.1f%%)，其中 AKI=%d\n",
      nrow(d), sum(d$aki), sum(stay72), 100*mean(stay72), sum(d$aki[stay72])))
  cat(sprintf("72h 前出科/死亡被剔除：%d，其中 AKI（出科前检出）=%d\n",
      sum(!stay72), sum(d$aki[!stay72])))
  dd <- d[stay72, ]
  for (tier in c("base", "extended")) {
    xs <- if (tier == "base") .pick(BASIC_X, dd)
          else c(.pick(BASIC_X, dd), .pick(EXTEND_X, dd))
    m <- glm(.fml("aki", c("delta_na", xs)), data = dd, family = binomial())
    co <- summary(m)$coeff["delta_na", ]
    cat(sprintf("  [%s] 72h在科子集 (n=%d, ev=%d): ΔNa %s\n", tier, nrow(dd),
        sum(dd$aki), fmt_or(co)))
    rows[[length(rows)+1]] <- data.frame(
      Database = tag, Model = if (tier == "base") "Base" else "Extended",
      `n (events)` = sprintf("%d (%d)", nrow(dd), sum(dd$aki)),
      `OR (95% CI)` = sprintf("%.3f (%.3f-%.3f)", exp(co[1]),
        exp(co[1]-1.96*co[2]), exp(co[1]+1.96*co[2])),
      `P value` = format.pval(co[4], digits = 3), check.names = FALSE,
      row.names = NULL)
  }
}
res <- do.call(rbind, rows)
write.csv(res, file.path(OUT, "stay72_sensitivity.csv"), row.names = FALSE)
cat("\nSAVED", file.path(OUT, "stay72_sensitivity.csv"), "\n")
