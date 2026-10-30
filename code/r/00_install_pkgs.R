# 一次性环境准备：安装本方案所需但机器上缺失的 R 包
# R 4.1.2 较老，CRAN 最新版可能要求更高 R，失败时回退到 archive 版本
options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 16)

need <- c("DBI", "duckdb", "mice", "cmprsk", "pROC", "MASS", "jsonlite", "ragg", "forestploter", "openxlsx")
have <- rownames(installed.packages())
todo <- setdiff(need, have)
cat("待安装:", if (length(todo)) paste(todo, collapse=", ") else "(无)", "\n")

install_one <- function(p) {
  ok <- tryCatch({ install.packages(p); requireNamespace(p, quietly=TRUE) },
                 error = function(e) FALSE, warning = function(w) FALSE)
  if (!ok) {
    cat("  最新版失败，回退 archive:", p, "\n")
    tryCatch(install.packages(p, repos="https://cran.r-project.org/src/contrib/Archive/"),
             error=function(e) cat("  archive 亦失败:", conditionMessage(e), "\n"))
  }
  requireNamespace(p, quietly=TRUE)
}

for (p in todo) {
  cat("=== 安装", p, "===\n")
  cat("  结果:", if (install_one(p)) "OK" else "失败", "\n")
}

cat("\n=== 最终核对 ===\n")
for (p in need) cat(sprintf("  %-10s %s\n", p, if (requireNamespace(p, quietly=TRUE)) "OK" else "仍缺失"))
