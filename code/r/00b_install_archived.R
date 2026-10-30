# R 4.1.2 兼容版安装：CRAN 最新版要求更高 R，改从 Archive 取旧版本
options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 16)
A <- "https://cran.r-project.org/src/contrib/Archive"

try_install <- function(pkg, ver, url) {
  if (requireNamespace(pkg, quietly = TRUE)) { cat("  已存在，跳过:", pkg, "\n"); return(TRUE) }
  cat("=== 安装", pkg, ver, "===\n")
  ok <- tryCatch({ install.packages(url, repos = NULL, type = "source"); TRUE },
                 error = function(e) { cat("  失败:", conditionMessage(e), "\n"); FALSE })
  requireNamespace(pkg, quietly = TRUE)
}

# 仅 duckdb 需要钉版（SQL 提取层引擎，版本与 sessionInfo.txt 一致）。
# 分析代码不调用 rms/Hmisc（样条为 models.R 手写基函数），无需安装。
try_install("duckdb","1.2.0", file.path(A, "duckdb","duckdb_1.2.0.tar.gz"))

cat("\n=== 核对 ===\n")
for (p in c("duckdb", "mice", "cmprsk", "pROC"))
  cat(sprintf("  %-10s %s\n", p, if (requireNamespace(p, quietly = TRUE)) "OK" else "仍缺失"))
