## =============================================================================
## tests/run_tests.R —— 测试入口（仓库根目录运行：Rscript tests/run_tests.R）
## 断言全部在 code/r/test_common.R 中维护（共 69 条）。common.R/cohort.R 等被测
## 文件用 --file= 参数自行定位同目录文件，因此本入口以子进程方式委派执行，
## 保证被测脚本的定位逻辑与其在分析机上直接运行时完全一致；子进程退出码
## 原样透传（非 0 即有断言失败）。
## =============================================================================
.args <- commandArgs(trailingOnly = FALSE)
.file <- sub("^--file=", "", .args[grep("^--file=", .args)])
.here <- if (length(.file) && nzchar(.file)) dirname(normalizePath(.file)) else "."

.target <- file.path(.here, "..", "code", "r", "test_common.R")
if (!file.exists(.target)) stop("未找到 code/r/test_common.R；请在仓库根目录运行 Rscript tests/run_tests.R")
.target <- normalizePath(.target)

## 被测脚本内部以裸相对路径 source 同目录文件（如 bias.R 的 source("cohort.R")），
## 故子进程工作目录须为 code/r/；test_common.R 自身经 --file= 绝对路径自定位，不受 cwd 影响。
setwd(dirname(.target))

Rbin <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
status <- system2(Rbin, shQuote(.target))
quit(status = if (is.na(status)) 1L else status)
