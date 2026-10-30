## =============================================================================
## _runEInline.R —— 把"在 build_manuscript.py 里就地拼出"的四张表翻成英文
##
## 表 3c / 3d / 14 / 16 不来自 R 侧的 CSV：它们由 results.json 直接拼出，因为
## 内容跨库或混用多个结果段，落 CSV 反而多一处可能不同步的副本。但英文版需要
## 同一批内容，所以 build_manuscript.py 在排版时把**与文档中完全相同的
## header/rows** 顺手写到 out/tables_inline/（同一对象写盘，不存在第二份真相）。
##
## 本脚本只做翻译与分发，不重算任何数字。
## =============================================================================

.here_in <- local({
  a <- commandArgs(trailingOnly = FALSE)
  s <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(s) && nzchar(s)) dirname(normalizePath(s)) else "."
})
source(file.path(.here_in, "extract.R"))   # OUT_DIR
source(file.path(.here_in, "common.R"))    # i18n 层（is_en / tr_df）

src <- file.path(OUT_DIR, "tables_inline")
if (!dir.exists(src)) {
  stop("找不到 ", src, "：先跑 build_manuscript.py，它会把就地拼出的表落在这里")
}

## 中文名 → 英文表名；落盘目录按内容归属分配：3c/3d 属主队列，14 属两库对照，
## 16 属 eICU 重复队列。英文模式下才写文件，中文版仍由 build_manuscript.py
## 就地排版，不产生多余副本。
PLAN <- list(
  "表3c_判别与校准"   = list(en = "Table3c_Discrimination",        dir = "tables_en"),
  "表3d_校准曲线"     = list(en = "Table3d_CalibrationCurve",      dir = "tables_en"),
  "表14_两库三档并排" = list(en = "Table14_TiersSideBySide",        dir = "tables_cross_en"),
  "表16_eICU严格表型" = list(en = "Table16_eICU_StrictPhenotype",   dir = "tables_eicu_en")
)

## 与 export.R 的 .write_csv() 同一写法：UTF-8 + BOM，Excel 双击不乱码
.write_csv <- function(d, name, dir) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  con <- file(file.path(dir, paste0(name, ".csv")), open = "wb")
  on.exit(close(con), add = TRUE)
  writeBin(charToRaw("\xEF\xBB\xBF"), con)
  utils::write.csv(d, con, row.names = FALSE, fileEncoding = "UTF-8")
  invisible(file.path(dir, paste0(name, ".csv")))
}

## 中文模式下这几张表已经在 docx 里排好了，不另写 CSV——否则 out/tables/ 会多出
## 四份没有人读的副本，还与 export_db 的产物混在一处。
if (!is_en()) {
  cat("OUT_LANG=zh：表 3c/3d/14/16 由 build_manuscript.py 直接排版，无需另存 CSV\n")
  quit(save = "no")
}

cat("就地表的英文版：\n")
for (zh in names(PLAN)) {
  fp <- file.path(src, paste0(zh, ".csv"))
  if (!file.exists(fp)) {
    cat(sprintf("  %-20s 未找到，跳过\n", zh))
    next
  }
  ## 全列读成字符：这几张表里的数字是**已经格式化好的呈现串**（"24.0"、"2,656"）。
  ## 若让 read.csv 猜类型，24.0 会被解析成数值、写回时变成 "24"，中文版印的是
  ## 24.0、英文版印 24，同一张表两副面孔。
  d <- utils::read.csv(fp, encoding = "UTF-8", check.names = FALSE,
                       stringsAsFactors = FALSE, colClasses = "character")
  p <- PLAN[[zh]]
  d <- tr_df(d)
  out_dir <- file.path(OUT_DIR, p$dir)
  path <- .write_csv(d, p$en, out_dir)
  cat(sprintf("  %-40s -> %s\n", basename(path), dirname(path)))
}
