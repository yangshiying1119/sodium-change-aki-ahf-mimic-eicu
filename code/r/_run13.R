## =============================================================================
## _run13.R —— 阴性对照暴露：把 ΔNa 换成 ΔK / ΔCl 重跑同一套模型
##
## 目的不是找另一个"阳性结果"，而是一个证伪检验：
##   · 若与钠无直接关系的电解质也给出同等量级的"保护"效应，说明观测到的
##     不是钠本身，而是"入院 24 h 内化验整体回归正常 / 病情好转"的伴随现象；
##   · 若只有钠阳性、钾与氯接近 1，则说明效应有一定特异性。
##
## 口径与主分析完全一致：同一个 0–24 h 暴露窗、同一个 delta_na() 构造
## （末次减首次）、同一组 BASIC_X、同一批完整病例。唯一的区别是暴露列。
## 钾 / 氯由 20_labs_*.sql 提取，**不进任何主分析或主模型**。
## =============================================================================
source("extract.R")
source("common.R")
source("models.R")
source("advance.R")

#' 只提取钾 / 氯两张长表（其余化验与队列结果复用已落盘产物）
extract_electrolytes <- function(db, root) {
  con <- dbConnect(duckdb())
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
  if (db == "mimic") {
    run_sql_file(con, "00_schema_mimic.sql", list(MIMIC_ROOT = root))
    run_sql_file(con, "20_labs_mimic.sql")
  } else {
    run_sql_file(con, "00_schema_eicu.sql", list(EICU_ROOT = root))
    run_sql_file(con, "20_labs_eicu.sql")
  }
  long <- function(t) {
    d <- dbGetQuery(con, sprintf("SELECT id, t_hours, value FROM %s", t))
    d$id <- as.character(d$id); d
  }
  list(potassium = long("potassium"), chloride = long("chloride"))
}

#' 用同一套 delta_na() 构造某个电解质的 0–24 h 变化，并把列名改成指定名
electrolyte_delta <- function(dt, ids, newname) {
  sub <- dt[dt$id %in% ids, , drop = FALSE]
  r <- delta_na(sub)
  if (is.null(r)) return(NULL)
  names(r)[names(r) == "delta_na"] <- newname
  r
}

#' 基础模型下某个暴露的效应量
exposure_effect <- function(d, expo) {
  xs <- c(expo, .pick(BASIC_X, d))
  m  <- glm(.fml("aki", xs), data = d, family = binomial())
  co <- summary(m)$coefficients
  if (!expo %in% rownames(co)) return(NULL)
  r  <- co[expo, ]
  data.frame(暴露 = expo, n = nobs(m), 事件数 = sum(d$aki),
             OR = exp(r[1]), 下限 = exp(r[1] - 1.96 * r[2]),
             上限 = exp(r[1] + 1.96 * r[2]), P = r[4], row.names = NULL)
}

DB <- list(
  list(名 = "MIMIC-IV", db = "mimic", root = MIMIC_ROOT,
       full = "mimic_B_full.rds", out = "mimic_negctl.rds"),
  list(名 = "eICU-CRD", db = "eicu",  root = EICU_ROOT,
       full = "eicu_C_full.rds", out = "eicu_C_negctl.rds")
)

for (s in DB) {
  cat("\n================ ", s$名, " ================\n", sep = "")
  a <- readRDS(file.path(OUT_DIR, s$full))$analysis
  a <- a[model_cases(a), , drop = FALSE]

  ee <- extract_electrolytes(s$db, s$root)
  cat(sprintf("钾 %d 行 / %d 个 stay；氯 %d 行 / %d 个 stay\n",
              nrow(ee$potassium), length(unique(ee$potassium$id)),
              nrow(ee$chloride),  length(unique(ee$chloride$id))))

  dk <- electrolyte_delta(ee$potassium, a$id, "delta_k")
  dc <- electrolyte_delta(ee$chloride,  a$id, "delta_cl")
  a$delta_k  <- dk$delta_k[match(a$id, dk$id)]
  a$delta_cl <- dc$delta_cl[match(a$id, dc$id)]

  ## 钠、钾、氯三者的变化在同一个 24 h 窗内有多相关——这是解释阴性对照的前提：
  ## 若相关性很低，三者本就是不同的信息，阴性对照的说服力更强。
  rho <- cor(a[, c("delta_na", "delta_k", "delta_cl")], use = "pairwise.complete.obs")
  cat("\n--- ΔNa / ΔK / ΔCl 的相关系数 ---\n"); print(round(rho, 3))
  cat(sprintf("有 ≥2 次钾者 %d 例，有 ≥2 次氯者 %d 例（分析人群 %d）\n",
              sum(!is.na(a$delta_k)), sum(!is.na(a$delta_cl)), nrow(a)))

  res <- lapply(c("delta_na", "delta_k", "delta_cl"), function(v) {
    d <- a[stats::complete.cases(a[, c("aki", v, .pick(BASIC_X, a)), drop = FALSE]), ]
    exposure_effect(d, v)
  })
  res <- do.call(rbind, res)

  ## 互相调整：钠与钾同时进模型。若钠的效应在调整钾之后基本不变，
  ## 说明两者携带的信息并不重合，阴性对照不构成对主结论的否定。
  m2 <- glm(.fml("aki", c("delta_na", "delta_k", .pick(BASIC_X, a))),
            data = a, family = binomial())
  c2 <- summary(m2)$coefficients
  res <- rbind(res, data.frame(
    暴露 = "delta_na（同时调整 delta_k）", n = nobs(m2), 事件数 = sum(a$aki),
    OR = exp(c2["delta_na", 1]),
    下限 = exp(c2["delta_na", 1] - 1.96 * c2["delta_na", 2]),
    上限 = exp(c2["delta_na", 1] + 1.96 * c2["delta_na", 2]),
    P = c2["delta_na", 4], row.names = NULL))

  ## ---- 钠与氯互调：正面回答"效应是否属于钠，而不是同期变化的氯" ----
  ## 两者同属一个 0–24 h 暴露窗且 r ≈ 0.6，只有放进同一个模型才能比较。
  ## 若钠保留而氯归零，说明氯的关联可由钠解释，反之不成立。
  ## 两种互调（每 1 mmol/L 与每 1 SD）共用同一批完整病例，n 与事件数才一致
  dj <- a[stats::complete.cases(
            a[, c("aki", "delta_na", "delta_cl", .pick(BASIC_X, a)), drop = FALSE]), ]
  m3 <- glm(.fml("aki", c("delta_na", "delta_cl", .pick(BASIC_X, dj))),
            data = dj, family = binomial())
  c3 <- summary(m3)$coefficients
  one3 <- function(v, lab) data.frame(
    暴露 = lab, n = nobs(m3), 事件数 = sum(dj$aki),
    OR = exp(c3[v, 1]),
    下限 = exp(c3[v, 1] - 1.96 * c3[v, 2]),
    上限 = exp(c3[v, 1] + 1.96 * c3[v, 2]),
    P = c3[v, 4], row.names = NULL)
  res <- rbind(res,
               one3("delta_na", "delta_na（同时调整 delta_cl）"),
               one3("delta_cl", "delta_cl（同时调整 delta_na）"))

  ## 按各库自身的标准差标准化后互调：两个化验的离散度不同，每 1 mmol/L 的
  ## OR 不能直接比大小，每 1 SD 才可以——这是回应"两者量级几乎相同"的正确尺度。
  sd_na <- stats::sd(dj$delta_na); sd_cl <- stats::sd(dj$delta_cl)
  dj$z_na <- dj$delta_na / sd_na; dj$z_cl <- dj$delta_cl / sd_cl
  m4 <- glm(.fml("aki", c("z_na", "z_cl", .pick(BASIC_X, dj))),
            data = dj, family = binomial())
  c4 <- summary(m4)$coefficients
  one4 <- function(v, lab) data.frame(
    暴露 = lab, n = nobs(m4), 事件数 = sum(dj$aki),
    OR = exp(c4[v, 1]),
    下限 = exp(c4[v, 1] - 1.96 * c4[v, 2]),
    上限 = exp(c4[v, 1] + 1.96 * c4[v, 2]),
    P = c4[v, 4], row.names = NULL)
  res <- rbind(res,
               one4("z_na", "delta_na（每 1 SD，调整 delta_cl）"),
               one4("z_cl", "delta_cl（每 1 SD，调整 delta_na）"))

  ## 两变量互调的 VIF 只需看相关系数；>5 才需要担心共线
  r_na_cl <- stats::cor(dj$delta_na, dj$delta_cl)
  vif <- 1 / (1 - r_na_cl^2)

  cat("\n--- 阴性对照（暴露：0–24 h 末次减首次，基础模型）---\n")
  for (i in seq_len(nrow(res)))
    cat(sprintf("  %-40s %4d/%4d  OR %.3f（%.3f–%.3f）  P = %s\n",
                res$暴露[i], res$事件数[i], res$n[i],
                res$OR[i], res$下限[i], res$上限[i],
                format.pval(res$P[i], digits = 2, eps = 1e-4)))
  cat(sprintf("  SD：ΔNa %.2f / ΔCl %.2f   r = %.3f   VIF = %.2f\n",
              sd_na, sd_cl, r_na_cl, vif))

  saveRDS(list(效应 = res, 相关 = rho,
               非缺失 = c(钾 = sum(!is.na(a$delta_k)),
                          氯 = sum(!is.na(a$delta_cl))),
               共线 = list(r = r_na_cl, VIF = vif,
                           SD_Na = sd_na, SD_Cl = sd_cl),
               n = nrow(a), 事件数 = sum(a$aki)),
          file.path(OUT_DIR, s$out))
  cat(sprintf("\n已存 %s\n", file.path(OUT_DIR, s$out)))
}
