## =============================================================================
## test_common.R —— common.R 的单元测试
## 构造已知答案的序列，断言判定输出。重点覆盖方案最容易被做错的两点：
##   (1) 比较值绝不取自事件之后（"事后更低值"陷阱）
##   (2) 界标边界：暴露窗 [0,24] 闭、结局窗 (24,72] 左开右闭
## 运行：Rscript test_common.R   （退出码非 0 表示有断言失败）
## =============================================================================

.args <- commandArgs(trailingOnly = FALSE)
.script <- sub("^--file=", "", .args[grep("^--file=", .args)])
.here <- if (length(.script) && nzchar(.script)) dirname(normalizePath(.script)) else "."
source(file.path(.here, "common.R"))

.pass <- 0L; .fail <- 0L
ok <- function(cond, label) {
  if (isTRUE(cond)) { .pass <<- .pass + 1L; cat("  PASS ", label, "\n") }
  else { .fail <<- .fail + 1L; cat("  FAIL ", label, "\n") }
}
near <- function(a, b, tol = 1e-9) !is.na(a) && !is.na(b) && abs(a - b) < tol
L <- function(id, t, v) data.frame(id = id, t_hours = t, value = v, stringsAsFactors = FALSE)

cat("\n== 1. 清洗边界 ==\n")
ok(is.na(clean_sodium(89))  && near(clean_sodium(90), 90) &&
   near(clean_sodium(200), 200) && is.na(clean_sodium(201)), "血钠 90–200 闭区间")
ok(is.na(clean_creat(0)) && near(clean_creat(0.1), 0.1) &&
   near(clean_creat(50), 50) && is.na(clean_creat(808)), "肌酐 0.1–50 闭区间")

cat("\n== 2. 同时间重复合并：奇数取中位数、偶数取较低者 ==\n")
m <- merge_same_time(L("a", c(1, 1.1, 1.2), c(140, 136, 138)))   # 3 条 → 中位数 138
ok(nrow(m) == 1L && near(m$value, 138), "3 条取中位数 = 138")
m <- merge_same_time(L("a", c(1, 1.1), c(140, 138)))             # 2 条 → 较低者 138
ok(nrow(m) == 1L && near(m$value, 138), "2 条取较低者 = 138（不是均值 139）")
m <- merge_same_time(L("a", c(1, 2), c(140, 136)))               # 间隔 1h > 15min → 不合并
ok(nrow(m) == 2L, "间隔 1 小时不合并")

cat("\n== 3. ΔNa = 0–24h 末次 − 首次 ==\n")
d <- delta_na(L("a", c(0, 12, 24), c(138, 140, 134)))
ok(nrow(d) == 1L && near(d$delta_na, -4) && near(d$na_first, 138) && near(d$na_last, 134),
   "ΔNa = 134 − 138 = −4")
d <- delta_na(L("a", c(0, 24, 30), c(138, 140, 130)))
ok(near(d$delta_na, 2), "t=30 在暴露窗外，不参与 ΔNa")
d <- delta_na(L("a", c(0, 12), c(138, 140)))
ok(nrow(d) == 1L, "2 次即可算 ΔNa")
ok(is.null(delta_na(L("a", 5, 138))), "仅 1 次 → 不可算 ΔNa")

cat("\n== 4. KDIGO：48h 内升 ≥0.3 ==\n")
r <- kdigo_first_event(L("a", c(0, 25, 30), c(1.0, 1.0, 1.4)))
ok(r$aki && near(r$t_event, 30) && r$stage_first == 1L, "t=30 检出 1 期，t_event=30")
ok(!r$aki_pre, "界标前无事件")

cat("\n== 5. 【关键反例】事后更低值不得作基线 ==\n")
## 峰值在 t=30，之后 t=40 掉到 0.5。若错误地把"事件后最低值"当基线，t=30 就检不出。
r <- kdigo_first_event(L("a", c(0, 25, 30, 40), c(1.0, 1.0, 1.4, 0.5)))
ok(r$aki && near(r$t_event, 30), "t=30 仍检出（未被 t=40 的低值抹掉）")
ok(near(r$low48_at_event, 1.0), "事件时点比较基线 = 1.0（严格取自 t<30）")

cat("\n== 6. KDIGO：1.5× 7 天基线（用院前负值）==\n")
## t=−100 的院前实测值落在 7 天窗内、48h 窗外 → 只能由 7 天标准检出
r <- kdigo_first_event(L("a", c(-100, 25), c(0.6, 1.0)))
ok(r$aki && r$stage_first == 1L && near(r$low7_at_event, 0.6),
   "院前实测值参与 7 天基线，t=25 检出")
ok(is.na(r$low48_at_event), "该时点 48h 窗内无比较值（−100 在 48h 之外）")

cat("\n== 7. 界标前已满足 KDIGO → aki_pre ==\n")
r <- kdigo_first_event(L("a", c(0, 10, 30), c(1.0, 1.4, 1.4)))
ok(r$aki_pre, "t=10 已达标 → aki_pre = TRUE（主分析排除）")

cat("\n== 8. 界标边界：暴露窗闭 [0,24]、结局窗开 (24,72] ==\n")
r <- kdigo_first_event(L("a", c(0, 24), c(1.0, 1.4)))
ok(r$aki_pre && !r$aki, "t=24 算界标前（aki_pre），不算结局事件")
r <- kdigo_first_event(L("a", c(0, 50, 72), c(1.0, 1.0, 1.4)))
ok(r$aki && near(r$t_event, 72), "t=72 在结局窗内，计入")
r <- kdigo_first_event(L("a", c(0, 50, 73), c(1.0, 1.0, 1.4)))
ok(!r$aki, "t=73 超出观察窗，不计入")

cat("\n== 9. 平稳序列不误报 ==\n")
r <- kdigo_first_event(L("a", c(0, 6, 12, 24, 36, 48, 72), rep(1.2, 7)))
ok(!r$aki && !r$aki_pre && r$stage_max == 0L, "完全平稳 → 无事件")

cat("\n== 10. 3 期与 RRT ==\n")
r <- kdigo_first_event(L("a", c(0, 25, 36), c(1.0, 1.0, 3.2)))
ok(r$aki && r$stage_first == 3L, "3.2 ≥ 3× 基线 → 3 期")
r <- kdigo_first_event(L("a", c(0, 25), c(1.0, 1.1)), rrt = "a")
ok(r$stage_max == 3L, "界标后 RRT → stage_max 提到 3")

cat("\n== 11. assemble() 的可判定性 ==\n")
cohort <- data.frame(id = c("a", "b", "c"), stringsAsFactors = FALSE)
na_tab <- L(c("a", "b", "c"), c(0, 0, 0), c(140, 140, 140))
na_tab <- rbind(na_tab, L(c("a", "b"), c(24, 24), c(136, 136)))   # c 只有 1 次 → 不可算
cr_tab <- L(c("a", "b", "b"), c(0, 0, 10), c(1.0, 1.0, 1.4))      # b 界标前已达标
res <- assemble(cohort, na_tab, rbind(cr_tab, L("a", 30, 1.4)))
res <- res[order(res$id), ]
ok(!is.na(res$delta_na[res$id == "a"]) && is.na(res$delta_na[res$id == "c"]),
   "c 血钠仅 1 次 → ΔNa 缺失")
ok(!res$evaluable[res$id == "b"], "b 界标前已达标 → 不可判定")
ok(res$evaluable[res$id == "a"] && res$aki[res$id == "a"], "a 可判定且检出事件")
ok(!res$has_delta_na[res$id == "c"], "has_delta_na 正确标记")

cat("\n== 12. 界标前判定窗口 pre_from（方案 4.2.6 原文取 [0,24]）==\n")
## 事件只发生在院前（t=-40），0–24h 内平稳：
##   默认 pre_from=0 → 不算"界标前已达标"，该例可进主分析
##   pre_from=-168 → 算界标前已达标，主分析排除（S 系列敏感性分析用）
pre_only <- L("a", c(-50, -40, 0, 12), c(1.0, 1.4, 1.4, 1.4))
r <- kdigo_first_event(pre_only)
ok(!r$aki_pre, "pre_from=0 时院前事件不计入 aki_pre")
r <- kdigo_first_event(pre_only, pre_from = -PRE_H)
ok(r$aki_pre, "pre_from=-168 时同一例被判为界标前已达标")
ok(kdigo_first_event(L("a", c(-50, -40), c(1.0, 1.4)))$n_creat_pre == 0L &&
   kdigo_first_event(L("a", c(-50, -40, 0), c(1.0, 1.4, 1.4)))$n_creat_pre == 1L,
   "n_creat_pre 只数界标前窗内的检测")
## pre_from 不改变比较基线：院前实测值始终参与 48h/7d 滚动窗
ok(near(kdigo_first_event(L("a", c(-100, 25), c(0.6, 1.0)))$low7_at_event, 0.6),
   "pre_from=0 下院前值仍参与 7 天比较基线")

cat("\n== 13. 浮点边界：恰好 0.30 / 恰好 1.5 倍必须达标（KDIGO 用 ≥）==\n")
## 1.60 + 0.3 在双精度下 = 1.90000000000000013 > 1.90，直接比较会漏检
ok(1.60 + 0.3 > 1.9, "前提确认：1.60 + 0.3 在双精度下确实大于 1.9")
r <- kdigo_first_event(L("a", c(0, 30), c(1.60, 1.90)))
ok(r$aki && r$stage_first == 1L && near(r$t_event, 30),
   "恰好升高 0.30 → 判 1 期（不是因浮点漏检）")
r <- kdigo_first_event(L("a", c(0, 30), c(1.20, 1.80)))
ok(r$aki && r$stage_first == 1L, "恰好 1.5 倍（1.20→1.80）→ 判 1 期")
r <- kdigo_first_event(L("a", c(0, 30), c(1.00, 3.00)))
ok(r$aki && r$stage_first == 3L, "恰好 3 倍 → 判 3 期")
r <- kdigo_first_event(L("a", c(0, 30), c(1.00, 2.00)))
ok(r$aki && r$stage_first == 2L, "恰好 2 倍（1.00→2.00）→ 判 2 期")
## 容差不得制造假阳性：差 0.01 仍不达标
r <- kdigo_first_event(L("a", c(0, 30), c(1.60, 1.89)))
ok(!r$aki, "升高 0.29 → 仍不达标（容差未越界）")

cat("\n== 14. n_creat_before：界标前（含院前）肌酐计数 ==\n")
## 方案 3.4④ 用"界标前可判定"筛人，院前实测值同样是有效比较基线，必须计入
r <- kdigo_first_event(L("a", c(-100, 5, 30), c(1.0, 1.0, 1.4)))
ok(r$n_creat_before == 2L, "t=−100 与 t=5 都算界标前（=2）")
ok(r$n_creat_pre == 1L, "n_creat_pre 只数 [0,24] 内（=1）")
ok(kdigo_first_event(L("a", c(40, 50), c(1.0, 1.4)))$n_creat_before == 0L,
   "界标前无肌酐 → 0（无法判定）")
ok(kdigo_first_event(L("a", c(24, 40), c(1.0, 1.4)))$n_creat_before == 1L,
   "t=24 恰在界标上 → 计入界标前（=1，仍不足以判定）")

cat("\n== 15. S7：48h 绝对升高 与 1.5 倍 两条标准的开关 ==\n")
## 关掉一条后，只由该条触发的事件必须消失，另一条仍照常生效。
## 要单独检验 1.5 倍这条，基线必须低到让 1.5 倍时的绝对升幅 <0.3，
## 否则 48h 规则会同时触发，测不出开关是否生效。
## 基线 1.00→1.50 的绝对升幅是 0.50 ≥ 0.3，两条都会亮，不能用。
## 改用 0.50→0.75（1.5 倍，升幅仅 0.25）。
r <- kdigo_first_event(L("a", c(0, 30), c(0.50, 0.75)))
ok(r$aki, "基准：1.5 倍（0.50→0.75，升幅 0.25）仅由倍率规则达标")
ok(!kdigo_first_event(L("a", c(0, 30), c(0.50, 0.75)),
                      use_ratio = FALSE)$aki,
   "use_ratio=FALSE → 1.5 倍不再触发（48h 规则本就不达标）")
r <- kdigo_first_event(L("a", c(0, 30), c(1.00, 1.30)))
ok(r$aki, "基准：48h 内升 0.30 达标")
ok(!kdigo_first_event(L("a", c(0, 30), c(1.00, 1.30)),
                      use_48h = FALSE)$aki,
   "use_48h=FALSE → 0.30 规则不再触发")
## 两条都关：任何升高都不再构成事件
ok(!kdigo_first_event(L("a", c(0, 30), c(1.00, 3.00)),
                      use_48h = FALSE, use_ratio = FALSE)$aki,
   "两条都关 → 3 倍升高也不再是 AKI（确认开关真的生效）")
## 只留 48h 时，分期上限受限于该规则（1.5 倍已达 2 期，关掉后只剩 1 期）
r <- kdigo_first_event(L("a", c(0, 30), c(1.00, 2.00)), use_ratio = FALSE)
ok(r$aki && r$stage_first == 1L,
   "只留 48h 规则：1.00→2.00 由 2 期降为 1 期")

cat("\n== 16. KDIGO 尿量标准（S8）==\n")
## 分母必须是**实测覆盖时长**而非名义窗长。两条相隔 10 小时、各 700 mL 的记录，
## 在 t=10 处 6 h 窗内只有后一条，但它代表的时长是 10 h：
## 速率 = 700/(70×10) = 1.0；误用名义窗长会算成 700/(70×6) = 1.667。
s <- .uo_rate_series(c(0, 10), c(700, 700), 70)
ok(near(s$rate6[2], 1.0), "分母用实测覆盖时长：700 mL / 70 kg / 10 h = 1.0（非 1.667）")
## 窗必须两端受限。前后尿量不同时，若窗只写下界，t=10 之后并不存在记录，
## 尚看不出差别；下面这条把它造出来：前段少尿、后段多尿，t=6 处不得被后段抬高。
s <- .uo_rate_series(0:12, c(rep(5, 7), rep(500, 6)), 70)
ok(s$rate6[7] < 0.5,
   "窗右端受限：t=6 只看 (0,6]，不被其后的多尿段抬高（否则会漏判少尿）")

W <- c(a = 70)
## 覆盖不足 → 不评估。首条记录的时长自入科起算，故两条记录都在 t≈1.5 h 时，
## 6/12/24 h 三个窗的覆盖时长都只有 1.5 h，全为 NA——这正是"没有记录"与
## "记录到的是正常尿量"必须分开的地方。
r <- kdigo_uo_first(L("a", c(1, 1.5), c(500, 500)), W)
ok(!r$uo_pre && !r$uo_aki && r$n_uo_pre_eval == 0L,
   "覆盖 1.5 h < 6 h → 三窗全 NA，判为不可评估而非正常")

## 25 mL/h × 70 kg = 0.357 mL/kg/h：低于 0.5（6/12 h 期）但不低于 0.3（24 h 期），
## 故最高只到 2 期。若 24 h 窗被判成 3 期，说明窗口长度或阈值取错了。
r <- kdigo_uo_first(L("a", 0:35, rep(25, 36)), W)
ok(r$uo_aki && r$uo_stage_max == 2L,
   "0.357 mL/kg/h → 6/12 h 窗达标（2 期），24 h 窗未达标（非 3 期）")

## 无尿 12 h → 3 期，与速率阈值无关
r <- kdigo_uo_first(L("a", 0:35, rep(0, 36)), W)
ok(r$uo_stage_max == 3L, "连续 12 h 尿量 0 → 无尿，3 期")

## 界标边界：t=24 属界标前，t=25 起才进结局窗。
## 0–24 h 每小时 5 mL（少尿）、之后每小时 600 mL（恢复正常）。
r <- kdigo_uo_first(L("a", c(0:24, 25:60), c(rep(5, 25), rep(600, 36))), W)
ok(r$uo_pre && !r$uo_aki,
   "0–24 h 少尿、之后恢复 → uo_pre=TRUE 且 uo_aki=FALSE（界标切分正确）")

## 结局窗右端：(24, 72] 闭于 72。造 66 h 前正常、67 h 起少尿的序列，
## 首次达标恰在 t=72；t=73 同样达标但已在窗外，不得被取作事件时刻。
r <- kdigo_uo_first(L("a", 1:73, c(rep(300, 66), rep(5, 7))), W)
ok(r$uo_aki && near(r$uo_t_event, 72),
   "事件恰落在 t=72（闭区间内），t=73 的达标不计入")

## 体重缺失不能崩、也不能默认成正常
r <- kdigo_uo_first(L("a", 0:35, rep(25, 36)), c(a = NA_real_))
ok(!r$uo_aki && r$n_uo_post_eval == 0L, "体重缺失 → 不评估（不是判为无 AKI）")

## 体重按 id 查表，不是按行序
r <- kdigo_uo_first(L("b", 0:35, rep(25, 36)), c(a = 70, b = 70))
ok(r$uo_aki, "多 id 时体重按 id 取到（b 的体重来自查表）")

cat("\n== 17. 队列组装：单用肌酐与复合口径两条路径都要能跑通 ==\n")
## 回归测试。曾经的写法是 `d$comp & d$n_uo_pre_eval > 0L`：单用肌酐时（uo = NULL）
## 那两列根本不存在，`NULL > 0L` 得 logical(0)，而 R 的 `logical(n) & logical(0)`
## **不报错**，直接塌成 logical(0)，一路带到下面的 `$<-` 才炸成
## "replacement has 0 rows, data has 8919"。错误点离病因隔了好几行，
## 所以这里直接把两种模式都跑一遍钉住。
source(file.path(.here, "cohort.R"))

## 两个模拟患者：s1 肌酐全程平稳（无 AKI），s2 在 t=48 从 1.0 升到 1.6（≥0.3 → AKI）。
## 两人界标前都有 3 次肌酐、界标后都有 2 次，故都"可判定"。
mk <- function(id) {
  list(sod = L(id, c(0, 12, 24), c(138, 140, 134)),
       cr  = L(id, c(0, 12, 24, 30, 48, 72), c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0)))
}
p1 <- mk("s1"); p2 <- mk("s2")
p2$cr$value[5] <- 1.6

fake_ex <- list(
  cohort = data.frame(id = c("s1", "s2"), subject_id = c(1L, 2L),
                      intime = as.POSIXct(c("2020-01-01 00:00", "2020-01-02 00:00"), tz = "UTC"),
                      esrd_excluded = c(0L, 0L), died_within_24h = c(0L, 0L),
                      still_in_icu = c(1L, 1L), stringsAsFactors = FALSE),
  sodium     = rbind(p1$sod, p2$sod),
  creatinine = rbind(p1$cr,  p2$cr),
  rrt        = data.frame(id = character(0), t_hours = numeric(0), value = numeric(0)),
  esrd_chart = character(0),
  attrition  = list(N0_icu_stays_all = 100L, N1_adult = 100L))

r1 <- try(run_mimic_cohort(ex = fake_ex), silent = TRUE)
ok(!inherits(r1, "try-error"), "单用肌酐（uo = NULL）能跑通，不报 replacement has 0 rows")
if (!inherits(r1, "try-error")) {
  ok(length(r1$data$not_judgeable) == nrow(r1$data) && !any(is.na(r1$data$not_judgeable)),
     "not_judgeable 与人群等长且无 NA（不是塌成 logical(0)）")
  ok(!any(r1$data$comp), "单用肌酐时 comp 全为 FALSE")
  ok(nrow(r1$analysis) == 2L && sum(r1$analysis$aki) == 1L &&
       r1$analysis$aki[r1$analysis$id == "s2"],
     "两人均可判定；仅 s2 在 t=48 满足 ≥0.3 mg/dL 而判为 AKI")
}

## 复合口径：每小时 100 mL / 70 kg = 1.43 mL/kg/h，三窗都不达少尿阈值，
## 故 uo 侧判为阴性，aki 应与单用肌酐完全一致——这样才说明复合口径在
## 尿量正常时不会凭空造出事件。
u70 <- c(s1 = 70, s2 = 70)
ur  <- rbind(data.frame(id = "s1", t_hours = 0:72, value = 100),
             data.frame(id = "s2", t_hours = 0:72, value = 100))
r2 <- try(run_mimic_cohort(ex = fake_ex, uo = list(urine = ur, weight = u70)), silent = TRUE)
ok(!inherits(r2, "try-error"), "复合口径（uo 给定）能跑通")
if (!inherits(r2, "try-error") && !inherits(r1, "try-error")) {
  ok(all(r2$data$comp), "复合口径下 comp 全为 TRUE")
  ok(identical(r2$analysis$aki, r1$analysis$aki) &&
       identical(r2$analysis$id,  r1$analysis$id),
     "尿量正常时复合口径与单用肌酐给出同一批事件（不凭空造事件）")
  ok(!any(r2$data$uo_pre), "1.43 mL/kg/h 不构成界标前少尿")
}

cat("\n== 18. 公式环境：非 data= 参数（如 weights）要能从调用帧取到 ==\n")
## 回归测试。`.fml` 原先用 as.formula 的默认环境，公式的环境是 .fml 自己的帧，
## 于是 glm(..., data = dd, weights = dd$..w) 里 glm 会按
## eval(extras, data, env) 先找模型框、再找 environment(formula)，两处都没有 dd，
## 报出 "object 'dd' not found"——错误信息里完全看不出是公式环境的问题。
source(file.path(.here, "models.R"))

set.seed(11)
n <- 200L
dd <- data.frame(
  aki  = rbinom(n, 1, 0.3),
  x    = rnorm(n),
  z    = rnorm(n))
dd$..w <- runif(n, 0.5, 2)

f <- .fml("aki", c("x", "z"))
ok(identical(environment(f), environment()),
   ".fml 生成的公式带调用帧环境，不再是 .fml 自己的帧")

m <- try(glm(f, data = dd, family = binomial(), weights = dd$..w), silent = TRUE)
ok(!inherits(m, "try-error") && isTRUE(m$prior.weights[1] > 0),
   "weights 写成 dd$..w 时能拟合（回归：object 'dd' not found）")

m2 <- glm(aki ~ x + z, data = dd, family = binomial(), weights = dd$..w)
ok(near(coef(m)[["x"]], coef(m2)[["x"]], 1e-9),
   "两种写法的系数一致，说明环境回退没有改变拟合结果")

## 固定效应仍应优先取 data= 里的列，而不是同名外部对象：
## 否则环境一变，模型用的是哪份数据就说不清了。
x <- rep(999, n)
m3 <- glm(f, data = dd, family = binomial())
ok(near(coef(m3)[["x"]], coef(glm(aki ~ x + z, data = dd, family = binomial()))[["x"]], 1e-9),
   "data= 里的 x 优先于调用帧中同名的 x")

cat("\n== 19. IPW 的目标人群：只差“界标后复测”这一道门 ==\n")
## 权重模型要在“除复测外全部入选条件均已满足”的人群上拟合，而不是全部表型人群。
## 四类反例各造一人，确认恰好只有该进来的那个被算进来。
source(file.path(.here, "bias.R"))

mk <- function(has_delta, aki_pre, not_judge) data.frame(
  has_delta_na = has_delta, aki_pre = aki_pre, not_judgeable = not_judge,
  stringsAsFactors = FALSE)
d19 <- rbind(
  mk(TRUE,  FALSE, FALSE),   # 唯一合格者：只差复测
  mk(FALSE, FALSE, FALSE),   # 没有 ΔNa
  mk(TRUE,  TRUE,  FALSE),   # 界标前已患 AKI
  mk(TRUE,  FALSE, TRUE))    # 界标前判不了
ok(identical(.eligible_except_retest(d19), c(TRUE, FALSE, FALSE, FALSE)),
   "只保留“除复测外条件全满足”的那一个")

## NA 不应被当成 TRUE 混进来（`%in% TRUE` 与裸判断的区别就在这里）
d19b <- data.frame(has_delta_na = c(NA, TRUE), aki_pre = c(FALSE, NA),
                   not_judgeable = c(FALSE, FALSE))
ok(identical(.eligible_except_retest(d19b), c(FALSE, FALSE)),
   "NA 一律判为不纳入，不因 NA 而误入目标人群")

cat(sprintf("\n==== 通过 %d，失败 %d ====\n", .pass, .fail))
if (.fail > 0L) quit(status = 1L)
