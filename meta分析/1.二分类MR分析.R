suppressPackageStartupMessages({
  library(meta)
})

# 1、读取并检查二分类Meta分析数据

# 输入文件路径
meta_file <- "meta分析二分类示例数据.csv"
required_cols <- c("Study", "p_i", "n_i", "p_c", "n_c")

if (!file.exists(meta_file)) stop("没有找到输入数据文件：", meta_file)

meta_data <- read.csv(meta_file, check.names = FALSE, stringsAsFactors = FALSE)
missing_cols <- setdiff(required_cols, colnames(meta_data))

if (length(missing_cols) > 0) {
  stop("输入数据缺少以下列：", paste(missing_cols, collapse = "、"))
}

meta_data <- meta_data[, required_cols]
meta_data$Study <- as.character(meta_data$Study)
meta_data$p_i <- suppressWarnings(as.numeric(meta_data$p_i))
meta_data$n_i <- suppressWarnings(as.numeric(meta_data$n_i))
meta_data$p_c <- suppressWarnings(as.numeric(meta_data$p_c))
meta_data$n_c <- suppressWarnings(as.numeric(meta_data$n_c))

if (any(is.na(meta_data$Study)) || any(trimws(meta_data$Study) == "")) {
  stop("Study列存在缺失值或空值。")
}

if (any(is.na(meta_data$p_i))) stop("p_i列存在非数值或缺失值。")
if (any(is.na(meta_data$n_i))) stop("n_i列存在非数值或缺失值。")
if (any(is.na(meta_data$p_c))) stop("p_c列存在非数值或缺失值。")
if (any(is.na(meta_data$n_c))) stop("n_c列存在非数值或缺失值。")

if (any(meta_data[, c("p_i", "n_i", "p_c", "n_c")] < 0)) {
  stop("事件数和非事件数不能小于0。")
}

meta_data$total_i <- meta_data$p_i + meta_data$n_i
meta_data$total_c <- meta_data$p_c + meta_data$n_c

if (any(meta_data$total_i <= 0)) stop("干预组总人数必须大于0。")
if (any(meta_data$total_c <= 0)) stop("对照组总人数必须大于0。")

if (nrow(meta_data) < 2) {
  stop("至少需要2项研究才能进行Meta分析和逐一剔除敏感性分析。")
}

print(head(meta_data, 10))


# 2、创建结果保存文件夹

# 结果保存文件夹
out_dir <- "二分类Meta分析"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}


# 3、设置参数并进行二分类Meta分析

# 效应量类型：RR、OR或RD
summary_measure <- "RR"

# 合并方法：MH、Inverse或GLMM
method_meta <- "MH"

# 是否计算固定效应模型和随机效应模型
common_model <- TRUE
random_model <- TRUE

# 随机效应置信区间方法：HK或classic
method_random_ci <- "HK"

if (!summary_measure %in% c("RR", "OR", "RD")) {
  stop("summary_measure只能设置为RR、OR或RD。")
}

if (!method_meta %in% c("MH", "Inverse", "GLMM")) {
  stop("method_meta只能设置为MH、Inverse或GLMM。")
}

if (!method_random_ci %in% c("HK", "classic")) {
  stop("method_random_ci只能设置为HK或classic。")
}

if (!common_model && !random_model) {
  stop("common_model和random_model不能同时设置为FALSE。")
}

if (method_meta == "GLMM" && summary_measure != "OR") {
  stop("GLMM方法只适用于OR，请将summary_measure设置为OR。")
}

m1 <- meta::metabin(
  event.e = p_i,
  n.e = total_i,
  event.c = p_c,
  n.c = total_c,
  studlab = Study,
  data = meta_data,
  sm = summary_measure,
  method = method_meta,
  common = common_model,
  random = random_model,
  method.random.ci = method_random_ci
)

print(m1)


# 4、进行逐一剔除敏感性分析

# 同时计算固定效应和随机效应模型时，敏感性分析优先使用随机效应模型
sensitivity_model <- if (random_model) "random" else "common"

m1_inf <- meta::metainf(
  m1,
  pooled = sensitivity_model
)

print(m1_inf)


# 5、进行偏倚检验和Trim-and-fill分析

egger_test <- tryCatch(
  meta::metabias(m1, method.bias = "Egger"),
  error = function(e) e
)

begg_test <- tryCatch(
  meta::metabias(m1, method.bias = "Begg"),
  error = function(e) e
)

tf <- tryCatch(
  meta::trimfill(m1),
  error = function(e) e
)


# 6、整理并导出主要Meta分析结果

# RR和OR存储在对数尺度，需要指数转换；RD直接使用原始效应值
if (summary_measure %in% c("RR", "OR")) {
  common_effect <- if (common_model) exp(m1$TE.common) else NA_real_
  common_lower <- if (common_model) exp(m1$lower.common) else NA_real_
  common_upper <- if (common_model) exp(m1$upper.common) else NA_real_
  
  random_effect <- if (random_model) exp(m1$TE.random) else NA_real_
  random_lower <- if (random_model) exp(m1$lower.random) else NA_real_
  random_upper <- if (random_model) exp(m1$upper.random) else NA_real_
} else {
  common_effect <- if (common_model) m1$TE.common else NA_real_
  common_lower <- if (common_model) m1$lower.common else NA_real_
  common_upper <- if (common_model) m1$upper.common else NA_real_
  
  random_effect <- if (random_model) m1$TE.random else NA_real_
  random_lower <- if (random_model) m1$lower.random else NA_real_
  random_upper <- if (random_model) m1$upper.random else NA_real_
}

main_result_table <- data.frame(
  Model = c("Common", "Random"),
  Effect = c(common_effect, random_effect),
  Lower95CI = c(common_lower, random_lower),
  Upper95CI = c(common_upper, random_upper),
  Pvalue = c(
    if (common_model) m1$pval.common else NA_real_,
    if (random_model) m1$pval.random else NA_real_
  ),
  stringsAsFactors = FALSE
)

meta_result <- data.frame(
  TE_common = common_effect,
  lower_common = common_lower,
  upper_common = common_upper,
  pval_common = if (common_model) m1$pval.common else NA_real_,
  TE_random = random_effect,
  lower_random = random_lower,
  upper_random = random_upper,
  pval_random = if (random_model) m1$pval.random else NA_real_,
  I2 = m1$I2,
  tau2 = m1$tau2,
  Q_pvalue = m1$pval.Q
)

print(main_result_table)

write.csv(
  main_result_table,
  file.path(out_dir, "main_result_table.csv"),
  row.names = FALSE
)

write.csv(
  meta_result,
  file.path(out_dir, "meta_result_summary.csv"),
  row.names = FALSE
)

write.csv(
  meta_data,
  file.path(out_dir, "meta_uploaded_data_checked.csv"),
  row.names = FALSE
)


# 7、生成并保存异质性结果解读

low_heterogeneity <- !is.na(m1$I2) &&
  !is.na(m1$pval.Q) &&
  m1$I2 < 50 &&
  m1$pval.Q > 0.10

heterogeneity_interpretation <- if (low_heterogeneity) {
  "解读：研究间异质性较低，可优先参考固定效应模型结果。"
} else {
  "解读：研究间存在一定异质性，建议优先参考随机效应模型结果。"
}

heterogeneity_text <- paste0(
  "Q检验P值：", signif(m1$pval.Q, 4), "\n",
  "I² = ", round(m1$I2, 2), "%\n",
  "Tau² = ", round(m1$tau2, 4), "\n\n",
  heterogeneity_interpretation
)

cat(heterogeneity_text, "\n")

writeLines(
  heterogeneity_text,
  file.path(out_dir, "异质性结果.txt")
)


# 8、生成并保存主要结果解读

if (random_model) {
  effect_value <- random_effect
  lower_value <- random_lower
  upper_value <- random_upper
  p_value <- m1$pval.random
  model_used <- "随机效应模型"
} else {
  effect_value <- common_effect
  lower_value <- common_lower
  upper_value <- common_upper
  p_value <- m1$pval.common
  model_used <- "固定效应模型"
}

if (summary_measure == "RR") {
  effect_text <- if (effect_value > 1) {
    "提示干预组（或暴露组）事件发生风险高于对照组。"
  } else if (effect_value < 1) {
    "提示干预组（或暴露组）事件发生风险低于对照组。"
  } else {
    "提示两组事件发生风险接近。"
  }
} else if (summary_measure == "OR") {
  effect_text <- if (effect_value > 1) {
    "提示干预组（或暴露组）事件优势比高于对照组。"
  } else if (effect_value < 1) {
    "提示干预组（或暴露组）事件优势比低于对照组。"
  } else {
    "提示两组优势比接近。"
  }
} else {
  effect_text <- if (effect_value > 0) {
    "提示干预组（或暴露组）风险差高于对照组。"
  } else if (effect_value < 0) {
    "提示干预组（或暴露组）风险差低于对照组。"
  } else {
    "提示两组风险差接近0。"
  }
}

significance_text <- if (!is.na(p_value) && p_value < 0.05) {
  "合并效应具有统计学意义。"
} else {
  "合并效应无统计学意义。"
}

heterogeneity_result_text <- if (low_heterogeneity) {
  "研究间异质性较低，结果一致性相对较好。"
} else {
  "研究间存在一定异质性，解释结果时需结合研究差异谨慎判断。"
}

interpretation_text <- paste0(
  "结果解读建议：\n",
  "1. 本次主要参考", model_used, "结果。\n",
  "2. 合并效应值 = ", round(effect_value, 4),
  "，95%CI：", round(lower_value, 4), " ~ ", round(upper_value, 4),
  "，P = ", signif(p_value, 4), "。\n",
  "3. ", effect_text, "\n",
  "4. ", significance_text, "\n",
  "5. ", heterogeneity_result_text, "\n",
  "6. 若逐一剔除任一研究后总体方向基本一致，说明结果稳健性较好；",
  "若效应方向或显著性明显改变，则提示该研究对总体结果影响较大。\n",
  "7. 漏斗图及Egger/Begg检验仅作为发表偏倚的辅助判断，",
  "不能单独作为最终结论依据。"
)

cat(interpretation_text, "\n")

writeLines(
  interpretation_text,
  file.path(out_dir, "结果解读.txt")
)


# 9、整理并保存偏倚检验结果

if (inherits(egger_test, "error")) {
  egger_text <- paste0(
    "Egger检验无法完成。原因：",
    egger_test$message
  )
} else {
  egger_text <- capture.output(print(egger_test))
}

if (inherits(begg_test, "error")) {
  begg_text <- paste0(
    "Begg检验无法完成。原因：",
    begg_test$message
  )
} else {
  begg_text <- capture.output(print(begg_test))
}

bias_test_text <- paste(
  c(
    "Egger检验结果：",
    egger_text,
    "",
    "Begg检验结果：",
    begg_text
  ),
  collapse = "\n"
)

cat(bias_test_text, "\n")

writeLines(
  bias_test_text,
  file.path(out_dir, "偏倚检验结果.txt")
)


# 10、生成并保存Trim-and-fill结果解读

if (inherits(tf, "error")) {
  trimfill_text <- paste0(
    "Trim-and-fill分析无法完成。原因：",
    tf$message
  )
} else {
  if (!is.null(tf$TE.random) && length(tf$TE.random) > 0) {
    trimfill_effect <- if (summary_measure %in% c("RR", "OR")) {
      exp(tf$TE.random)
    } else {
      tf$TE.random
    }
  } else {
    trimfill_effect <- NA_real_
  }
  
  trimfill_text <- paste0(
    "Trim-and-fill分析说明：\n",
    "1. 该分析用于作为发表偏倚的补充探索，而不是核心结论依据。\n",
    "2. 原始合并效应约为：", round(effect_value, 4), "\n",
    "3. Trim-and-fill校正后随机效应合并值约为：",
    round(trimfill_effect, 4), "\n",
    "4. 若校正前后效应变化较小，通常提示发表偏倚对总体结论影响有限；",
    "若变化较大，则需更谨慎解释。"
  )
}

cat(trimfill_text, "\n")

writeLines(
  trimfill_text,
  file.path(out_dir, "Trim-and-fill结果解读.txt")
)


# 11、设置森林图参数并保存森林图

forest_common <- common_model
forest_random <- random_model
forest_print_tau2 <- TRUE
forest_print_I2 <- TRUE
forest_print_pvalQ <- TRUE

forest_col_diamond <- "#E64B35FF"
forest_col_square <- "#4DBBD5FF"
forest_col_study <- "#00A087FF"

forest_xlab <- switch(
  summary_measure,
  RR = "Risk Ratio (RR)",
  OR = "Odds Ratio (OR)",
  RD = "Risk Difference (RD)"
)

forest_width <- 11
forest_height <- 8
forest_file_name <- "1.森林图.pdf"

pdf(
  file.path(out_dir, forest_file_name),
  width = forest_width,
  height = forest_height
)

meta::forest(
  m1,
  sortvar = m1$TE,
  common = forest_common,
  random = forest_random,
  print.tau2 = forest_print_tau2,
  print.I2 = forest_print_I2,
  print.pval.Q = forest_print_pvalQ,
  col.diamond = forest_col_diamond,
  col.square = forest_col_square,
  col.study = forest_col_study,
  xlab = forest_xlab,
  leftcols = c(
    "studlab",
    "event.e",
    "n.e",
    "event.c",
    "n.c"
  ),
  leftlabs = c(
    "Study",
    "E.event",
    "E.total",
    "C.event",
    "C.total"
  )
)

dev.off()


# 12、设置敏感性分析图参数并保存逐一剔除森林图

metainf_col_diamond <- "#E64B35FF"
metainf_col_square <- "#4DBBD5FF"

metainf_xlab <- switch(
  summary_measure,
  RR = "RR after omitting one study",
  OR = "OR after omitting one study",
  RD = "RD after omitting one study"
)

metainf_width <- 11
metainf_height <- 8
metainf_file_name <- "2.敏感性分析图.pdf"

pdf(
  file.path(out_dir, metainf_file_name),
  width = metainf_width,
  height = metainf_height
)

meta::forest(
  m1_inf,
  xlab = metainf_xlab,
  col.diamond = metainf_col_diamond,
  col.bg = metainf_col_square
)

dev.off()


# 13、设置Baujat图参数并保存Baujat图

baujat_color <- "#E64B35FF"
baujat_pch <- 19

baujat_width <- 9
baujat_height <- 7
baujat_file_name <- "3.Baujat图.pdf"

pdf(
  file.path(out_dir, baujat_file_name),
  width = baujat_width,
  height = baujat_height
)

meta::baujat(
  m1,
  col = baujat_color,
  pch = baujat_pch
)

dev.off()


# 14、设置漏斗图参数并保存漏斗图

funnel_studlab <- TRUE
funnel_color <- "#E64B35FF"
funnel_pch <- 19

funnel_xlab <- switch(
  summary_measure,
  RR = "Risk Ratio (RR)",
  OR = "Odds Ratio (OR)",
  RD = "Risk Difference (RD)"
)

funnel_ylab <- "Standard Error"

funnel_width <- 9
funnel_height <- 7
funnel_file_name <- "4.漏斗图.pdf"

pdf(
  file.path(out_dir, funnel_file_name),
  width = funnel_width,
  height = funnel_height
)

meta::funnel(
  m1,
  xlab = funnel_xlab,
  ylab = funnel_ylab,
  studlab = funnel_studlab,
  col = funnel_color,
  pch = funnel_pch
)

dev.off()


# 15、记录并保存本次分析参数

parameter_text <- paste0(
  "运行流程说明：\n",
  "1. 使用普通R脚本读取二分类结局Meta分析数据。\n",
  "2. 使用meta包的metabin()进行二分类Meta分析，",
  "合并效应量可选择RR、OR或RD。\n",
  "3. 使用metainf()进行逐一剔除敏感性分析。\n",
  "4. 使用forest()绘制森林图和逐一剔除敏感性分析图。\n",
  "5. 使用baujat()绘制Baujat图。\n",
  "6. 使用funnel()绘制漏斗图。\n",
  "7. 使用metabias()进行Egger检验与Begg检验。\n",
  "8. 使用trimfill()进行Trim-and-fill分析。\n\n",
  "本次参数总结：\n",
  "- 效应量类型（sm）：", summary_measure, "\n",
  "- 合并方法（method）：", method_meta, "\n",
  "- 固定效应模型（common）：", common_model, "\n",
  "- 随机效应模型（random）：", random_model, "\n",
  "- 随机效应置信区间方法（method.random.ci）：",
  method_random_ci, "\n",
  "- 敏感性分析使用模型：", sensitivity_model
)

writeLines(
  parameter_text,
  file.path(out_dir, "meta_binary_parameters.txt")
)

message("二分类Meta分析完成。")
message(
  "结果保存位置：",
  normalizePath(
    out_dir,
    winslash = "/",
    mustWork = FALSE
  )
)