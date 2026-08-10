# 1. 读取暴露到结局的MR结果

# 暴露到结局结果文件
exposure_outcome_file <- "暴露到结局.csv"

Medium1 <- read.csv(
  exposure_outcome_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


# 2. 读取暴露到中介的MR结果

# 暴露到中介结果文件
exposure_mediator_file <- "暴露到中介.csv"

Medium2 <- read.csv(
  exposure_mediator_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


# 3. 读取中介到结局的MR结果

# 中介到结局结果文件
mediator_outcome_file <- "中介到结局.csv"

Medium3 <- read.csv(
  mediator_outcome_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


# 4. 检查三个输入文件

if (nrow(Medium1) == 0 || nrow(Medium2) == 0 || nrow(Medium3) == 0) {
  stop("至少有一个输入文件为空，请检查输入文件。")
}

required_cols1 <- c("id.exposure", "b", "method")
required_cols2 <- c("id.outcome", "b", "se")
required_cols3 <- c("id.outcome", "b", "se")

missing_cols1 <- setdiff(required_cols1, colnames(Medium1))
missing_cols2 <- setdiff(required_cols2, colnames(Medium2))
missing_cols3 <- setdiff(required_cols3, colnames(Medium3))

if (length(missing_cols1) > 0) {
  stop("暴露到结局文件缺少以下必要列：", paste(missing_cols1, collapse = "、"))
}

if (length(missing_cols2) > 0) {
  stop("暴露到中介文件缺少以下必要列：", paste(missing_cols2, collapse = "、"))
}

if (length(missing_cols3) > 0) {
  stop("中介到结局文件缺少以下必要列：", paste(missing_cols3, collapse = "、"))
}

if (length(unique(c(nrow(Medium1), nrow(Medium2), nrow(Medium3)))) != 1) {
  stop(
    "三个输入文件的行数不一致：暴露到结局为", nrow(Medium1),
    "行，暴露到中介为", nrow(Medium2),
    "行，中介到结局为", nrow(Medium3), "行。"
  )
}


# 5. 检查三个文件中的MR方法顺序

if ("method" %in% colnames(Medium2)) {
  method1 <- trimws(as.character(Medium1$method))
  method2 <- trimws(as.character(Medium2$method))
  
  if (!identical(method1, method2)) {
    stop("暴露到结局与暴露到中介文件的method列或方法顺序不一致。")
  }
}

if ("method" %in% colnames(Medium3)) {
  method1 <- trimws(as.character(Medium1$method))
  method3 <- trimws(as.character(Medium3$method))
  
  if (!identical(method1, method3)) {
    stop("暴露到结局与中介到结局文件的method列或方法顺序不一致。")
  }
}


# 6. 将计算所需列转换为数值型

Medium1$b <- suppressWarnings(as.numeric(Medium1$b))
Medium2$b <- suppressWarnings(as.numeric(Medium2$b))
Medium2$se <- suppressWarnings(as.numeric(Medium2$se))
Medium3$b <- suppressWarnings(as.numeric(Medium3$b))
Medium3$se <- suppressWarnings(as.numeric(Medium3$se))

if (anyNA(Medium1$b)) {
  stop("暴露到结局文件的b列存在缺失值或无法转换为数值的内容。")
}

if (anyNA(Medium2$b) || anyNA(Medium2$se)) {
  stop("暴露到中介文件的b列或se列存在缺失值或无法转换为数值的内容。")
}

if (anyNA(Medium3$b) || anyNA(Medium3$se)) {
  stop("中介到结局文件的b列或se列存在缺失值或无法转换为数值的内容。")
}

if (any(Medium2$se < 0) || any(Medium3$se < 0)) {
  stop("标准误se不能小于0，请检查输入数据。")
}


# 7. 计算中介效应

# 中介效应 = 暴露到中介的效应 × 中介到结局的效应
MIbeta <- Medium2$b * Medium3$b


# 8. 计算中介效应占总效应的比例

# 中介占比 = 中介效应 ÷ 暴露到结局的总效应
Medium_pro <- MIbeta / Medium1$b

zero_total_rows <- Medium1$b == 0

if (any(zero_total_rows)) {
  Medium_pro[zero_total_rows] <- NA_real_
  warning("部分暴露到结局的总效应为0，对应行的中介效应占比已设置为NA。")
}


# 9. 使用Delta法计算中介效应的标准误

S <- sqrt(
  Medium2$b^2 * Medium3$se^2 +
    Medium3$b^2 * Medium2$se^2
)


# 10. 计算Z值和P值

Z <- MIbeta / S

zero_se_rows <- S == 0

if (any(zero_se_rows)) {
  Z[zero_se_rows] <- NA_real_
  warning("部分中介效应标准误为0，对应行的Z值和P值已设置为NA。")
}

P <- 2 * pnorm(abs(Z), lower.tail = FALSE)


# 11. 整理中介孟德尔随机化结果

results <- data.frame(
  exposure = as.character(Medium1$id.exposure),
  Medium = as.character(Medium2$id.outcome),
  outcome = as.character(Medium3$id.outcome),
  method = as.character(Medium1$method),
  MIbeta = MIbeta,
  Total = Medium1$b,
  Medium_pro = Medium_pro,
  S = S,
  Z = Z,
  P = P,
  stringsAsFactors = FALSE
)


# 12. 创建结果文件夹

output_folder <- "中介MR结果"

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}


# 13. 设置结果文件名

if (
  is.na(results$exposure[1]) || results$exposure[1] == "" ||
  is.na(results$Medium[1]) || results$Medium[1] == "" ||
  is.na(results$outcome[1]) || results$outcome[1] == ""
) {
  stop("第一行的暴露ID、中介ID或结局ID缺失，无法生成结果文件名。")
}

file_name <- paste(
  results$exposure[1],
  results$Medium[1],
  results$outcome[1],
  sep = "___"
)

file_path <- file.path(output_folder, paste0(file_name, ".csv"))


# 14. 保存中介孟德尔随机化结果

write.csv(results, file = file_path, row.names = FALSE)

message("中介孟德尔随机化分析完成，结果已保存至：", file_path)