# 1. 加载R包

suppressPackageStartupMessages({
  library(MRPRESSO)
})


# 2. 读取孟德尔随机化汇总结果

# 输入循环寻找暴露或循环寻找结局生成的汇总结果
summary_file <- "final_results(exposure).csv"

results1 <- read.csv(summary_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

if (!any(c("id.exposure", "id.outcome") %in% colnames(results1))) {
  stop("汇总结果中至少需要包含id.exposure或id.outcome列。")
}


# 3. 读取每个暴露或结局文件夹中的协调后数据

# 存放各个dat_results.csv文件的主目录
data_dir <- "原始数据"

if (!dir.exists(data_dir)) {
  stop("指定的数据目录不存在：", data_dir)
}

sub_dirs <- list.dirs(data_dir, recursive = FALSE, full.names = TRUE)
data_list <- list()

for (sub_dir in sub_dirs) {
  file_path <- file.path(sub_dir, "dat_results.csv")
  
  if (file.exists(file_path)) {
    var_name <- basename(sub_dir)
    data_list[[var_name]] <- read.csv(file_path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  }
}

if (length(data_list) == 0) {
  stop("在数据目录的子文件夹中没有找到dat_results.csv文件。")
}

message("共读取到 ", length(data_list), " 个dat_results.csv文件。")


# 4. 逐个进行MR-PRESSO分析

# MR-PRESSO使用的列名
beta_outcome_col <- "beta.outcome"
beta_exposure_col <- "beta.exposure"
se_outcome_col <- "se.outcome"
se_exposure_col <- "se.exposure"

# MR-PRESSO分析参数
OUTLIERtest <- TRUE
DISTORTIONtest <- TRUE
NbDistribution <- 1000
SignifThreshold <- 0.05
seed <- 1234

required_cols <- c(beta_outcome_col, beta_exposure_col, se_outcome_col, se_exposure_col)
results_list <- list()

for (var_name in names(data_list)) {
  cat("\n", var_name, " 开始进行MR-PRESSO分析\n", sep = "")
  
  data <- data_list[[var_name]]
  missing_cols <- setdiff(required_cols, colnames(data))
  
  if (length(missing_cols) > 0) {
    message(var_name, " 缺少以下必要列，已跳过：", paste(missing_cols, collapse = "、"))
    next
  }
  
  result <- tryCatch(
    suppressWarnings(
      mr_presso(
        BetaOutcome = beta_outcome_col,
        BetaExposure = beta_exposure_col,
        SdOutcome = se_outcome_col,
        SdExposure = se_exposure_col,
        OUTLIERtest = OUTLIERtest,
        DISTORTIONtest = DISTORTIONtest,
        data = data,
        NbDistribution = NbDistribution,
        SignifThreshold = SignifThreshold,
        seed = seed
      )
    ),
    error = function(e) {
      message(var_name, " 的MR-PRESSO分析失败，已跳过。错误信息：", conditionMessage(e))
      NULL
    }
  )
  
  if (!is.null(result)) {
    results_list[[var_name]] <- result
    cat(var_name, " MR-PRESSO分析完成\n", sep = "")
  }
}

message("\nMR-PRESSO分析全部结束。")
message("成功完成分析的数据集数量：", length(results_list))


# 5. 将Global Test的P值添加到汇总结果

results1$MR.PRESSO <- NA

for (name in names(results_list)) {
  mr_presso_results <- results_list[[name]]$`MR-PRESSO results`
  
  if (is.null(mr_presso_results) || is.null(mr_presso_results$`Global Test`)) {
    message(name, " 没有Global Test结果，无法回填P值。")
    next
  }
  
  global_test <- mr_presso_results$`Global Test`
  
  if (!"Pvalue" %in% colnames(global_test) || nrow(global_test) == 0) {
    message(name, " 的Global Test中没有Pvalue。")
    next
  }
  
  p_value <- global_test$Pvalue[1]
  matching_rows <- rep(FALSE, nrow(results1))
  
  if ("id.exposure" %in% colnames(results1)) {
    matching_rows <- matching_rows | as.character(results1$id.exposure) == name
  }
  
  if ("id.outcome" %in% colnames(results1)) {
    matching_rows <- matching_rows | as.character(results1$id.outcome) == name
  }
  
  if (any(matching_rows, na.rm = TRUE)) {
    results1$MR.PRESSO[matching_rows] <- p_value
  } else {
    message(name, " 在汇总结果的id.exposure和id.outcome中均未找到匹配项。")
  }
}

# 只保留成功获得MR-PRESSO Global Test结果的行
results1 <- results1[!is.na(results1$MR.PRESSO), , drop = FALSE]


# 6. 判断结果属于循环暴露还是循环结局

exposure_match <- FALSE
outcome_match <- FALSE

if ("id.exposure" %in% colnames(results1)) {
  exposure_match <- any(as.character(results1$id.exposure) %in% names(results_list))
}

if ("id.outcome" %in% colnames(results1)) {
  outcome_match <- any(as.character(results1$id.outcome) %in% names(results_list))
}

if (exposure_match) {
  file_suffix <- "exposure"
} else if (outcome_match) {
  file_suffix <- "outcome"
} else {
  file_suffix <- "default"
}

summary_output_file <- paste0("MRPRESSO_result_", file_suffix, ".csv")
write.csv(results1, summary_output_file, row.names = FALSE)


# 7. 整理MR-PRESSO原始分析结果

MRPRESSO_Original <- data.frame(
  Name = character(),
  Exposure = character(),
  MR_Analysis = character(),
  Causal_Estimate = numeric(),
  Sd = numeric(),
  Tstat = numeric(),
  Pvalue = numeric(),
  RSSobs = numeric(),
  GlobalTestPvalue = numeric(),
  stringsAsFactors = FALSE
)

for (name in names(results_list)) {
  main_mr_results <- results_list[[name]]$`Main MR results`
  mr_presso_results <- results_list[[name]]$`MR-PRESSO results`
  
  if (is.null(main_mr_results) || nrow(main_mr_results) == 0) {
    message(name, " 没有Main MR results，已跳过原始结果整理。")
    next
  }
  
  if (is.null(mr_presso_results) || is.null(mr_presso_results$`Global Test`)) {
    message(name, " 没有Global Test，已跳过原始结果整理。")
    next
  }
  
  global_test <- mr_presso_results$`Global Test`
  
  if (!all(c("RSSobs", "Pvalue") %in% colnames(global_test)) || nrow(global_test) == 0) {
    message(name, " 的Global Test缺少RSSobs或Pvalue，已跳过原始结果整理。")
    next
  }
  
  rss_obs <- global_test$RSSobs[1]
  global_test_pvalue <- global_test$Pvalue[1]
  
  for (i in seq_len(nrow(main_mr_results))) {
    temp_df <- data.frame(
      Name = name,
      Exposure = main_mr_results$Exposure[i],
      MR_Analysis = main_mr_results$`MR Analysis`[i],
      Causal_Estimate = main_mr_results$`Causal Estimate`[i],
      Sd = main_mr_results$Sd[i],
      Tstat = main_mr_results$`T-stat`[i],
      Pvalue = main_mr_results$`P-value`[i],
      RSSobs = rss_obs,
      GlobalTestPvalue = global_test_pvalue,
      stringsAsFactors = FALSE
    )
    
    MRPRESSO_Original <- rbind(MRPRESSO_Original, temp_df)
  }
}

write.csv(MRPRESSO_Original, "MRPRESSO_Original.csv", row.names = FALSE)


# 8. 保留MR-PRESSO完整结果列表

MR.PRESSO_results <- results_list

message("\n结果保存完成：")
message("1. 汇总结果：", summary_output_file)
message("2. MR-PRESSO原始结果：MRPRESSO_Original.csv")