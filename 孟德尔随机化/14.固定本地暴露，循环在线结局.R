suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
  library(ggplot2)
})

# 1、读取结局ID文件

# 结局ID文件路径，第一列为OpenGWAS结局ID，不包含表头
outcome_id_file <- "结局ID.csv"

if (!file.exists(outcome_id_file)) {
  stop("没有找到结局ID文件：", outcome_id_file)
}

outcome_id_data <- read.csv(
  outcome_id_file,
  header = FALSE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

outcomes <- trimws(as.character(outcome_id_data[[1]]))
outcomes <- outcomes[!is.na(outcomes) & outcomes != ""]

if (length(outcomes) == 0) {
  stop("结局ID文件第一列中没有有效的结局ID。")
}


# 2、读取过滤后的暴露SNP数据

# 已经格式化并完成筛选的暴露数据CSV文件
exposure_data_file <- "过滤后的暴露SNPs.csv"

if (!file.exists(exposure_data_file)) {
  stop("没有找到暴露SNP数据文件：", exposure_data_file)
}

exposure_dat <- as.data.frame(
  data.table::fread(exposure_data_file)
)

required_exposure_cols <- c(
  "SNP",
  "beta.exposure",
  "se.exposure",
  "effect_allele.exposure",
  "other_allele.exposure",
  "eaf.exposure",
  "exposure",
  "id.exposure"
)

missing_exposure_cols <- setdiff(
  required_exposure_cols,
  colnames(exposure_dat)
)

if (length(missing_exposure_cols) > 0) {
  stop(
    "暴露数据缺少以下必要列：",
    paste(missing_exposure_cols, collapse = "、")
  )
}

if (nrow(exposure_dat) == 0) {
  stop("暴露数据中没有可用于分析的SNP。")
}


# 3、设置暴露和结局数据协调参数

# 处理回文SNP和等位基因方向的严格程度，可设置为1、2或3
action_value <- 2

if (!action_value %in% 1:3) {
  stop("action_value只能设置为1、2或3。")
}


# 4、设置结果保存参数

# 每个结局详细结果的保存文件夹
output_folder <- "原始数据"

# 汇总结果文件
final_result_file <- "final_results(outcome).csv"

# 错误日志文件
error_log_file <- "error_log.csv"

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}


# 5、初始化汇总结果

final_results <- data.frame(
  id.exposure = character(),
  nSNP = numeric(),
  P1 = numeric(),
  P2 = numeric(),
  P3 = numeric(),
  P4 = numeric(),
  P5 = numeric(),
  b1 = numeric(),
  b2 = numeric(),
  b3 = numeric(),
  b4 = numeric(),
  b5 = numeric(),
  id.outcome = character(),
  heterogeneity1 = numeric(),
  heterogeneity2 = numeric(),
  pleiotropy = numeric(),
  or = numeric(),
  or_lci95 = numeric(),
  or_uci95 = numeric(),
  outcome_name = character(),
  stringsAsFactors = FALSE
)


# 6、初始化错误日志

error_log <- data.frame(
  outcome_id = character(),
  error_message = character(),
  stringsAsFactors = FALSE
)


# 7、循环处理每个结局ID

for (i in seq_along(outcomes)) {
  outcome_id <- outcomes[i]
  
  message(
    "[", i, "/", length(outcomes),
    "] 正在处理结局ID：", outcome_id
  )
  
  processing_error <- tryCatch({
    
    # 8、创建当前结局的结果文件夹
    
    folder_name <- file.path(output_folder, outcome_id)
    
    if (!dir.exists(folder_name)) {
      dir.create(folder_name, recursive = TRUE)
    }
    
    
    # 9、在线提取当前结局GWAS数据
    
    outcome_dat <- TwoSampleMR::extract_outcome_data(
      snps = exposure_dat$SNP,
      outcomes = outcome_id
    )
    
    if (is.null(outcome_dat) || nrow(outcome_dat) == 0) {
      stop("没有提取到当前结局对应的GWAS数据。")
    }
    
    
    # 10、协调暴露数据和结局数据
    
    dat <- TwoSampleMR::harmonise_data(
      exposure_dat = exposure_dat,
      outcome_dat = outcome_dat,
      action = action_value
    )
    
    write.csv(
      dat,
      file = file.path(folder_name, "dat_results.csv"),
      row.names = FALSE
    )
    
    if (is.null(dat) || nrow(dat) == 0) {
      stop("暴露和结局数据协调后没有剩余SNP。")
    }
    
    if (
      "mr_keep" %in% colnames(dat) &&
      !any(dat$mr_keep, na.rm = TRUE)
    ) {
      stop("数据协调后没有可用于MR分析的SNP。")
    }
    
    
    # 11、执行孟德尔随机化分析
    
    res <- TwoSampleMR::mr(dat)
    
    if (is.null(res) || nrow(res) == 0) {
      stop("孟德尔随机化分析没有返回结果。")
    }
    
    res1 <- TwoSampleMR::generate_odds_ratios(res)
    
    write.csv(
      res1,
      file = file.path(folder_name, "mr-result.csv"),
      row.names = FALSE
    )
    
    print(res)
    
    
    # 12、执行异质性检验
    
    het <- TwoSampleMR::mr_heterogeneity(dat)
    
    write.csv(
      het,
      file = file.path(folder_name, "heterogeneity.csv"),
      row.names = FALSE
    )
    
    
    # 13、执行水平多效性检验
    
    pleio <- TwoSampleMR::mr_pleiotropy_test(dat)
    
    write.csv(
      pleio,
      file = file.path(folder_name, "pleiotropy_test.csv"),
      row.names = FALSE
    )
    
    
    # 14、提取五种MR分析方法的结果
    
    p_values <- rep(NA_real_, 5)
    b_values <- rep(NA_real_, 5)
    
    result_number <- min(5, nrow(res))
    
    if (result_number > 0) {
      p_values[seq_len(result_number)] <- res$pval[seq_len(result_number)]
      b_values[seq_len(result_number)] <- res$b[seq_len(result_number)]
    }
    
    
    # 15、提取异质性检验结果
    
    heterogeneity_values <- rep(NA_real_, 2)
    
    if (!is.null(het) && nrow(het) > 0) {
      heterogeneity_number <- min(2, nrow(het))
      
      heterogeneity_values[seq_len(heterogeneity_number)] <-
        het$Q_pval[seq_len(heterogeneity_number)]
    }
    
    
    # 16、提取水平多效性检验结果
    
    pleiotropy_value <- NA_real_
    
    if (!is.null(pleio) && nrow(pleio) > 0) {
      pleiotropy_value <- pleio$pval[1]
    }
    
    
    # 17、提取逆方差加权法的OR结果
    
    ivw_row <- which(
      res1$method == "Inverse variance weighted"
    )
    
    ivw_or <- NA_real_
    ivw_or_lci95 <- NA_real_
    ivw_or_uci95 <- NA_real_
    
    if (length(ivw_row) > 0) {
      ivw_row <- ivw_row[1]
      ivw_or <- res1$or[ivw_row]
      ivw_or_lci95 <- res1$or_lci95[ivw_row]
      ivw_or_uci95 <- res1$or_uci95[ivw_row]
    }
    
    
    # 18、汇总当前结局的主要分析结果
    
    new_row <- data.frame(
      id.exposure = res$id.exposure[1],
      nSNP = res$nsnp[1],
      P1 = p_values[1],
      P2 = p_values[2],
      P3 = p_values[3],
      P4 = p_values[4],
      P5 = p_values[5],
      b1 = b_values[1],
      b2 = b_values[2],
      b3 = b_values[3],
      b4 = b_values[4],
      b5 = b_values[5],
      id.outcome = res$id.outcome[1],
      heterogeneity1 = heterogeneity_values[1],
      heterogeneity2 = heterogeneity_values[2],
      pleiotropy = pleiotropy_value,
      or = ivw_or,
      or_lci95 = ivw_or_lci95,
      or_uci95 = ivw_or_uci95,
      outcome_name = res1$outcome[1],
      stringsAsFactors = FALSE
    )
    
    final_results <- rbind(final_results, new_row)
    
    write.csv(
      final_results,
      file = final_result_file,
      row.names = FALSE
    )
    
    
    # 19、执行逐个剔除检验和单SNP分析
    
    res_loo <- TwoSampleMR::mr_leaveoneout(dat)
    res_single <- TwoSampleMR::mr_singlesnp(dat)
    
    
    # 20、生成敏感性分析图
    
    P1_plot <- suppressWarnings(
      TwoSampleMR::mr_leaveoneout_plot(res_loo)
    )
    
    P2_plot <- TwoSampleMR::mr_scatter_plot(
      mr_results = res,
      dat = dat
    )
    
    P3_plot <- suppressWarnings(
      TwoSampleMR::mr_forest_plot(res_single)
    )
    
    P4_plot <- TwoSampleMR::mr_funnel_plot(res_single)
    
    
    # 21、保存逐个剔除检验图
    
    if (length(P1_plot) > 0 && !is.null(P1_plot[[1]])) {
      ggplot2::ggsave(
        filename = file.path(folder_name, "剔除检验.pdf"),
        plot = P1_plot[[1]],
        width = 8,
        height = 7
      )
    }
    
    
    # 22、保存效应散点图
    
    if (length(P2_plot) > 0 && !is.null(P2_plot[[1]])) {
      ggplot2::ggsave(
        filename = file.path(folder_name, "效应散点图.pdf"),
        plot = P2_plot[[1]],
        width = 8,
        height = 7
      )
    }
    
    
    # 23、保存效应森林图
    
    if (length(P3_plot) > 0 && !is.null(P3_plot[[1]])) {
      ggplot2::ggsave(
        filename = file.path(folder_name, "效应森林图.pdf"),
        plot = P3_plot[[1]],
        width = 8,
        height = 7
      )
    }
    
    
    # 24、保存漏斗图
    
    if (length(P4_plot) > 0 && !is.null(P4_plot[[1]])) {
      ggplot2::ggsave(
        filename = file.path(folder_name, "漏斗图.pdf"),
        plot = P4_plot[[1]],
        width = 8,
        height = 7
      )
    }
    
    message("结局ID处理完成：", outcome_id)
    
    NA_character_
    
  }, error = function(e) {
    conditionMessage(e)
  })
  
  
  # 25、记录当前结局的错误信息
  
  if (!is.na(processing_error)) {
    error_info <- data.frame(
      outcome_id = outcome_id,
      error_message = processing_error,
      stringsAsFactors = FALSE
    )
    
    error_log <- rbind(error_log, error_info)
    
    message(
      "处理结局ID失败：", outcome_id,
      "；错误信息：", processing_error
    )
    
    write.csv(
      error_log,
      file = error_log_file,
      row.names = FALSE
    )
  }
}


# 26、保存最终汇总结果和错误日志

write.csv(
  final_results,
  file = final_result_file,
  row.names = FALSE
)

write.csv(
  error_log,
  file = error_log_file,
  row.names = FALSE
)


# 27、输出结果说明

message(
  "分析完成。\n",
  "final_results(outcome).csv中的P1-P5和b1-b5依次对应：\n",
  "1、MR Egger回归；只有一个SNP时为Wald ratio\n",
  "2、加权中位数法\n",
  "3、逆方差加权法\n",
  "4、简单模式法\n",
  "5、加权模式法\n",
  "主要结果通常以逆方差加权法为准。"
)