suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
  library(ggplot2)
})

# 1、设置批量分析参数

# Clumping结果文件夹路径
clump_folder <- "clump_file"

# 每个暴露分析结果的保存文件夹
output_folder <- "原始数据"

# 在线结局GWAS的ID
outcome_id <- "ieu-a-1165"

# 暴露与结局数据合并时的严格程度，可设置为1、2或3
action_value <- 2

if (!dir.exists(clump_folder)) {
  stop("Clump文件夹不存在：", clump_folder)
}

if (!action_value %in% 1:3) {
  stop("action_value只能设置为1、2或3。")
}

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}


# 2、初始化汇总结果和错误日志

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
  exposure_name = character(),
  stringsAsFactors = FALSE
)

error_log <- data.frame(
  exposure_id = character(),
  error_message = character(),
  stringsAsFactors = FALSE
)


# 3、获取Clump文件夹中的CSV和TXT文件

file_paths <- list.files(
  path = clump_folder,
  pattern = "\\.csv$|\\.txt$",
  full.names = TRUE
)

if (length(file_paths) == 0) {
  stop("Clump文件夹中没有找到CSV或TXT文件：", clump_folder)
}


# 4、依次处理每个暴露数据文件

for (i in seq_along(file_paths)) {
  file_path <- file_paths[i]
  exposure_id <- basename(file_path)
  
  message("[", i, "/", length(file_paths), "] 正在处理：", exposure_id)
  
  processing_error <- tryCatch({
    
    # 5、创建当前暴露的结果保存文件夹
    
    folder_name <- file.path(output_folder, exposure_id)
    
    if (!dir.exists(folder_name)) {
      dir.create(folder_name, recursive = TRUE)
    }
    
    
    # 6、读取暴露数据
    
    if (grepl("\\.csv$", file_path)) {
      exposure_dat <- data.table::fread(file_path)
    } else if (grepl("\\.txt$", file_path)) {
      exposure_dat <- TwoSampleMR::read_exposure_data(
        filename = file_path,
        sep = "\t",
        snp_col = "SNP",
        beta_col = "beta",
        se_col = "se",
        effect_allele_col = "A1",
        other_allele_col = "A2",
        eaf_col = "eaf",
        pval_col = "p",
        samplesize_col = "n"
      )
    }
    
    required_cols <- c(
      "SNP",
      "effect_allele.exposure",
      "other_allele.exposure",
      "beta.exposure",
      "se.exposure",
      "pval.exposure",
      "samplesize.exposure",
      "exposure",
      "eaf.exposure"
    )
    
    missing_cols <- setdiff(required_cols, colnames(exposure_dat))
    
    if (length(missing_cols) > 0) {
      stop("文件缺少必要的列：", paste(missing_cols, collapse = "、"))
    }
    
    if (nrow(exposure_dat) == 0) {
      stop("暴露数据中没有可用于分析的SNP。")
    }
    
    
    # 7、计算R2和F值并筛选强工具变量
    
    if (all(is.na(exposure_dat$eaf.exposure))) {
      message("eaf.exposure列全部为NA，跳过F值计算：", exposure_id)
    } else {
      exposure_dat$R2 <- 2 * (1 - exposure_dat$eaf.exposure) *
        exposure_dat$eaf.exposure *
        (exposure_dat$beta.exposure)^2
      
      exposure_dat$F <- exposure_dat$R2 / (1 - exposure_dat$R2) *
        (exposure_dat$samplesize.exposure - 2)
      
      exposure_dat <- subset(exposure_dat, F > 10)
      
      write.csv(
        exposure_dat,
        file = file.path(folder_name, "exposure_dat-SNPs-F.csv"),
        row.names = FALSE
      )
    }
    
    if (nrow(exposure_dat) == 0) {
      stop("F值大于10筛选后没有剩余SNP。")
    }
    
    
    # 8、在线提取结局数据
    
    outcome_dat <- TwoSampleMR::extract_outcome_data(
      snps = exposure_dat$SNP,
      outcomes = outcome_id
    )
    
    if (is.null(outcome_dat) || nrow(outcome_dat) == 0) {
      stop("未提取到对应的结局GWAS数据。")
    }
    
    
    # 9、协调暴露数据和结局数据
    
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
    
    if (nrow(dat) == 0) {
      stop("暴露与结局数据协调后没有剩余SNP。")
    }
    
    
    # 10、执行孟德尔随机化分析
    
    res <- TwoSampleMR::mr(dat)
    res1 <- TwoSampleMR::generate_odds_ratios(res)
    
    write.csv(
      res1,
      file = file.path(folder_name, "mr-result.csv"),
      row.names = FALSE
    )
    
    print(res)
    
    if (nrow(res) == 0) {
      stop("孟德尔随机化分析未返回结果。")
    }
    
    
    # 11、执行异质性检验
    
    het <- TwoSampleMR::mr_heterogeneity(dat)
    
    write.csv(
      het,
      file = file.path(folder_name, "heterogeneity.csv"),
      row.names = FALSE
    )
    
    
    # 12、执行水平多效性检验
    
    pleio <- TwoSampleMR::mr_pleiotropy_test(dat)
    
    write.csv(
      pleio,
      file = file.path(folder_name, "pleiotropy_test.csv"),
      row.names = FALSE
    )
    
    
    # 13、汇总当前暴露的主要分析结果
    
    new_row <- data.frame(
      id.exposure = exposure_id,
      nSNP = res$nsnp[1],
      P1 = res$pval[1],
      P2 = res$pval[2],
      P3 = res$pval[3],
      P4 = res$pval[4],
      P5 = res$pval[5],
      b1 = res$b[1],
      b2 = res$b[2],
      b3 = res$b[3],
      b4 = res$b[4],
      b5 = res$b[5],
      id.outcome = res$id.outcome[1],
      heterogeneity1 = het$Q_pval[1],
      heterogeneity2 = het$Q_pval[2],
      pleiotropy = pleio$pval[1],
      or = res1$or[3],
      or_lci95 = res1$or_lci95[3],
      or_uci95 = res1$or_uci95[3],
      exposure_name = res$exposure[1],
      stringsAsFactors = FALSE
    )
    
    final_results <- rbind(final_results, new_row)
    
    write.csv(
      final_results,
      file = "final_results(exposure).csv",
      row.names = FALSE
    )
    
    
    # 14、执行逐个剔除和单SNP分析
    
    res_loo <- TwoSampleMR::mr_leaveoneout(dat)
    res_single <- TwoSampleMR::mr_singlesnp(dat)
    
    
    # 15、生成敏感性分析图
    
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
    
    
    # 16、保存逐个剔除检验图
    
    if (length(P1_plot) > 0) {
      ggplot2::ggsave(
        filename = file.path(folder_name, "剔除检验.pdf"),
        plot = P1_plot[[1]],
        width = 8,
        height = 7
      )
    }
    
    
    # 17、保存效应散点图
    
    if (length(P2_plot) > 0) {
      ggplot2::ggsave(
        filename = file.path(folder_name, "效应散点图.pdf"),
        plot = P2_plot[[1]],
        width = 8,
        height = 7
      )
    }
    
    
    # 18、保存效应森林图
    
    if (length(P3_plot) > 0) {
      ggplot2::ggsave(
        filename = file.path(folder_name, "效应森林图.pdf"),
        plot = P3_plot[[1]],
        width = 8,
        height = 7
      )
    }
    
    
    # 19、保存漏斗图
    
    if (length(P4_plot) > 0) {
      ggplot2::ggsave(
        filename = file.path(folder_name, "漏斗图.pdf"),
        plot = P4_plot[[1]],
        width = 8,
        height = 7
      )
    }
    
    NA_character_
    
  }, error = function(e) {
    conditionMessage(e)
  })
  
  
  # 20、记录当前文件的错误信息
  
  if (!is.na(processing_error)) {
    error_info <- data.frame(
      exposure_id = exposure_id,
      error_message = processing_error,
      stringsAsFactors = FALSE
    )
    
    error_log <- rbind(error_log, error_info)
    
    message(
      "处理暴露数据失败：", exposure_id,
      "；错误信息：", processing_error
    )
    
    write.csv(
      error_log,
      file = "error_log.csv",
      row.names = FALSE
    )
  }
}


# 21、保存最终汇总结果和错误日志

write.csv(
  final_results,
  file = "final_results(exposure).csv",
  row.names = FALSE
)

write.csv(
  error_log,
  file = "error_log.csv",
  row.names = FALSE
)


# 22、输出结果说明

message(
  "final_results(exposure).csv中的P1-P5和b1-b5依次对应：\n",
  "1、MR Egger回归；只有一个SNP时为Wald ratio\n",
  "2、加权中位数法\n",
  "3、逆方差加权法\n",
  "4、简单模式法\n",
  "5、加权模式法\n",
  "主要结果通常以逆方差加权法为准。"
)