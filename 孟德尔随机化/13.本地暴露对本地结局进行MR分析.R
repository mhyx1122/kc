suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
  library(ggplot2)
})

# 1、读取本地结局数据

# 已经格式化为TwoSampleMR结局数据的CSV文件
outcome_data_file <- "outcome_data/outcome_data_BRCA.csv"

if (!file.exists(outcome_data_file)) {
  stop("没有找到结局数据文件：", outcome_data_file)
}

outcome_dat <- as.data.frame(
  data.table::fread(outcome_data_file)
)

required_outcome_cols <- c(
  "SNP",
  "beta.outcome",
  "se.outcome",
  "effect_allele.outcome",
  "other_allele.outcome",
  "eaf.outcome",
  "outcome",
  "id.outcome"
)

missing_outcome_cols <- setdiff(
  required_outcome_cols,
  colnames(outcome_dat)
)

if (length(missing_outcome_cols) > 0) {
  stop(
    "结局数据缺少以下必要列：",
    paste(missing_outcome_cols, collapse = "、")
  )
}


# 2、设置Clump暴露数据文件夹

# 存放Clump后暴露数据的文件夹
clump_folder <- "clump_file"

if (!dir.exists(clump_folder)) {
  stop("Clump文件夹不存在：", clump_folder)
}

file_paths <- list.files(
  path = clump_folder,
  pattern = "\\.(csv|txt)$",
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(file_paths) == 0) {
  stop("Clump文件夹中没有找到CSV或TXT文件：", clump_folder)
}


# 3、设置结果保存参数

# 每个暴露分析结果的保存文件夹
output_folder <- "原始数据"

# 暴露和结局数据协调时的严格程度，可设置为1、2或3
action_value <- 2

# 汇总结果保存文件
final_result_file <- "final_results(exposure).csv"

# 错误日志保存文件
error_log_file <- "error_log.csv"

if (!action_value %in% 1:3) {
  stop("action_value只能设置为1、2或3。")
}

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}


# 4、初始化汇总结果

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


# 5、初始化错误日志

error_log <- data.frame(
  exposure_id = character(),
  error_message = character(),
  stringsAsFactors = FALSE
)


# 6、循环处理每个暴露数据文件

for (i in seq_along(file_paths)) {
  file_path <- file_paths[i]
  exposure_id <- basename(file_path)
  exposure_name <- tools::file_path_sans_ext(exposure_id)
  
  message(
    "[", i, "/", length(file_paths),
    "] 正在处理：", exposure_id
  )
  
  processing_error <- tryCatch({
    
    # 7、创建当前暴露的结果文件夹
    
    folder_name <- file.path(
      output_folder,
      exposure_id
    )
    
    if (!dir.exists(folder_name)) {
      dir.create(folder_name, recursive = TRUE)
    }
    
    
    # 8、读取暴露数据
    
    file_extension <- tolower(
      tools::file_ext(file_path)
    )
    
    if (file_extension == "csv") {
      exposure_dat <- as.data.frame(
        data.table::fread(file_path)
      )
    }
    
    if (file_extension == "txt") {
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
      
      exposure_dat$exposure <- exposure_name
      exposure_dat$id.exposure <- exposure_name
    }
    
    cat("\n现在运行到：", file_path, "\n")
    
    
    # 9、检查暴露数据必要列
    
    required_exposure_cols <- c(
      "SNP",
      "effect_allele.exposure",
      "other_allele.exposure",
      "beta.exposure",
      "se.exposure",
      "pval.exposure",
      "samplesize.exposure",
      "exposure",
      "id.exposure",
      "eaf.exposure"
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
      stop("当前暴露文件中没有可用于分析的SNP。")
    }
    
    
    # 10、计算R2和F值
    
    all_eaf_missing <- all(
      is.na(exposure_dat$eaf.exposure)
    )
    
    all_samplesize_missing <- all(
      is.na(exposure_dat$samplesize.exposure)
    )
    
    if (all_eaf_missing || all_samplesize_missing) {
      message(
        "eaf.exposure或samplesize.exposure全部为NA，",
        "跳过F值计算和筛选：", exposure_id
      )
    } else {
      exposure_dat$R2 <- 2 *
        (1 - exposure_dat$eaf.exposure) *
        exposure_dat$eaf.exposure *
        (exposure_dat$beta.exposure)^2
      
      exposure_dat$F <- exposure_dat$R2 /
        (1 - exposure_dat$R2) *
        (exposure_dat$samplesize.exposure - 2)
      
      if (any(!is.na(exposure_dat$F))) {
        exposure_dat <- subset(
          exposure_dat,
          !is.na(F) & F > 10
        )
      } else {
        message(
          "没有能够成功计算的F值，",
          "跳过F值筛选：", exposure_id
        )
      }
    }
    
    if (nrow(exposure_dat) == 0) {
      stop("F值大于10筛选后没有剩余SNP。")
    }
    
    write.csv(
      exposure_dat,
      file = file.path(
        folder_name,
        "exposure_dat-SNPs-F.csv"
      ),
      row.names = FALSE
    )
    
    
    # 11、从本地结局数据中提取共有SNP
    
    exposure_snps <- exposure_dat$SNP
    
    common_snps <- intersect(
      outcome_dat$SNP,
      exposure_snps
    )
    
    if (length(common_snps) == 0) {
      stop("暴露数据与结局数据之间没有共有SNP。")
    }
    
    outcome_datkkk <- outcome_dat[
      outcome_dat$SNP %in% common_snps,
      ,
      drop = FALSE
    ]
    
    
    # 12、协调暴露数据和结局数据
    
    dat <- TwoSampleMR::harmonise_data(
      exposure_dat = exposure_dat,
      outcome_dat = outcome_datkkk,
      action = action_value
    )
    
    write.csv(
      dat,
      file = file.path(
        folder_name,
        "dat_results.csv"
      ),
      row.names = FALSE
    )
    
    if (nrow(dat) == 0) {
      stop("数据协调后没有剩余SNP。")
    }
    
    if (
      "mr_keep" %in% colnames(dat) &&
      !any(dat$mr_keep, na.rm = TRUE)
    ) {
      stop("数据协调后没有可用于MR分析的SNP。")
    }
    
    
    # 13、执行孟德尔随机化分析
    
    res <- TwoSampleMR::mr(dat)
    
    if (is.null(res) || nrow(res) == 0) {
      stop("孟德尔随机化分析没有返回结果。")
    }
    
    res1 <- TwoSampleMR::generate_odds_ratios(res)
    
    write.csv(
      res1,
      file = file.path(
        folder_name,
        "mr-result.csv"
      ),
      row.names = FALSE
    )
    
    print(res)
    
    
    # 14、执行异质性检验
    
    het <- TwoSampleMR::mr_heterogeneity(dat)
    
    write.csv(
      het,
      file = file.path(
        folder_name,
        "heterogeneity.csv"
      ),
      row.names = FALSE
    )
    
    
    # 15、执行水平多效性检验
    
    pleio <- TwoSampleMR::mr_pleiotropy_test(dat)
    
    write.csv(
      pleio,
      file = file.path(
        folder_name,
        "pleiotropy_test.csv"
      ),
      row.names = FALSE
    )
    
    
    # 16、提取五种MR分析方法的结果
    
    p_values <- rep(
      NA_real_,
      5
    )
    
    b_values <- rep(
      NA_real_,
      5
    )
    
    result_number <- min(
      5,
      nrow(res)
    )
    
    p_values[seq_len(result_number)] <- res$pval[
      seq_len(result_number)
    ]
    
    b_values[seq_len(result_number)] <- res$b[
      seq_len(result_number)
    ]
    
    
    # 17、提取异质性和多效性结果
    
    heterogeneity_values <- rep(
      NA_real_,
      2
    )
    
    if (!is.null(het) && nrow(het) > 0) {
      heterogeneity_number <- min(
        2,
        nrow(het)
      )
      
      heterogeneity_values[
        seq_len(heterogeneity_number)
      ] <- het$Q_pval[
        seq_len(heterogeneity_number)
      ]
    }
    
    pleiotropy_value <- NA_real_
    
    if (!is.null(pleio) && nrow(pleio) > 0) {
      pleiotropy_value <- pleio$pval[1]
    }
    
    
    # 18、提取逆方差加权法OR结果
    
    ivw_or <- NA_real_
    ivw_or_lci95 <- NA_real_
    ivw_or_uci95 <- NA_real_
    
    if (nrow(res1) >= 3) {
      ivw_or <- res1$or[3]
      ivw_or_lci95 <- res1$or_lci95[3]
      ivw_or_uci95 <- res1$or_uci95[3]
    }
    
    
    # 19、汇总当前暴露的分析结果
    
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
      exposure_name = res$exposure[1],
      stringsAsFactors = FALSE
    )
    
    final_results <- rbind(
      final_results,
      new_row
    )
    
    write.csv(
      final_results,
      file = final_result_file,
      row.names = FALSE
    )
    
    
    # 20、执行逐个剔除检验和单SNP分析
    
    res_loo <- TwoSampleMR::mr_leaveoneout(dat)
    
    res_single <- TwoSampleMR::mr_singlesnp(dat)
    
    
    # 21、生成逐个剔除检验图
    
    P1_plot <- suppressWarnings(
      TwoSampleMR::mr_leaveoneout_plot(
        res_loo
      )
    )
    
    
    # 22、生成效应散点图
    
    P2_plot <- TwoSampleMR::mr_scatter_plot(
      mr_results = res,
      dat = dat
    )
    
    
    # 23、生成效应森林图
    
    P3_plot <- suppressWarnings(
      TwoSampleMR::mr_forest_plot(
        res_single
      )
    )
    
    
    # 24、生成漏斗图
    
    P4_plot <- TwoSampleMR::mr_funnel_plot(
      res_single
    )
    
    
    # 25、保存逐个剔除检验图
    
    if (length(P1_plot) > 0) {
      ggplot2::ggsave(
        filename = file.path(
          folder_name,
          "剔除检验.pdf"
        ),
        plot = P1_plot[[1]],
        width = 8,
        height = 7
      )
    }
    
    
    # 26、保存效应散点图
    
    if (length(P2_plot) > 0) {
      ggplot2::ggsave(
        filename = file.path(
          folder_name,
          "效应散点图.pdf"
        ),
        plot = P2_plot[[1]],
        width = 8,
        height = 7
      )
    }
    
    
    # 27、保存效应森林图
    
    if (length(P3_plot) > 0) {
      ggplot2::ggsave(
        filename = file.path(
          folder_name,
          "效应森林图.pdf"
        ),
        plot = P3_plot[[1]],
        width = 8,
        height = 7
      )
    }
    
    
    # 28、保存漏斗图
    
    if (length(P4_plot) > 0) {
      ggplot2::ggsave(
        filename = file.path(
          folder_name,
          "漏斗图.pdf"
        ),
        plot = P4_plot[[1]],
        width = 8,
        height = 7
      )
    }
    
    message(
      "暴露数据处理完成：",
      exposure_id
    )
    
    NA_character_
    
  }, error = function(e) {
    conditionMessage(e)
  })
  
  
  # 29、记录当前暴露的错误信息
  
  if (!is.na(processing_error)) {
    error_info <- data.frame(
      exposure_id = exposure_id,
      error_message = processing_error,
      stringsAsFactors = FALSE
    )
    
    error_log <- rbind(
      error_log,
      error_info
    )
    
    message(
      "处理暴露数据失败：", exposure_id,
      "；错误信息：", processing_error
    )
    
    write.csv(
      error_log,
      file = error_log_file,
      row.names = FALSE
    )
  }
}


# 30、保存最终汇总结果和错误日志

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


# 31、输出结果说明

message(
  "分析完成。\n",
  "final_results(exposure).csv中的P1-P5和b1-b5依次对应：\n",
  "1、MR Egger回归；只有一个SNP时为Wald ratio\n",
  "2、加权中位数法\n",
  "3、逆方差加权法\n",
  "4、简单模式法\n",
  "5、加权模式法\n",
  "主要结果通常以逆方差加权法为准。"
)