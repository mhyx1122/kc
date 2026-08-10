suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
  library(ggplot2)
})

# 1、读取暴露ID文件

# 暴露ID文件路径，第一列为OpenGWAS暴露ID，不包含表头
exposure_id_file <- "暴露ID.csv"

if (!file.exists(exposure_id_file)) {
  stop("没有找到暴露ID文件：", exposure_id_file)
}

exp_data <- read.csv(
  exposure_id_file,
  header = FALSE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

exposure <- trimws(as.character(exp_data[[1]]))
exposure <- exposure[!is.na(exposure) & exposure != ""]

if (length(exposure) == 0) {
  stop("暴露ID文件第一列中没有有效的暴露ID。")
}


# 2、读取本地结局数据

# 已通过TwoSampleMR格式化的结局数据CSV文件
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


# 3、设置在线提取暴露数据的参数

# 暴露工具变量的P值阈值
p1 <- 1e-5

# 是否进行Clumping
clump <- TRUE

# LD Clumping的R²阈值
r2 <- 0.001

# LD Clumping的窗口大小，单位为kb
kb <- 10000


# 4、设置结果保存路径

# 每个暴露对应的结果保存总文件夹
output_folder <- "原始数据"

# 汇总结果文件
final_result_file <- "final_results(exposure).csv"

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
  exposure_name = character(),
  stringsAsFactors = FALSE
)


# 6、初始化错误日志

error_log <- data.frame(
  exposure_id = character(),
  error_message = character(),
  stringsAsFactors = FALSE
)


# 7、循环处理每个暴露ID

for (i in seq_along(exposure)) {
  exposure_id <- exposure[i]
  
  message(
    "[", i, "/", length(exposure),
    "] 正在处理暴露ID：", exposure_id
  )
  
  processing_error <- tryCatch({
    
    # 8、创建当前暴露的结果文件夹
    
    folder_name <- file.path(output_folder, exposure_id)
    
    if (!dir.exists(folder_name)) {
      dir.create(folder_name, recursive = TRUE)
    }
    
    
    # 9、在线提取暴露GWAS工具变量
    
    exposure_dat <- suppressWarnings(
      TwoSampleMR::extract_instruments(
        outcomes = exposure_id,
        p1 = p1,
        clump = clump,
        r2 = r2,
        kb = kb
      )
    )
    
    if (is.null(exposure_dat) || nrow(exposure_dat) == 0) {
      stop("没有提取到满足条件的暴露工具变量。")
    }
    
    required_exposure_cols <- c(
      "SNP",
      "beta.exposure",
      "se.exposure",
      "effect_allele.exposure",
      "other_allele.exposure",
      "eaf.exposure",
      "samplesize.exposure",
      "exposure",
      "id.exposure"
    )
    
    missing_exposure_cols <- setdiff(
      required_exposure_cols,
      colnames(exposure_dat)
    )
    
    if (length(missing_exposure_cols) > 0) {
      stop(
        "提取的暴露数据缺少以下必要列：",
        paste(missing_exposure_cols, collapse = "、")
      )
    }
    
    
    # 10、计算R2和F值
    
    exposure_dat$R2 <- 2 *
      (1 - exposure_dat$eaf.exposure) *
      exposure_dat$eaf.exposure *
      (exposure_dat$beta.exposure)^2
    
    exposure_dat$F <- exposure_dat$R2 /
      (1 - exposure_dat$R2) *
      (exposure_dat$samplesize.exposure - 2)
    
    
    # 11、筛选F值大于10的工具变量
    
    if (any(!is.na(exposure_dat$F))) {
      exposure_dat <- subset(exposure_dat, F > 10)
    } else {
      message(
        "当前暴露缺少有效的样本量或等位基因频率，",
        "跳过F值筛选：", exposure_id
      )
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
    
    
    # 12、从本地结局数据中提取共有SNP
    
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
    
    
    # 13、协调暴露数据和结局数据
    
    dat <- TwoSampleMR::harmonise_data(
      exposure_dat = exposure_dat,
      outcome_dat = outcome_datkkk,
      action = 2
    )
    
    write.csv(
      dat,
      file = file.path(
        folder_name,
        "dat_results.csv"
      ),
      row.names = FALSE
    )
    
    if (
      nrow(dat) == 0 ||
      !"mr_keep" %in% colnames(dat) ||
      !any(dat$mr_keep, na.rm = TRUE)
    ) {
      stop("数据协调后没有可用于MR分析的SNP。")
    }
    
    
    # 14、执行孟德尔随机化分析
    
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
    
    
    # 15、执行异质性检验
    
    het <- TwoSampleMR::mr_heterogeneity(dat)
    
    write.csv(
      het,
      file = file.path(
        folder_name,
        "heterogeneity.csv"
      ),
      row.names = FALSE
    )
    
    
    # 16、执行水平多效性检验
    
    pleio <- TwoSampleMR::mr_pleiotropy_test(dat)
    
    write.csv(
      pleio,
      file = file.path(
        folder_name,
        "pleiotropy_test.csv"
      ),
      row.names = FALSE
    )
    
    
    # 17、按照原代码顺序提取五种MR方法结果
    
    p_values <- res$pval[1:5]
    b_values <- res$b[1:5]
    heterogeneity_values <- het$Q_pval[1:2]
    
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
      pleiotropy = pleio$pval[1],
      or = res1$or[3],
      or_lci95 = res1$or_lci95[3],
      or_uci95 = res1$or_uci95[3],
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
    
    
    # 18、执行逐个剔除检验和单SNP分析
    
    res_loo <- TwoSampleMR::mr_leaveoneout(dat)
    res_single <- TwoSampleMR::mr_singlesnp(dat)
    
    
    # 19、生成敏感性分析图
    
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
    
    P4_plot <- TwoSampleMR::mr_funnel_plot(
      res_single
    )
    
    
    # 20、保存逐个剔除检验图
    
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
    
    
    # 21、保存效应散点图
    
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
    
    
    # 22、保存效应森林图
    
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
    
    
    # 23、保存漏斗图
    
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
      "暴露ID处理完成：",
      exposure_id
    )
    
    NA_character_
    
  }, error = function(e) {
    conditionMessage(e)
  })
  
  
  # 24、记录当前暴露的错误信息
  
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
      "处理暴露ID失败：", exposure_id,
      "；错误信息：", processing_error
    )
    
    write.csv(
      error_log,
      file = error_log_file,
      row.names = FALSE
    )
  }
}


# 25、保存最终汇总结果和错误日志

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


# 26、输出结果说明

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