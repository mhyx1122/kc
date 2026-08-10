# 1. 加载R包

suppressPackageStartupMessages({
  library(TwoSampleMR)
})


# 2. 读取暴露ID文件

# 暴露ID文件路径，CSV文件不设置表头，第一列为暴露ID
exposure_file <- "暴露ID.csv"

exp_data <- read.csv(exposure_file, header = FALSE, check.names = FALSE, stringsAsFactors = FALSE)

if (ncol(exp_data) < 1) {
  stop("暴露ID文件中没有可用的列。")
}

exposure <- trimws(as.character(exp_data[, 1]))
exposure <- exposure[!is.na(exposure) & exposure != ""]

if (length(exposure) == 0) {
  stop("暴露ID文件中没有可用的暴露ID。")
}


# 3. 设置结局数据及暴露工具变量筛选参数

# 结局ID
outcomes <- "ieu-a-1165"

# 暴露工具变量的P值阈值
p1 <- 1e-5

# 是否进行Clumping
clump <- TRUE

# LD Clumping的R2阈值
r2 <- 0.001

# LD Clumping的窗口范围，单位为kb
kb <- 10000


# 4. 设置输出路径并初始化结果表

# 每个暴露的详细结果保存目录
output_root <- "原始数据"

# 汇总结果文件
final_results_file <- "final_results(exposure).csv"

# 错误日志文件
error_log_file <- "error_log.csv"

if (!dir.exists(output_root)) {
  dir.create(output_root, recursive = TRUE)
}

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


# 5. 循环分析每一个暴露ID

for (exposure_id in exposure) {
  
  message("\n开始处理暴露：", exposure_id)
  
  processing_error <- tryCatch({
    
    # 5.1 创建当前暴露的结果文件夹
    
    folder_name <- file.path(output_root, exposure_id)
    
    if (!dir.exists(folder_name)) {
      dir.create(folder_name, recursive = TRUE)
    }
    
    
    # 5.2 从OpenGWAS提取暴露工具变量
    
    exposure_dat <- suppressWarnings(
      extract_instruments(outcomes = exposure_id, p1 = p1, clump = clump, r2 = r2, kb = kb)
    )
    
    if (is.null(exposure_dat) || nrow(exposure_dat) == 0) {
      stop("没有提取到符合条件的暴露SNP。")
    }
    
    
    # 5.3 计算每个SNP的R2和F统计量
    
    f_required_columns <- c("eaf.exposure", "beta.exposure", "samplesize.exposure")
    
    if (all(f_required_columns %in% colnames(exposure_dat))) {
      
      exposure_dat$R2 <- 2 * (1 - exposure_dat$eaf.exposure) *
        exposure_dat$eaf.exposure * exposure_dat$beta.exposure^2
      
      exposure_dat$F <- exposure_dat$R2 / (1 - exposure_dat$R2) *
        (exposure_dat$samplesize.exposure - 2)
      
      if (any(!is.na(exposure_dat$F))) {
        exposure_dat <- subset(exposure_dat, F > 10)
      } else {
        message("暴露 ", exposure_id, " 缺少可用样本量，跳过F值过滤。")
      }
      
    } else {
      
      exposure_dat$R2 <- NA_real_
      exposure_dat$F <- NA_real_
      
      message("暴露 ", exposure_id, " 缺少计算F值所需的列，跳过F值计算和过滤。")
    }
    
    if (nrow(exposure_dat) == 0) {
      stop("经过F > 10过滤后没有剩余SNP。")
    }
    
    write.csv(
      exposure_dat,
      file.path(folder_name, "exposure_dat-SNPs-F.csv"),
      row.names = FALSE
    )
    
    
    # 5.4 提取结局数据
    
    outcome_dat <- extract_outcome_data(snps = exposure_dat$SNP, outcomes = outcomes)
    
    if (is.null(outcome_dat) || nrow(outcome_dat) == 0) {
      stop("没有提取到对应的结局数据。")
    }
    
    
    # 5.5 协调暴露数据和结局数据
    
    dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
    
    write.csv(dat, file.path(folder_name, "dat_results.csv"), row.names = FALSE)
    
    if (!"mr_keep" %in% colnames(dat) || sum(dat$mr_keep, na.rm = TRUE) == 0) {
      stop("数据协调后没有可用于MR分析的SNP。")
    }
    
    
    # 5.6 进行孟德尔随机化分析
    
    res <- suppressWarnings(mr(dat))
    
    if (is.null(res) || nrow(res) == 0) {
      stop("孟德尔随机化分析没有返回有效结果。")
    }
    
    res_or <- generate_odds_ratios(res)
    
    write.csv(res_or, file.path(folder_name, "mr-result.csv"), row.names = FALSE)
    
    print(res)
    
    
    # 5.7 进行异质性检验
    
    het <- tryCatch(
      suppressWarnings(mr_heterogeneity(dat)),
      error = function(e) {
        message("暴露 ", exposure_id, " 的异质性检验无法完成：", conditionMessage(e))
        data.frame()
      }
    )
    
    if (is.null(het)) {
      het <- data.frame()
    }
    
    write.csv(het, file.path(folder_name, "heterogeneity.csv"), row.names = FALSE)
    
    
    # 5.8 进行水平多效性检验
    
    pleio <- tryCatch(
      suppressWarnings(mr_pleiotropy_test(dat)),
      error = function(e) {
        message("暴露 ", exposure_id, " 的水平多效性检验无法完成：", conditionMessage(e))
        data.frame()
      }
    )
    
    if (is.null(pleio)) {
      pleio <- data.frame()
    }
    
    write.csv(pleio, file.path(folder_name, "pleiotropy_test.csv"), row.names = FALSE)
    
    
    # 5.9 按照方法名称提取五种MR方法的位置
    
    egger_index <- match("MR Egger", res$method)
    wald_index <- match("Wald ratio", res$method)
    weighted_median_index <- match("Weighted median", res$method)
    ivw_index <- match("Inverse variance weighted", res$method)
    simple_mode_index <- match("Simple mode", res$method)
    weighted_mode_index <- match("Weighted mode", res$method)
    
    if (!is.na(egger_index)) {
      method1_index <- egger_index
    } else {
      method1_index <- wald_index
    }
    
    
    # 5.10 提取五种方法的P值
    
    pval_1 <- if (!is.na(method1_index)) res$pval[method1_index] else NA_real_
    pval_2 <- if (!is.na(weighted_median_index)) res$pval[weighted_median_index] else NA_real_
    pval_3 <- if (!is.na(ivw_index)) res$pval[ivw_index] else NA_real_
    pval_4 <- if (!is.na(simple_mode_index)) res$pval[simple_mode_index] else NA_real_
    pval_5 <- if (!is.na(weighted_mode_index)) res$pval[weighted_mode_index] else NA_real_
    
    
    # 5.11 提取五种方法的效应值
    
    beta_1 <- if (!is.na(method1_index)) res$b[method1_index] else NA_real_
    beta_2 <- if (!is.na(weighted_median_index)) res$b[weighted_median_index] else NA_real_
    beta_3 <- if (!is.na(ivw_index)) res$b[ivw_index] else NA_real_
    beta_4 <- if (!is.na(simple_mode_index)) res$b[simple_mode_index] else NA_real_
    beta_5 <- if (!is.na(weighted_mode_index)) res$b[weighted_mode_index] else NA_real_
    
    
    # 5.12 提取MR-Egger和IVW异质性检验P值
    
    heterogeneity_egger <- NA_real_
    heterogeneity_ivw <- NA_real_
    
    if (nrow(het) > 0 && all(c("method", "Q_pval") %in% colnames(het))) {
      
      het_egger_index <- match("MR Egger", het$method)
      het_ivw_index <- match("Inverse variance weighted", het$method)
      
      if (!is.na(het_egger_index)) {
        heterogeneity_egger <- het$Q_pval[het_egger_index]
      }
      
      if (!is.na(het_ivw_index)) {
        heterogeneity_ivw <- het$Q_pval[het_ivw_index]
      }
    }
    
    
    # 5.13 提取水平多效性检验P值
    
    pleiotropy_pvalue <- NA_real_
    
    if (nrow(pleio) > 0 && "pval" %in% colnames(pleio)) {
      pleiotropy_pvalue <- pleio$pval[1]
    }
    
    
    # 5.14 提取IVW或Wald ratio的OR和95%置信区间
    
    or_ivw_index <- match("Inverse variance weighted", res_or$method)
    or_wald_index <- match("Wald ratio", res_or$method)
    
    if (!is.na(or_ivw_index)) {
      main_or_index <- or_ivw_index
    } else {
      main_or_index <- or_wald_index
    }
    
    if (!is.na(main_or_index)) {
      main_or <- res_or$or[main_or_index]
      main_or_lci95 <- res_or$or_lci95[main_or_index]
      main_or_uci95 <- res_or$or_uci95[main_or_index]
    } else {
      main_or <- NA_real_
      main_or_lci95 <- NA_real_
      main_or_uci95 <- NA_real_
    }
    
    
    # 5.15 将当前暴露的结果添加到汇总表
    
    new_row <- data.frame(
      id.exposure = as.character(res$id.exposure[1]),
      nSNP = as.numeric(res$nsnp[1]),
      P1 = pval_1,
      P2 = pval_2,
      P3 = pval_3,
      P4 = pval_4,
      P5 = pval_5,
      b1 = beta_1,
      b2 = beta_2,
      b3 = beta_3,
      b4 = beta_4,
      b5 = beta_5,
      id.outcome = as.character(res$id.outcome[1]),
      heterogeneity1 = heterogeneity_egger,
      heterogeneity2 = heterogeneity_ivw,
      pleiotropy = pleiotropy_pvalue,
      or = main_or,
      or_lci95 = main_or_lci95,
      or_uci95 = main_or_uci95,
      exposure_name = as.character(res$exposure[1]),
      stringsAsFactors = FALSE
    )
    
    final_results <- rbind(final_results, new_row)
    
    write.csv(final_results, final_results_file, row.names = FALSE)
    
    
    # 5.16 进行逐个剔除检验并保存图形
    
    leaveoneout_error <- tryCatch({
      
      res_loo <- suppressWarnings(mr_leaveoneout(dat))
      leaveoneout_plot_list <- suppressWarnings(mr_leaveoneout_plot(res_loo))
      
      if (length(leaveoneout_plot_list) > 0) {
        pdf(file.path(folder_name, "剔除检验.pdf"), width = 8, height = 7)
        print(leaveoneout_plot_list[[1]])
        dev.off()
      }
      
      NULL
      
    }, error = function(e) {
      
      if (dev.cur() > 1) {
        dev.off()
      }
      
      e
    })
    
    if (inherits(leaveoneout_error, "error")) {
      message("暴露 ", exposure_id, " 的逐个剔除检验或绘图失败：",
              conditionMessage(leaveoneout_error))
    }
    
    
    # 5.17 保存MR效应散点图
    
    scatter_error <- tryCatch({
      
      scatter_plot_list <- suppressWarnings(mr_scatter_plot(res, dat))
      
      if (length(scatter_plot_list) > 0) {
        pdf(file.path(folder_name, "效应散点图.pdf"), width = 8, height = 7)
        print(scatter_plot_list[[1]])
        dev.off()
      }
      
      NULL
      
    }, error = function(e) {
      
      if (dev.cur() > 1) {
        dev.off()
      }
      
      e
    })
    
    if (inherits(scatter_error, "error")) {
      message("暴露 ", exposure_id, " 的效应散点图生成失败：",
              conditionMessage(scatter_error))
    }
    
    
    # 5.18 进行单SNP分析
    
    single_snp_result <- tryCatch(
      suppressWarnings(mr_singlesnp(dat)),
      error = function(e) e
    )
    
    
    # 5.19 保存单SNP森林图
    
    if (!inherits(single_snp_result, "error")) {
      
      forest_error <- tryCatch({
        
        forest_plot_list <- suppressWarnings(mr_forest_plot(single_snp_result))
        
        if (length(forest_plot_list) > 0) {
          pdf(file.path(folder_name, "效应森林图.pdf"), width = 8, height = 7)
          print(forest_plot_list[[1]])
          dev.off()
        }
        
        NULL
        
      }, error = function(e) {
        
        if (dev.cur() > 1) {
          dev.off()
        }
        
        e
      })
      
      if (inherits(forest_error, "error")) {
        message("暴露 ", exposure_id, " 的效应森林图生成失败：",
                conditionMessage(forest_error))
      }
      
      
      # 5.20 保存漏斗图
      
      funnel_error <- tryCatch({
        
        funnel_plot_list <- suppressWarnings(mr_funnel_plot(single_snp_result))
        
        if (length(funnel_plot_list) > 0) {
          pdf(file.path(folder_name, "漏斗图.pdf"), width = 8, height = 7)
          print(funnel_plot_list[[1]])
          dev.off()
        }
        
        NULL
        
      }, error = function(e) {
        
        if (dev.cur() > 1) {
          dev.off()
        }
        
        e
      })
      
      if (inherits(funnel_error, "error")) {
        message("暴露 ", exposure_id, " 的漏斗图生成失败：",
                conditionMessage(funnel_error))
      }
      
    } else {
      
      message("暴露 ", exposure_id, " 的单SNP分析失败：",
              conditionMessage(single_snp_result))
    }
    
    message("暴露 ", exposure_id, " 处理完成。")
    
    NULL
    
  }, error = function(e) {
    e
  })
  
  
  # 5.21 记录当前暴露的错误信息
  
  if (inherits(processing_error, "error")) {
    
    error_info <- data.frame(
      exposure_id = as.character(exposure_id),
      error_message = conditionMessage(processing_error),
      stringsAsFactors = FALSE
    )
    
    error_log <- rbind(error_log, error_info)
    
    write.csv(error_log, error_log_file, row.names = FALSE)
    
    message("暴露 ", exposure_id, " 处理失败：", conditionMessage(processing_error))
  }
}


# 6. 保存最终汇总结果和错误日志

write.csv(final_results, final_results_file, row.names = FALSE)
write.csv(error_log, error_log_file, row.names = FALSE)


# 7. 输出结果说明

message(
  "\n分析完成。\n",
  "final_results(exposure).csv中的P1-P5和b1-b5分别对应：\n",
  "1. MR Egger；只有一个SNP时对应Wald ratio\n",
  "2. Weighted median\n",
  "3. Inverse variance weighted\n",
  "4. Simple mode\n",
  "5. Weighted mode\n",
  "主要结果通常以Inverse variance weighted方法为主；",
  "只有一个SNP时使用Wald ratio结果。\n",
  "详细结果保存在：", output_root
)