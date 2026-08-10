# 1. 加载 R 包
library(data.table)
library(dplyr)
library(stringr)
library(survival)
library(limma)
library(survminer)
library(future)
library(future.apply)
library(timeROC)

options(shiny.maxRequestSize = 1 * 1024^3)


# 2. 模块一：表达矩阵与临床数据合并

# 2.1 模块一参数：表达矩阵文件
exp_file <- "after_group_TCGA.csv"

# 2.2 模块一参数：临床数据文件
cli_file <- "2.处理好的生存数据(因子型).csv"

# 2.3 模块一参数：是否进行 log2(x + 1) 转换
apply_log2 <- TRUE

# 2.4 模块一参数：数据类型
# 可选：
# "TCGA肿瘤样本"
# "其它肿瘤数据"
data_type <- "TCGA肿瘤样本"

# 2.5 模块一参数：输出目录
output_dir <- "2.生存分析筛选"

if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}


# 2.6 读取表达数据
message("模块一：开始读取表达数据")

if (grepl("\\.csv$", exp_file, ignore.case = TRUE)) {
  rt2 <- fread(exp_file, data.table = FALSE)
} else {
  rt2 <- fread(exp_file, sep = "\t", data.table = FALSE)
}

rownames(rt2) <- rt2[, 1]
rt2 <- rt2[, -1, drop = FALSE]

if (apply_log2) {
  rt2 <- log2(rt2 + 1)
}


# 2.7 读取临床数据
message("模块一：开始读取临床数据")

if (grepl("\\.csv$", cli_file, ignore.case = TRUE)) {
  cli <- fread(cli_file, data.table = FALSE)
} else {
  cli <- fread(cli_file, sep = "\t", data.table = FALSE)
}

rownames(cli) <- cli[, 1]
cli <- cli[, -1, drop = FALSE]


# 2.8 合并表达数据和临床数据
message("模块一：开始合并表达数据和临床数据")

if (data_type == "TCGA肿瘤样本") {
  
  # 2.8.1 提取 TCGA 肿瘤样本，删除正常样品
  exp_data_T <- rt2 %>%
    select(which(as.numeric(str_sub(colnames(rt2), 14, 15)) < 10))
  
  tumorData <- as.matrix(exp_data_T)
  tumorData <- t(tumorData)
  rownames(tumorData) <- substr(rownames(tumorData), 1, 12)
  
  data <- avereps(tumorData)
  
  rownames(cli) <- gsub("-", ".", rownames(cli))
  
  # 2.8.2 合并表达数据和生存数据
  sameSample <- intersect(row.names(data), row.names(cli))
  data <- data[sameSample, , drop = FALSE]
  cli <- cli[sameSample, , drop = FALSE]
  
  rt3 <- cbind(cli, data)
  
  fwrite(
    data.frame(Sample = rownames(rt3), rt3, check.names = FALSE),
    na = "NA",
    file = file.path(output_dir, "合并后的标准肿瘤分析数据.csv")
  )
  
} else {
  
  # 2.8.3 其它肿瘤数据直接转置表达矩阵
  exp_data_T <- rt2
  
  tumorData <- as.matrix(exp_data_T)
  tumorData <- t(tumorData)
  
  data <- avereps(tumorData)
  
  # 2.8.4 合并表达数据和生存数据
  sameSample <- intersect(row.names(data), row.names(cli))
  data <- data[sameSample, , drop = FALSE]
  cli <- cli[sameSample, , drop = FALSE]
  
  rt3 <- cbind(cli, data)
  
  fwrite(
    data.frame(Sample = rownames(rt3), rt3, check.names = FALSE),
    na = "NA",
    file = file.path(output_dir, "合并后的标准肿瘤分析数据.csv")
  )
}

message("模块一完成：合并后的标准肿瘤分析数据.csv 已保存")


# 3. 模块二：生存分析筛选

# 3.1 模块二参数：输入文件
merged_file <- file.path(output_dir, "合并后的标准肿瘤分析数据.csv")

# 3.2 模块二参数：基因起始列
start_col_index <- 9

# 3.3 模块二参数：并行核心数
numCores <- 6

# 3.4 模块二参数：低标准差基因过滤阈值
sd_threshold <- 0.1

# 3.5 模块二参数：P 值筛选阈值
surv_cutoff <- 0.05
pval_cutoff <- surv_cutoff

# 3.6 模块二参数：输出目录
output_dir <- "2.生存分析筛选"

if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}


# 3.7 读取合并后的生存分析数据
message("模块二：开始读取合并后的生存分析数据")

plan(multisession, workers = numCores)

if (grepl("\\.csv$", merged_file, ignore.case = TRUE)) {
  rt2 <- fread(merged_file, data.table = FALSE)
} else {
  rt2 <- fread(merged_file, sep = "\t", data.table = FALSE)
}

rownames(rt2) <- rt2[, 1]
rt2 <- rt2[, -1, drop = FALSE]

message("计算过程较慢，耐心等待")


# 3.8 最佳截断值 KM 分析
message("模块二：开始最佳截断值 KM 分析")

outTab <- future_lapply(colnames(rt2[, start_col_index:ncol(rt2), drop = FALSE]), function(gene) {
  
  if (sd(rt2[, gene], na.rm = TRUE) < sd_threshold) return(NULL)
  
  tryCatch({
    cutpoint <- surv_cutpoint(rt2, time = "futime", event = "fustat", variables = gene)
    rt2$group <- ifelse(rt2[[gene]] > cutpoint$cutpoint$cutpoint, "High", "Low")
    
    surv_object <- Surv(rt2$futime, rt2$fustat)
    fit <- survfit(surv_object ~ group, data = rt2)
    diff_test <- survdiff(surv_object ~ group, data = rt2)
    p_value <- 1 - pchisq(diff_test$chisq, length(diff_test$n) - 1)
    
    rt2$group <- NULL
    
    cbind(gene = gene, KM = p_value)
  }, error = function(e) {
    NULL
  })
}, future.seed = TRUE)

outTab <- do.call(rbind, outTab)
outTab <- as.data.frame(outTab)
outTab$KM <- as.numeric(as.character(outTab$KM))

write.csv(
  outTab,
  file = file.path(output_dir, "survival(最佳截断值全结果).csv"),
  row.names = FALSE,
  quote = FALSE
)

filtered_outTabKM <- outTab[, c(1, 2), drop = FALSE]
filtered_outTabKM <- filtered_outTabKM[filtered_outTabKM$KM < pval_cutoff, , drop = FALSE]

write.csv(
  filtered_outTabKM,
  file = file.path(output_dir, "survival(KM最佳截断值,P值过滤后).csv"),
  row.names = FALSE,
  quote = FALSE
)

message("已完成：最佳截断值 KM 显著基因筛选")


# 3.9 中位数分组 KM 分析
message("模块二：开始中位数分组 KM 分析")

outTab1 <- future_lapply(colnames(rt2[, start_col_index:ncol(rt2), drop = FALSE]), function(gene) {
  
  if (sd(rt2[, gene], na.rm = TRUE) < sd_threshold) return(NULL)
  
  tryCatch({
    median_value <- median(rt2[[gene]], na.rm = TRUE)
    rt2$group <- ifelse(rt2[[gene]] > median_value, "High", "Low")
    
    surv_object <- Surv(rt2$futime, rt2$fustat)
    fit <- survfit(surv_object ~ group, data = rt2)
    diff_test <- survdiff(surv_object ~ group, data = rt2)
    p_value <- 1 - pchisq(diff_test$chisq, length(diff_test$n) - 1)
    
    rt2$group <- NULL
    
    cbind(gene = gene, KM = p_value)
  }, error = function(e) {
    NULL
  })
}, future.seed = TRUE)

outTab1 <- do.call(rbind, outTab1)
outTab1 <- as.data.frame(outTab1)
outTab1$KM <- as.numeric(as.character(outTab1$KM))

write.csv(
  outTab1,
  file = file.path(output_dir, "survival(中位数全结果).csv"),
  row.names = FALSE,
  quote = FALSE
)

filtered_outTabKM1 <- outTab1[, c(1, 2), drop = FALSE]
filtered_outTabKM1 <- filtered_outTabKM1[filtered_outTabKM1$KM < pval_cutoff, , drop = FALSE]

write.csv(
  filtered_outTabKM1,
  file = file.path(output_dir, "survival(KM中位数分组,P值过滤后).csv"),
  row.names = FALSE,
  quote = FALSE
)

message("已完成：中位数分组 KM 显著基因筛选")


# 3.10 单因素 Cox 分析
message("模块二：开始单因素 Cox 分析")

uniCoxResults <- future_lapply(colnames(rt2[, start_col_index:ncol(rt2), drop = FALSE]), function(i) {
  
  if (sd(rt2[, i], na.rm = TRUE) < sd_threshold) return(NULL)
  
  tryCatch({
    cox <- suppressWarnings(coxph(Surv(futime, fustat) ~ rt2[, i], data = rt2))
    coxSummary <- summary(cox)
    
    data.frame(
      cox_gene = i,
      HR = coxSummary$conf.int[, "exp(coef)"],
      HR.95L = coxSummary$conf.int[, "lower .95"],
      HR.95H = coxSummary$conf.int[, "upper .95"],
      pvalue = coxSummary$coefficients[, "Pr(>|z|)"]
    )
  }, error = function(e) {
    NULL
  })
}, future.seed = TRUE)

uniCoxResults <- do.call(rbind, uniCoxResults)
uniCoxResults <- as.data.frame(uniCoxResults)

write.csv(
  uniCoxResults,
  file = file.path(output_dir, "uniCoxResults(单基因COX).csv"),
  row.names = FALSE
)

filtered_uniCoxResults <- uniCoxResults[uniCoxResults$pvalue < pval_cutoff, , drop = FALSE]

write.csv(
  filtered_uniCoxResults,
  file = file.path(output_dir, "uniCoxResults(单基因COX,P值过滤后).csv"),
  row.names = FALSE,
  quote = FALSE
)

plan(sequential)
suppressWarnings(gc())

message("已完成：单因素 Cox 显著基因筛选")


# 3.11 汇总 KM 和 Cox 共同显著基因
message("模块二：开始汇总共同显著基因")

common_genes_km_cut_cox <- intersect(
  filtered_outTabKM$gene,
  filtered_uniCoxResults$cox_gene
)

common_result_km_cut_cox <- merge(
  filtered_outTabKM,
  filtered_uniCoxResults,
  by.x = "gene",
  by.y = "cox_gene"
)

write.csv(
  common_result_km_cut_cox,
  file = file.path(output_dir, "KM最佳截断值+COX显著的基因.csv"),
  row.names = FALSE,
  quote = FALSE
)

common_genes_km_median_cox <- intersect(
  filtered_outTabKM1$gene,
  filtered_uniCoxResults$cox_gene
)

common_result_km_median_cox <- merge(
  filtered_outTabKM1,
  filtered_uniCoxResults,
  by.x = "gene",
  by.y = "cox_gene"
)

write.csv(
  common_result_km_median_cox,
  file = file.path(output_dir, "common_genes_KM中位数分组_and_COX.csv"),
  row.names = FALSE,
  quote = FALSE
)

plan(sequential)
suppressWarnings(gc())

message("模块二完成")
message(paste0("最佳截断值 KM + Cox 共同显著基因数：", length(common_genes_km_cut_cox)))
message(paste0("中位数 KM + Cox 共同显著基因数：", length(common_genes_km_median_cox)))


# 4. 模块三：timeROC 生存 ROC 分析

# 4.1 模块三参数：输入文件
ROC_file <- file.path(output_dir, "合并后的标准肿瘤分析数据.csv")

# 4.2 模块三参数：基因起始列
start_col <- 9

# 4.3 模块三参数：ROC 时间点，单位为年
times <- c(1, 3, 5)

# 4.4 模块三参数：输出目录
output_dir <- "2.生存分析筛选"

if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}


# 4.5 读取 ROC 输入数据
message("模块三：开始读取 timeROC 输入数据")

datacgxr410 <- fread(ROC_file, header = TRUE, data.table = FALSE)

ngenestar <- start_col

rownames(datacgxr410) <- datacgxr410[[1]]
datacgxr410 <- as.data.frame(datacgxr410)

data <- datacgxr410[, -1]

fustat <- data$fustat
futime <- data$futime

futime_years <- futime / 365


# 4.6 逐基因进行 timeROC 分析
message("模块三：开始逐基因 timeROC 分析")

auc_results <- data.frame(
  Gene = character(),
  AUC_1yr = numeric(),
  AUC_3yr = numeric(),
  AUC_5yr = numeric(),
  Total_AUC = numeric(),
  stringsAsFactors = FALSE
)

total_genes <- ncol(data) - ngenestar + 1
counter <- 0

for (i in ngenestar:ncol(data)) {
  
  counter <- counter + 1
  
  gene_expression <- data[[i]]
  
  roc_result <- tryCatch({
    timeROC(
      T = futime_years,
      delta = fustat,
      marker = gene_expression,
      cause = 1,
      weighting = "marginal",
      times = times,
      iid = TRUE
    )
  }, error = function(e) {
    NULL
  })
  
  if (!is.null(roc_result)) {
    
    auc_1yr <- roc_result$AUC[which(roc_result$times == 1)]
    auc_3yr <- roc_result$AUC[which(roc_result$times == 3)]
    auc_5yr <- roc_result$AUC[which(roc_result$times == 5)]
    
    auc_value <- roc_result$AUC[1]
    
    gene_name <- colnames(data)[i]
    
    auc_results <- rbind(
      auc_results,
      data.frame(
        Gene = gene_name,
        AUC_1yr = auc_1yr,
        AUC_3yr = auc_3yr,
        AUC_5yr = auc_5yr,
        Total_AUC = auc_value
      )
    )
  }
  
  if (counter %% 500 == 0) {
    cat("Processed", counter, "genes...\n")
  }
}

auc_results$AUC_1yr <- as.numeric(as.character(auc_results$AUC_1yr))
auc_results$AUC_3yr <- as.numeric(as.character(auc_results$AUC_3yr))
auc_results$AUC_5yr <- as.numeric(as.character(auc_results$AUC_5yr))
auc_results$Total_AUC <- as.numeric(as.character(auc_results$Total_AUC))


# 4.7 按 AUC 方向区分风险基因和保护基因
message("模块三：开始区分风险基因和保护基因")

risk_auc_results <- auc_results[auc_results$Total_AUC >= 0.5, , drop = FALSE]
risk_auc_results$Type <- "Risk"

protective_auc_results <- auc_results[auc_results$Total_AUC < 0.5, , drop = FALSE]

protective_auc_results$AUC_1yr <- 1 - protective_auc_results$AUC_1yr
protective_auc_results$AUC_3yr <- 1 - protective_auc_results$AUC_3yr
protective_auc_results$AUC_5yr <- 1 - protective_auc_results$AUC_5yr
protective_auc_results$Total_AUC <- 1 - protective_auc_results$Total_AUC
protective_auc_results$Type <- "Protective"

risk_auc_results <- risk_auc_results[order(-risk_auc_results$Total_AUC), , drop = FALSE]
protective_auc_results <- protective_auc_results[order(-protective_auc_results$Total_AUC), , drop = FALSE]


# 4.8 保存 ROC 结果
write.csv(
  auc_results,
  file.path(output_dir, "gene_auc_results.csv"),
  row.names = FALSE
)

write.csv(
  risk_auc_results,
  file.path(output_dir, "gene_auc_results_risk.csv"),
  row.names = FALSE
)

write.csv(
  protective_auc_results,
  file.path(output_dir, "gene_auc_results_protective.csv"),
  row.names = FALSE
)

message("模块三完成：timeROC 结果已保存")
message("全部分析完成")