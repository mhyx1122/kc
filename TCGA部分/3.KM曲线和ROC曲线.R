# 1. 加载必要 R 包
library(data.table)
library(survival)
library(survminer)
library(gridExtra)
library(timeROC)

# 2. 参数设置

# 2.1 输入文件参数
input_file <- "合并后的标准肿瘤分析数据.csv"

# 2.2 KM 曲线参数
gene_names <- unlist(strsplit("FUT4, PIWIL4, CST1", ",\\s*"))

width1124 <- 8
height1124 <- 8

output_dir <- "6.单基因分析KM曲线"

# 2.3 ROC 曲线参数
gene_role <- "risk"
# gene_role <- "protective"

color_1 <- "blue"
color_2 <- "red"
color_3 <- "green"

roc_folder_name <- "7.ROC曲线"

# 3. 读取数据

data11221 <- fread(
  input_file,
  header = TRUE,
  sep = ",",
  check.names = FALSE,
  data.table = FALSE
)

rownames(data11221) <- data11221[, 1]
data11221 <- data11221[, -1]

# 4. KM 曲线分析与保存

if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

data <- data11221

for (gene_name in gene_names) {
  tryCatch({
    
    message("开始处理基因: ", gene_name)
    
    # 4.1 计算最佳截断点分组
    cutpoint <- surv_cutpoint(
      data,
      time = "futime",
      event = "fustat",
      variables = gene_name
    )
    
    cutpoint_value <- cutpoint$cutpoint$cutpoint
    message("最佳截断点计算完成，截断值: ", cutpoint_value)
    
    # 4.2 计算中位数分组
    median_value <- median(data[[gene_name]], na.rm = TRUE)
    message("中位数计算完成，值: ", median_value)
    
    # 4.3 生成两种分组方式
    data$cutpoint_group <- ifelse(data[[gene_name]] > cutpoint_value, "High", "Low")
    data$median_group <- ifelse(data[[gene_name]] > median_value, "High", "Low")
    
    # 4.4 KM 生存分析：最佳截断点
    surv_object <- Surv(data$futime, data$fustat)
    
    fit_cutpoint <- survfit(
      surv_object ~ cutpoint_group,
      data = data
    )
    
    km_plot_cutpoint <- ggsurvplot(
      fit_cutpoint,
      data = data,
      pval = TRUE,
      pval.method = TRUE,
      risk.table = TRUE,
      conf.int = TRUE,
      conf.int.style = "ribbon",
      title = paste(gene_name, " (Best Cutpoint)")
    )
    
    # 4.5 KM 生存分析：中位数分组
    fit_median <- survfit(
      surv_object ~ median_group,
      data = data
    )
    
    km_plot_median <- ggsurvplot(
      fit_median,
      data = data,
      pval = TRUE,
      pval.method = TRUE,
      risk.table = TRUE,
      conf.int = TRUE,
      conf.int.style = "ribbon",
      title = paste(gene_name, " (Median Split)")
    )
    
    # 4.6 保存最佳截断点 KM 曲线
    pdf_path_cutpoint <- file.path(
      output_dir,
      paste0(gene_name, "_KM_survival_Cutpoint.pdf")
    )
    
    pdf(
      pdf_path_cutpoint,
      width = width1124,
      height = height1124
    )
    
    grid.arrange(
      km_plot_cutpoint$plot,
      km_plot_cutpoint$table,
      nrow = 2,
      heights = c(3, 1)
    )
    
    dev.off()
    
    message("最佳截断点 K-M 生存曲线图已保存: ", pdf_path_cutpoint)
    
    # 4.7 保存中位数分组 KM 曲线
    pdf_path_median <- file.path(
      output_dir,
      paste0(gene_name, "_KM_survival_Median.pdf")
    )
    
    pdf(
      pdf_path_median,
      width = width1124,
      height = height1124
    )
    
    grid.arrange(
      km_plot_median$plot,
      km_plot_median$table,
      nrow = 2,
      heights = c(3, 1)
    )
    
    dev.off()
    
    message("中位数 K-M 生存曲线图已保存: ", pdf_path_median)
    
  }, error = function(e) {
    message("处理基因 ", gene_name, " 时发生错误: ", e$message, "\n")
  })
}

# 5. ROC 曲线分析与保存

target_genes <- gene_names
folder_name <- roc_folder_name
data <- data11221

if (!dir.exists(folder_name)) {
  dir.create(folder_name)
}

fustat <- data$fustat
futime <- data$futime

# 5.1 将生存时间从天转换为年
futime_years <- futime / 365

# 5.2 设置 ROC 预测时间点
times <- c(1, 3, 5)

# 5.3 循环绘制每个基因的 timeROC 曲线
for (gene in target_genes) {
  tryCatch({
    
    gene_expression <- data[[gene]]
    
    # 5.4 根据基因类型决定 marker 方向
    marker_use <- if (gene_role == "protective") {
      -gene_expression
    } else {
      gene_expression
    }
    
    # 5.5 计算 timeROC
    roc_result <- timeROC(
      T = futime_years,
      delta = fustat,
      marker = marker_use,
      cause = 1,
      weighting = "marginal",
      times = times,
      iid = TRUE
    )
    
    # 5.6 保存 ROC 曲线 PDF
    output_pdf <- paste0(folder_name, "/", gene, "_ROC_Curve.pdf")
    
    pdf(
      output_pdf,
      width = width1124,
      height = height1124
    )
    
    plot(
      roc_result,
      time = 1,
      col = color_1,
      title = paste("ROC Curves for", gene, "Gene Expression"),
      lwd = 2,
      xlim = c(0, 1),
      ylim = c(0, 1)
    )
    
    lines(
      roc_result$FP[, 2],
      roc_result$TP[, 2],
      col = color_2,
      lwd = 2
    )
    
    lines(
      roc_result$FP[, 3],
      roc_result$TP[, 3],
      col = color_3,
      lwd = 2
    )
    
    legend(
      "bottomright",
      legend = c(
        sprintf("1 Year AUC: %.2f", roc_result$AUC[1]),
        sprintf("3 Years AUC: %.2f", roc_result$AUC[2]),
        sprintf("5 Years AUC: %.2f", roc_result$AUC[3])
      ),
      col = c(color_1, color_2, color_3),
      lwd = 2,
      bty = "n",
      xjust = 1,
      yjust = 0,
      inset = c(0.05, 0.05)
    )
    
    dev.off()
    
    message("ROC生存曲线图已保存: ", output_pdf, "\n")
    
  }, error = function(e) {
    message("处理基因: ", gene, " 时发生错误，跳过到下一个\n")
  })
}
