library(ggplot2)
library(aPEAR)

# 1. 设置输入参数
input_file <- "your_enrichment_result.csv"
pvalue_cutoff <- 0.05
analysis_type <- "KEGG"
color_low <- "#E64B35FF"
color_high <- "#4DBBD5FF"
output_dir <- "定制图片"
plot_width <- 8
plot_height <- 6

# 2. 创建输出文件夹
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# 3. 读取CSV文件并筛选P值
data <- read.csv(input_file,header = TRUE)
data_sig <- data[data$pvalue <= pvalue_cutoff,]

# 4. 生成KEGG富集网络图
if (analysis_type == "KEGG") {
  p1 <- enrichmentNetwork(
    data_sig,
    colorBy = "pvalue",
    colorType = c("pval"),
    nodeSize = "Count",
    fontSize = 4,
    drawEllipses = TRUE,
    pCutoff = -20,
    verbose = TRUE
  )
  p1 <- p1 +
    theme(legend.key = element_blank()) +
    scale_color_gradient(low = color_low,high = color_high) +
    scale_fill_gradient(low = color_low,high = color_high)
  output_file_p1 <- file.path(output_dir,"KEGG_enrichment_network.pdf")
  ggsave(output_file_p1,plot = p1,width = plot_width,height = plot_height)
  message("KEGG富集网络图已保存到：",output_file_p1)
}

# 5. 生成GSEA富集网络图
if (analysis_type == "GSEA") {
  p2 <- enrichmentNetwork(
    data_sig,
    colorBy = "nes",
    colorType = c("nes"),
    nodeSize = "setSize",
    fontSize = 4,
    drawEllipses = TRUE,
    pCutoff = -20,
    verbose = TRUE
  )
  p2 <- p2 +
    theme(legend.key = element_blank()) +
    scale_color_gradient(low = color_low,high = color_high) +
    scale_fill_gradient(low = color_low,high = color_high)
  output_file_p2 <- file.path(output_dir,"GSEA_enrichment_network.pdf")
  ggsave(output_file_p2,plot = p2,width = plot_width,height = plot_height)
  message("GSEA富集网络图已保存到：",output_file_p2)
}