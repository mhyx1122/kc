library(ggVolcano)
library(ggplot2)
library(RColorBrewer)

# 1. 随机生成颜色的函数
randomColor <- function() {
  paste0("#", paste0(sample(c(0:9, letters[1:6]), 6, replace = TRUE), collapse = ""))
}

# 2. 设置输入文件路径
file_path <- "差异分析结果文件.csv"
genename_file <- "感兴趣基因文件.csv"

# 3. 是否生成随机颜色
是否生成随机颜色 <- FALSE

# 4. 火山图参数
log2FC_name <- "logFC"
fdr_name <- "P.Value"
log2FC <- 1
fdr <- 0.05
colors <- "#39489f, #39bbec, #f9ed36, #f38466, #b81f25"
legend_title <- "P-value"
y_lab <- "-log10 P-value"
legend_position <- "UR"
output_dir <- "定制的美化图片"
pdf_width <- 8
pdf_height <- 6

# 5. 如需生成随机颜色，则执行
if (是否生成随机颜色) {
  random_colors <- replicate(100, randomColor())
  print(random_colors)
}

# 6. 读取差异分析结果文件
deg_data <- read.csv(file_path, row.names = 1)

# 7. 读取感兴趣基因文件
gene_data <- read.csv(genename_file, header = FALSE, stringsAsFactors = FALSE)
genes_to_label <- as.character(gene_data[, 1])

# 8. 添加调控分组信息
data <- add_regulate(
  deg_data,
  log2FC_name = log2FC_name,
  fdr_name = fdr_name,
  log2FC = log2FC,
  fdr = fdr
)

data$geneName <- rownames(data)

# 9. 处理自定义颜色
custom_colors <- colorRampPalette(unlist(strsplit(colors, ",\\s*")))(100)

# 10. 创建火山图
p32 <- gradual_volcano(
  data,
  x = "log2FoldChange",
  y = "padj",
  fills = custom_colors,
  colors = custom_colors,
  label = "row",
  log2FC_cut = log2FC,
  custom_label = genes_to_label,
  output = FALSE,
  legend_title = legend_title,
  y_lab = y_lab,
  legend_position = legend_position
)

# 11. 创建输出目录
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 12. 保存火山图为 PDF
pdf_path <- file.path(output_dir, "volcano_plot.pdf")
ggsave(pdf_path, plot = p32, width = pdf_width, height = pdf_height)

# 13. 输出结果
print(p32)
cat("火山图已保存为 PDF 文件：", pdf_path, "\n")
