library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(ggprism)
library(grid)

options(shiny.maxRequestSize = 1 * 1024^3)

# 1. 定义生成随机颜色的函数
randomColor <- function() {
  paste0("#", paste0(sample(c(0:9, letters[1:6]), 6, replace = TRUE), collapse = ""))
}

# 2. 设置输入文件路径
input_file <- "差异分析结果文件.csv"
gene_names_file <- "感兴趣基因文件.csv"

# 3. 是否生成随机颜色
是否生成随机颜色 <- FALSE

# 4. 随机颜色种子
random_seed <- 1234

# 5. 火山图参数设置
alpha_value <- 0.85
point_size <- 1.5
down_color <- "steelblue"
no_change_color <- "gray"
up_color <- "brown"
vline_position <- 1
x_limit <- 11
x_label <- "Log Fold Change"
hline_pvalue <- 0.05
pvalue_column <- "P.Value"
y_label <- "-Log10 p-value"
plot_title <- "DEG names"
pdf_width <- 10
pdf_height <- 8
pdf_file_name <- "GeneName_volcano.pdf"
segment_color <- "black"
label_size <- 3
box_padding <- 0.5
point_padding <- 0.8
legend_position <- "right"

# 6. 如需生成随机颜色，则执行
if (是否生成随机颜色) {
  set.seed(random_seed)
  random_colors <- replicate(100, randomColor())
  print(random_colors)
}

# 7. 读取数据
BRCA_Match_DEG <- read.csv(input_file, stringsAsFactors = FALSE)
gene_data <- read.csv(gene_names_file, header = FALSE, stringsAsFactors = FALSE)
gene_names <- as.character(gene_data[, 1])

# 8. 给感兴趣基因添加标签
BRCA_Match_DEG <- BRCA_Match_DEG %>%
  mutate(label = ifelse(X %in% gene_names, X, ""))

# 9. 动态选择 p-value 列
p_value_col <- pvalue_column

# 10. 绘制火山图
p <- ggplot(BRCA_Match_DEG, aes(x = logFC, y = -log10(.data[[p_value_col]]), colour = change)) +
  geom_point(alpha = alpha_value, size = point_size) +
  scale_color_manual(values = c(down_color, no_change_color, up_color)) +
  xlim(c(-x_limit, x_limit)) +
  geom_vline(xintercept = c(-vline_position, vline_position), lty = 4, col = "black", lwd = 0.8) +
  geom_hline(yintercept = -log10(hline_pvalue), lty = 4, col = "black", lwd = 0.8) +
  labs(x = x_label, y = y_label) +
  ggtitle(plot_title) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = legend_position,
    legend.title = element_blank()
  ) +
  geom_label_repel(
    aes(label = label),
    size = label_size,
    box.padding = unit(box_padding, "lines"),
    point.padding = unit(point_padding, "lines"),
    segment.color = segment_color,
    show.legend = FALSE,
    max.overlaps = 20000
  ) +
  theme_prism(border = TRUE)

# 11. 创建输出目录
if (!dir.exists("1.差异分析/指定基因的可视化")) {
  dir.create("1.差异分析/指定基因的可视化", recursive = TRUE)
}

# 12. 设置输出路径
output_path <- file.path("1.差异分析/指定基因的可视化", pdf_file_name)

# 13. 保存图表到 PDF
ggsave(output_path, plot = p, width = pdf_width, height = pdf_height, device = "pdf")

# 14. 输出参数摘要
summary_info <- list(
  Alpha = alpha_value,
  PointSize = point_size,
  Colors = c(down_color, no_change_color, up_color),
  FileName = input_file,
  PDFName = pdf_file_name,
  PlotTitle = plot_title,
  VlinePositions = c(-vline_position, vline_position),
  HlinePosition = -log10(hline_pvalue),
  PValueType = pvalue_column,
  XLimit = x_limit,
  GeneNames = gene_names,
  RandomSeed = random_seed
)

print(summary_info)
cat("数据处理完成。图表已保存为PDF。\n")