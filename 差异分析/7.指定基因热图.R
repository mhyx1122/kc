library(ggplot2)
library(pheatmap)
library(reshape2)
library(fs)

# 1. 随机生成颜色的函数
randomColor <- function() {
  paste0("#", paste0(sample(c(0:9, letters[1:6]), 6, replace = TRUE), collapse = ""))
}

# 2. 设置输入文件路径
file_path <- "基因表达CSV文件.csv"
group_file_path <- "分组文件.csv"
gene_file_path <- "感兴趣基因CSV文件.csv"

# 3. 是否生成随机颜色
是否生成随机颜色 <- FALSE

# 4. 热图参数设置
plot_width <- 6
plot_height <- 5

scale_method <- "row"
cluster_cols <- FALSE
cluster_rows <- TRUE
show_rownames <- TRUE
show_colnames <- TRUE

color_1 <- "red"
color_2 <- "blue"

color_3 <- "navy"
color_4 <- "white"
color_5 <- "firebrick"

fontsize <- 8
fontsize_row <- 6
fontsize_col <- 6
angle_col <- 45

# 5. 输出文件夹名称
output_folder <- "定制的美化图片"

# 6. 如需生成随机颜色，则执行
if (是否生成随机颜色) {
  random_colors <- replicate(100, randomColor())
  print(random_colors)
}

# 7. 检查输出文件夹是否存在，不存在则创建
if (!dir_exists(output_folder)) {
  dir_create(output_folder)
}

# 8. 读取数据
rt <- read.csv(file_path, header = TRUE, sep = ",", row.names = 1)
ann <- read.csv(group_file_path, header = TRUE, sep = ",", row.names = 1, check.names = FALSE)

# 9. 读取感兴趣基因文件
gene_data <- read.csv(gene_file_path, header = FALSE, stringsAsFactors = FALSE)
gene_names <- as.character(gene_data[, 1])

# 10. 提取感兴趣基因的数据
rt_subset <- rt[gene_names, , drop = FALSE]

# 11. 获取分组名称并设置颜色
unique_groups <- unique(ann$group)

cloor_use <- c(color_1, color_2)

# 12. 检查颜色数量是否足够
if (length(cloor_use) < length(unique_groups)) {
  stop("The number of colors in 'cloor_use' is less than the number of unique groups.")
}

# 13. 生成 annotation_colors
annotation_colors <- list(
  group = setNames(cloor_use[1:length(unique_groups)], unique_groups)
)

# 14. 设置输出文件路径
output_file_path <- file.path(output_folder, "heatmap.pdf")

# 15. 绘制并保存热图
pdf(output_file_path, width = plot_width, height = plot_height)
pheatmap(
  rt_subset,
  annotation_col = ann,
  annotation_colors = annotation_colors,
  cluster_cols = cluster_cols,
  cluster_rows = cluster_rows,
  color = colorRampPalette(c(color_3, color_4, color_5))(50),
  scale = scale_method,
  border_color = NA,
  fontsize = fontsize,
  fontsize_row = fontsize_row,
  fontsize_col = fontsize_col,
  angle_col = angle_col,
  show_rownames = show_rownames,
  show_colnames = show_colnames
)
dev.off()

# 16. 输出完成提示
cat("热图已生成并保存为 PDF 文件：", output_file_path, "\n")

