library(dplyr)
library(tidyr)
library(ggplot2)
library(fs)

# 1. 设置输入文件路径
ssgsea_file <- "基因表达.csv"
group_file <- "分组匹配信息.csv"
gene_file <- "感兴趣基因文件.csv"

# 2. 设置颜色参数
color1 <- "#e5451d"
color2 <- "#9084bd"

# 3. 设置统计检验方法
# 可选："t.test" 或 "wilcox.test"
test_method <- "t.test"

# 4. 是否添加散点图层
plot_points <- TRUE

# 5. 图形透明度参数
violin_alpha <- 0.6
scatter_alpha <- 0.6

# 6. PDF 输出大小
pdf_width <- 10
pdf_height <- 15

# 7. p 值显示方式
# 可选："numeric" 或 "stars"
p_display <- "numeric"

# 8. 输出文件夹名称
output_folder <- "定制的美化图片"

# 9. 检查输出文件夹是否存在，不存在则创建
if (!dir_exists(output_folder)) {
  dir_create(output_folder)
}

# 10. 读取表达矩阵
data <- read.csv(
  ssgsea_file,
  header = TRUE,
  row.names = 1,
  stringsAsFactors = FALSE
)

# 11. 读取感兴趣基因文件
gene_data <- read.csv(gene_file, header = FALSE, stringsAsFactors = FALSE)
gene_names <- as.character(gene_data[, 1])

# 12. 提取感兴趣基因
data <- data[gene_names, , drop = FALSE]

# 13. 判断是否需要 log2 转换
expMatrix <- data
ex <- expMatrix
qx <- as.numeric(quantile(ex, c(0.00, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))
LogC <- (qx[6] > 100) || (qx[6] - qx[1] > 50)

if (LogC) {
  expMatrix <- log2(ex + 1)
  print("log2 transform finished")
} else {
  print("log2 transform not needed")
}

# 14. 转置表达矩阵
data <- t(expMatrix)

# 15. 如存在 P_value 行则去除
if ("P_value" %in% rownames(data)) {
  data <- data[rownames(data) != "P_value", ]
}

data <- as.data.frame(data)

# 16. 读取分组匹配信息
new_dbssgsea <- read.csv(
  group_file,
  header = TRUE,
  row.names = 1,
  stringsAsFactors = FALSE
)

data$group <- new_dbssgsea[rownames(data), "group"]

# 17. 转成长格式
cols_to_include <- setdiff(colnames(data), "group")
data_long <- pivot_longer(
  data,
  cols = all_of(cols_to_include),
  names_to = "Metric",
  values_to = "Value"
)

# 18. 计算每个基因的 p 值
if (test_method == "t.test") {
  pvalues <- data_long %>%
    group_by(Metric) %>%
    summarise(p = t.test(Value ~ group)$p.value, .groups = "drop")
} else if (test_method == "wilcox.test") {
  pvalues <- data_long %>%
    group_by(Metric) %>%
    summarise(p = wilcox.test(Value ~ group)$p.value, .groups = "drop")
}

# 19. 设置 p 值显示标签
if (p_display == "numeric") {
  pvalues <- pvalues %>%
    mutate(label = paste0("p = ", signif(p, digits = 3)))
} else {
  pvalues <- pvalues %>%
    mutate(label = case_when(
      p < 0.001 ~ "***",
      p < 0.01 ~ "**",
      p < 0.05 ~ "*",
      TRUE ~ "ns"
    ))
}

# 20. 构造 facet 标签
new_labels <- setNames(paste0(pvalues$Metric, " (", pvalues$label, ")"), pvalues$Metric)

# 21. 设置颜色
custom_colors <- c(color1, color2)

# 22. 绘制图形
p1 <- ggplot(data_long, aes(x = group, y = Value, fill = group)) +
  geom_violin(trim = TRUE, alpha = violin_alpha, color = NA) +
  {if (plot_points) geom_jitter(
    shape = 16,
    position = position_jitter(0.2),
    aes(color = group),
    alpha = scatter_alpha,
    size = 2
  ) else NULL} +
  geom_boxplot(width = 0.5, color = "black", fill = NA, outlier.shape = NA) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  facet_wrap(~ Metric, ncol = 3, scales = "free_y", labeller = as_labeller(new_labels)) +
  labs(title = " ", x = NULL, y = " ") +
  theme_minimal(base_size = 15) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "grey", color = "black", linewidth = 0.5),
    strip.text = element_text(size = 10),
    axis.ticks.y = element_line(color = "black")
  )

# 23. 保存 PDF 文件
output_file_path <- file.path(output_folder, "gene_violin.pdf")
ggsave(
  output_file_path,
  plot = p1,
  width = pdf_width,
  height = pdf_height,
  dpi = 300
)
