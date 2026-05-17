# 1. 随机生成颜色的函数
randomColor <- function() {
  paste0("#", paste0(sample(c(0:9, letters[1:6]), 6, replace = TRUE), collapse = ""))
}

# 2. 设置输入文件路径
bulk_sample_file <- "选择的输入文件.csv"

# 3. 分组名称和样本数量输入
group1_name <- "T"
group1_count <- 13

group2_name <- "N"
group2_count <- 18

# 4. 数据处理参数（样本聚类时）
distance_method <- "euclidean"
hclust_method <- "mcquitty"

# 5. 数据处理参数（样本相关性分析时）
correlation_method <- "spearman"

# 6. 作图参数（聚类图）
hclust_pdf_width <- 10
hclust_pdf_height <- 8
dend_branch_lwd <- 2
dend_labels_cex <- 0.8
plot_margins <- c(5, 3, 2, 2)

# 7. 作图参数（样本相关性图）
correlation_pdf_width <- 10
correlation_pdf_height <- 8
heatmap_font_size <- 0.5
heatmap_fill_color <- "white"
heatmap_border <- TRUE
show_dendrogram <- TRUE

# 8. 作图参数（PCA图）
color_range <- c("#FF5733", "#E69F00", "#3357FF", "#56B4E9", "#009E73", "#FF33FF")
pca_pdf_width <- 6
pca_pdf_height <- 5
pca_title <- "PCA of Data"
pca_title_size <- 15
pca_axis_title_size <- 15

# 9. 读取输入文件
data_cluster <- read.csv(bulk_sample_file, header = TRUE, row.names = 1, check.names = FALSE)

# 10. 生成分组信息
group <- factor(c(
  rep(group1_name, group1_count),
  rep(group2_name, group2_count)
))

# 11. 计算距离
distance <- dist(t(data_cluster), method = distance_method)

# 12. 执行层次聚类
hc <- hclust(distance, method = hclust_method)
plot(hc)
dend <- as.dendrogram(hc)
dend <- set(dend, "branches_lwd", dend_branch_lwd)  # 设置所有树枝的默认线宽
dend <- set(dend, "labels_cex", dend_labels_cex)

# 13. 绘制聚类树，设置PDF设备，并指定长宽
pdf('样本分析结果/样本聚类树.pdf', width = hclust_pdf_width, height = hclust_pdf_height)
# 调整边距
par(mar=plot_margins)  # 设置边距
# 绘制聚类树
plot(dend)
# 关闭PDF设备
dev.off()
print("层次聚类树绘制完成")

# 14. 执行相关性分析
dat_cormatrix <- cor(data_cluster, method = correlation_method)

# 15. 绘制相关性热图
library(ComplexHeatmap)
p <- Heatmap(
  matrix = dat_cormatrix,
  name = correlation_method,
  border = heatmap_border,
  cluster_rows = show_dendrogram,  # 是否显示行聚类树
  cluster_columns = show_dendrogram,  # 是否显示列聚类树
  clustering_distance_rows = distance_method,  # 行聚类的距离计算方法
  clustering_method_rows = hclust_method,  # 行聚类的方法
  clustering_distance_columns = distance_method,  # 列聚类的距离计算方法
  clustering_method_columns = hclust_method,  # 列聚类的方法
  row_dend_side = 'left',  # 聚类树显示在左边
  column_dend_side = 'top',  # 聚类树显示在顶部
  row_names_gp = gpar(fontsize = 12),  # 行名字体大小
  column_names_gp = gpar(fontsize = 12),  # 列名字体大小
  cell_fun = function(j, i, x, y, width, height, fill) {  # 自定义单元格绘制
    grid.rect(x, y, width, height, gp = gpar(fill = fill, col = heatmap_fill_color))
    grid.text(sprintf('%.2f', dat_cormatrix[i, j]), x, y, gp = gpar(col = 'black', fontface = 'bold', cex = heatmap_font_size))
  }
)

# 16. 绘制相关性热图，设置PDF设备，并指定长宽
pdf('样本分析结果/样本相关性.pdf', width = correlation_pdf_width, height = correlation_pdf_height)
# 绘制热图
plot(p)
# 关闭PDF设备
dev.off()
print("相关性分析完成")

# 17. PCA分析
df_t <- t(data_cluster)
# 过滤掉方差为0的列
var_nonzero_pre <- apply(df_t, 2, var) != 0
df_filtered_pre <- df_t[, var_nonzero_pre]

# 进行PCA分析
pca_res_pre <- prcomp(df_filtered_pre, scale. = TRUE)
# 获取解释度信息
explained_var <- summary(pca_res_pre)$importance[2, 1:2] * 100
# 准备绘图数据（提取前两个主成分）
pca_scores_pre <- as.data.frame(pca_res_pre$x[, 1:2])
pca_scores_pre$group <- group  # group是事先定义的

# 指定颜色范围
library(randomcoloR)
group_levels <- levels(group)
group_colors <- setNames(distinctColorPalette(length(group_levels)), group_levels)

# 18. 使用ggplot2绘制PCA图
library(ggplot2)
pca_plot_pre <- ggplot(pca_scores_pre, aes(x = PC1, y = PC2, color = group, fill = group)) +
  geom_point(size = 4, alpha = 0.9, shape = 16) +
  stat_ellipse(
    geom = "polygon",
    level = 0.95,
    alpha = 0.15,
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  labs(
    title = "PCA Analysis",
    x = sprintf("PC1 (%.1f%%)", explained_var[1]),
    y = sprintf("PC2 (%.1f%%)", explained_var[2]),
    color = "Group",
    fill = "Group"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14, face = "bold", color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    axis.line = element_line(linewidth = 0.8, color = "black"),
    axis.ticks = element_line(linewidth = 0.8, color = "black"),
    axis.ticks.length = unit(0.2, "cm"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.position = "right",
    legend.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(12, 15, 10, 10)
  )
pca_plot_pre

# 19. 保存PCA图形到PDF
pdf('样本分析结果/PCA分析图.pdf', width = pca_pdf_width, height = pca_pdf_height)
plot(pca_plot_pre)
dev.off()  # 关闭PDF设备
print("PCA分析完成")