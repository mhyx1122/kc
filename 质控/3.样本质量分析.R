# 1. 加载R包

suppressPackageStartupMessages({
  library(ggplot2)
  library(dendextend)
  library(ComplexHeatmap)
})


# 2. 读取表达矩阵

# 第一列为基因名或特征名，其余列为样本
input_file <- "表达矩阵.csv"

data_cluster <- read.csv(input_file, header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

if (nrow(data_cluster) < 2 || ncol(data_cluster) < 2) {
  stop("表达矩阵至少需要包含2个特征和2个样本。")
}

if (anyDuplicated(colnames(data_cluster))) {
  stop("表达矩阵中存在重复的样本名称，请先处理重复列名。")
}

data_cluster[] <- lapply(data_cluster, function(x) suppressWarnings(as.numeric(x)))

if (anyNA(data_cluster)) {
  stop("表达矩阵中存在缺失值或无法转换为数值的内容，请检查输入文件。")
}

data10 <- as.matrix(data_cluster)
storage.mode(data10) <- "numeric"

if (any(!is.finite(data10))) {
  stop("表达矩阵中存在Inf或-Inf，请检查输入数据。")
}


# 3. 设置样本分组

# 分组名称和每组样本数必须一一对应
# 分组顺序按照表达矩阵中的样本列顺序依次设置
group_names <- c("T", "N")
group_counts <- c(13, 18)

# 三组示例
# group_names <- c("Control", "Treatment1", "Treatment2")
# group_counts <- c(10, 10, 10)

# 四组示例
# group_names <- c("A", "B", "C", "D")
# group_counts <- c(8, 8, 8, 8)

if (length(group_names) < 2) {
  stop("至少需要设置两个分组。")
}

if (length(group_names) != length(group_counts)) {
  stop("group_names和group_counts的长度不一致。")
}

if (anyNA(group_names) || any(trimws(group_names) == "") || anyDuplicated(group_names)) {
  stop("分组名称不能包含缺失值、空名称或重复名称。")
}

if (!is.numeric(group_counts) || anyNA(group_counts) || any(group_counts <= 0) || any(group_counts %% 1 != 0)) {
  stop("每组样本数量必须为大于0的整数。")
}

if (sum(group_counts) != ncol(data10)) {
  stop(
    "各组样本数之和为", sum(group_counts),
    "，但表达矩阵包含", ncol(data10),
    "个样本，请检查group_counts。"
  )
}

group <- factor(rep(group_names, times = group_counts), levels = group_names)
names(group) <- colnames(data10)


# 4. 创建结果文件夹

output_folder <- "样本分析结果"

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}


# 5. 设置层次聚类参数

distance_method <- "euclidean"
hclust_method <- "mcquitty"

hclust_pdf_width <- 10
hclust_pdf_height <- 8
dend_branch_lwd <- 2
dend_labels_cex <- 0.8

# 顺序为bottom、left、top、right
plot_margins <- c(5, 3, 2, 2)


# 6. 进行样本层次聚类

distance <- dist(t(data10), method = distance_method)
hc <- hclust(distance, method = hclust_method)

dend <- as.dendrogram(hc)
dend <- dendextend::set(dend, "branches_lwd", dend_branch_lwd)
dend <- dendextend::set(dend, "labels_cex", dend_labels_cex)

pdf(file.path(output_folder, "样本聚类树.pdf"), width = hclust_pdf_width, height = hclust_pdf_height)
par(mar = plot_margins)
plot(dend)
dev.off()

message("层次聚类树绘制完成。")


# 7. 设置样本相关性分析参数

correlation_method <- "spearman"

correlation_pdf_width <- 10
correlation_pdf_height <- 8

heatmap_font_size <- 0.5

# 该参数用于热图单元格的边框颜色
heatmap_fill_color <- "white"

# 是否显示热图边框
heatmap_border <- TRUE

# TRUE时进行聚类并显示聚类树，FALSE时不进行聚类
show_dendrogram <- TRUE


# 8. 计算样本相关性

dat_cormatrix <- cor(data10, method = correlation_method)

if (anyNA(dat_cormatrix)) {
  stop("样本相关性矩阵中出现NA，可能存在数值完全相同的样本。")
}

cell_border_color <- if (heatmap_border) heatmap_fill_color else NA


# 9. 绘制样本相关性热图

correlation_heatmap <- Heatmap(
  matrix = dat_cormatrix,
  name = correlation_method,
  border = heatmap_border,
  cluster_rows = show_dendrogram,
  cluster_columns = show_dendrogram,
  clustering_distance_rows = distance_method,
  clustering_method_rows = hclust_method,
  clustering_distance_columns = distance_method,
  clustering_method_columns = hclust_method,
  row_dend_side = "left",
  column_dend_side = "top",
  row_names_gp = grid::gpar(fontsize = 12),
  column_names_gp = grid::gpar(fontsize = 12),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid::grid.rect(x, y, width, height, gp = grid::gpar(fill = fill, col = cell_border_color))
    grid::grid.text(
      sprintf("%.2f", dat_cormatrix[i, j]),
      x,
      y,
      gp = grid::gpar(col = "black", fontface = "bold", cex = heatmap_font_size)
    )
  }
)

pdf(file.path(output_folder, "样本相关性.pdf"), width = correlation_pdf_width, height = correlation_pdf_height)
ComplexHeatmap::draw(correlation_heatmap)
dev.off()

message("相关性分析完成。")


# 10. 设置PCA分析参数

pca_pdf_width <- 6
pca_pdf_height <- 5

pca_title <- "PCA of Data"
pca_title_size <- 15
pca_axis_title_size <- 15

# 按照group_names的顺序设置分组颜色
color_range <- c("#FF5733", "#E69F00", "#3357FF", "#56B4E9", "#009E73", "#FF33FF")

if (length(color_range) < length(group_names)) {
  color_range <- grDevices::hcl.colors(length(group_names), palette = "Dark 3")
} else {
  color_range <- color_range[seq_along(group_names)]
}

group_colors <- setNames(color_range, group_names)


# 11. 过滤零方差特征

df_t <- t(data10)
var_nonzero_pre <- apply(df_t, 2, var) > 0
df_filtered_pre <- df_t[, var_nonzero_pre, drop = FALSE]

if (ncol(df_filtered_pre) < 2) {
  stop("过滤零方差特征后，剩余特征少于2个，无法进行PCA分析。")
}


# 12. 进行PCA分析

pca_res_pre <- prcomp(df_filtered_pre, center = TRUE, scale. = TRUE)

if (ncol(pca_res_pre$x) < 2) {
  stop("当前数据无法提取前两个主成分。")
}

explained_var <- summary(pca_res_pre)$importance[2, 1:2] * 100

pca_scores_pre <- as.data.frame(pca_res_pre$x[, 1:2, drop = FALSE])
pca_scores_pre$sample <- rownames(pca_scores_pre)
pca_scores_pre$group <- group


# 13. 绘制PCA图

pca_plot_pre <- ggplot(pca_scores_pre, aes(x = PC1, y = PC2, color = group, fill = group)) +
  geom_point(size = 4, alpha = 0.9, shape = 16) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  labs(
    title = pca_title,
    x = sprintf("PC1 (%.1f%%)", explained_var[1]),
    y = sprintf("PC2 (%.1f%%)", explained_var[2]),
    color = "Group",
    fill = "Group"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, size = pca_title_size, face = "bold"),
    axis.title = element_text(size = pca_axis_title_size, face = "bold", color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    axis.line = element_line(linewidth = 0.8, color = "black"),
    axis.ticks = element_line(linewidth = 0.8, color = "black"),
    axis.ticks.length = grid::unit(0.2, "cm"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.position = "right",
    legend.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(12, 15, 10, 10)
  )


# 14. 为样本数不少于3的分组添加95%置信椭圆

group_sample_counts <- table(pca_scores_pre$group)
ellipse_groups <- names(group_sample_counts[group_sample_counts >= 3])

if (length(ellipse_groups) > 0) {
  ellipse_data <- pca_scores_pre[pca_scores_pre$group %in% ellipse_groups, , drop = FALSE]
  
  pca_plot_pre <- pca_plot_pre +
    stat_ellipse(
      data = ellipse_data,
      geom = "polygon",
      level = 0.95,
      alpha = 0.15,
      linewidth = 0.8,
      show.legend = FALSE
    )
}


# 15. 保存PCA图

ggsave(
  filename = file.path(output_folder, "PCA分析图.pdf"),
  plot = pca_plot_pre,
  width = pca_pdf_width,
  height = pca_pdf_height
)

message("PCA分析完成。")


# 16. 输出分析流程说明

cat(
  "主要流程，写作参考：\n",
  "1. 层次聚类分析：计算样本间距离矩阵，并基于指定的聚类方法进行层次聚类。\n",
  "2. 相关性分析：计算样本间相关系数，并使用热图展示样本相关性结构。\n",
  "3. 主成分分析：过滤零方差特征后进行PCA，并展示不同分组样本在前两个主成分中的分布。\n",
  sep = ""
)

message("样本质量控制与可视化分析完成，结果已保存至：", output_folder)