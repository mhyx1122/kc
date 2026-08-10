# 1. 加载必要 R 包
library(ConsensusClusterPlus)
library(pheatmap)
library(ggplot2)
library(Rtsne)
library(umap)
library(randomcoloR)

# 2. 设置输入文件
cluster_input_file <- "exp_data.csv"

# 3. 设置一致性聚类输出文件夹
output_folder <- "ConsensusCluster"

# 4. 设置 ConsensusClusterPlus 参数
maxK <- 9
reps <- 100
pItem <- 0.8
pFeature <- 1
clusterAlg <- "km"
distance <- "euclidean"
seed <- 123456

# 5. 设置一致性聚类热图渐变色
tmyPal_input <- c("white", "blue")

# 6. 设置指定 K 的一致性矩阵热图参数
kCluster <- 4

clusterColors <- c(
  "#E64B35FF",
  "#4DBBD5FF",
  "#00A087FF",
  "#3C5488FF",
  "#F39B7FFF",
  "#8491B4FF",
  "#91D1C2FF",
  "#DC0000FF",
  "#7E6148FF"
)

clustering_distance_cols <- "correlation"
clustering_method <- "average"
heatmap_width <- 8
heatmap_height <- 8

start_color <- "white"
end_color <- "steelblue"

# 7. 设置 PCA / t-SNE / UMAP 参数
k_numb <- 4

color_range <- c(
  "#E64B35FF",
  "#4DBBD5FF",
  "#00A087FF",
  "#3C5488FF",
  "#F39B7FFF",
  "#8491B4FF",
  "#91D1C2FF",
  "#DC0000FF",
  "#7E6148FF"
)

point_size <- 3
point_alpha <- 0.8

pca_pdf_width <- 8
pca_pdf_height <- 8
pca_title <- "PCA of Data"
pca_title_size <- 15
pca_axis_title_size <- 15

tsne_dims <- 2
tsne_perplexity <- 30
tsne_max_iter <- 500

umap_n_neighbors <- 15
umap_min_dist <- 0.1

# 8. 创建输出文件夹
if (!dir.exists(output_folder)) {
  dir.create(output_folder)
}

# 9. 读取聚类输入矩阵
df <- read.csv(
  cluster_input_file,
  row.names = 1,
  check.names = FALSE
)

df_matrix <- as.matrix(df)
storage.mode(df_matrix) <- "numeric"

if (anyNA(df_matrix)) {
  stop("输入矩阵中存在 NA。请检查是否有非数字字符、空值或异常符号。")
}

# 10. 提示输入矩阵方向
message("当前代码按：基因/特征在行，样本在列 进行一致性聚类。")

# 11. 生成一致性聚类热图调色板
tmyPal <- colorRampPalette(tmyPal_input)(100)

# 12. 运行 ConsensusClusterPlus
ccres <- ConsensusClusterPlus(
  df_matrix,
  maxK = maxK,
  reps = reps,
  pItem = pItem,
  pFeature = pFeature,
  tmyPal = tmyPal,
  title = output_folder,
  clusterAlg = clusterAlg,
  distance = distance,
  seed = seed,
  plot = "pdf"
)

# 13. 检查指定 K 是否可用
if (!kCluster %in% seq_along(ccres)) {
  stop("kCluster 不在 ccres 可用范围内。请检查 kCluster 是否小于等于 maxK。")
}

# 14. 提取指定 K 的一致性矩阵
ConsensusMatrix <- data.frame(
  ccres[[kCluster]][["consensusMatrix"]]
)

ConsensusMatrix <- ConsensusMatrix[
  ccres[[kCluster]]$consensusTree$order,
  ccres[[kCluster]]$consensusTree$order
]

# 15. 构建一致性矩阵热图注释
annCol <- data.frame(
  ccres = paste0(
    "Cluster",
    ccres[[kCluster]][["consensusClass"]][
      ccres[[kCluster]]$consensusTree$order
    ]
  ),
  row.names = colnames(ConsensusMatrix)
)

sorted_order <- order(annCol$ccres)
annCol <- annCol[sorted_order, , drop = FALSE]

clusterNames <- unique(annCol$ccres)

if (length(clusterColors) < length(clusterNames)) {
  stop("clusterColors 的颜色数量少于聚类数量，请增加颜色。")
}

annColors <- list(
  ccres = setNames(
    clusterColors[seq_along(clusterNames)],
    clusterNames
  )
)

# 16. 绘制并保存指定 K 的一致性矩阵热图
heatmapFile <- paste0(
  output_folder,
  "/K=",
  kCluster,
  "_Heatmap.pdf"
)

pdf(
  heatmapFile,
  width = heatmap_width,
  height = heatmap_height
)

Heatmap <- pheatmap(
  ConsensusMatrix,
  color = colorRampPalette(c(start_color, end_color))(100),
  clustering_distance_cols = clustering_distance_cols,
  clustering_method = clustering_method,
  border_color = NA,
  annotation_col = annCol,
  annotation_colors = annColors,
  show_colnames = FALSE,
  show_rownames = FALSE
)

print(Heatmap)

dev.off()

# 17. 检查 PCA / t-SNE / UMAP 使用的 K 是否可用
if (!k_numb %in% seq_along(ccres)) {
  stop("k_numb 不在 ccres 可用范围内。请检查 k_numb 是否小于等于 maxK。")
}

# 18. 提取指定 K 的聚类结果
Cluster <- ccres[[k_numb]]

AA <- Cluster$consensusClass
original_names <- names(AA)

cluster_output <- data.frame(
  Sample = names(AA),
  Cluster = paste0("Cluster", AA)
)

cluster_file_name <- paste0(
  output_folder,
  "/cluster_output(K=",
  k_numb,
  ").csv"
)

write.csv(
  cluster_output,
  file = cluster_file_name,
  row.names = FALSE
)

AA <- paste0("Cluster", AA)
names(AA) <- original_names
AA <- as.factor(AA)

# 19. 转置表达矩阵，使样本为行、基因为列
df_t <- t(df_matrix)

if (!all(names(AA) %in% rownames(df_t))) {
  stop("聚类结果中的样本名和转置后的输入矩阵行名不一致。请检查输入矩阵方向。")
}

df_t <- df_t[names(AA), , drop = FALSE]

# 20. 过滤方差为 0 的特征
var_nonzero_pre <- apply(df_t, 2, var) != 0

df_filtered_pre <- df_t[, var_nonzero_pre, drop = FALSE]

if (ncol(df_filtered_pre) < 2) {
  stop("过滤方差为 0 的特征后，剩余特征数少于 2，无法进行 PCA / t-SNE / UMAP。")
}

# 21. 根据聚类数量准备颜色
group_levels <- levels(AA)

num_colors_needed <- length(group_levels) - length(color_range)

if (num_colors_needed > 0) {
  additional_colors <- randomcoloR::distinctColorPalette(num_colors_needed)
  colors <- c(color_range, additional_colors)
} else {
  colors <- color_range[seq_along(group_levels)]
}

names(colors) <- group_levels

# 22. 进行 PCA 分析
pca_res_pre <- prcomp(
  df_filtered_pre,
  scale. = TRUE
)

explained_var <- summary(pca_res_pre)$importance[2, ] * 100

pca_scores_pre <- as.data.frame(
  pca_res_pre$x[, 1:2, drop = FALSE]
)

pca_scores_pre$group <- AA

pca_file_name <- paste0(
  output_folder,
  "/PCA_scores(K=",
  k_numb,
  ").csv"
)

write.csv(
  pca_scores_pre,
  file = pca_file_name,
  row.names = TRUE
)

# 23. 绘制 PCA 图
pca_plot_pre <- ggplot(
  pca_scores_pre,
  aes(x = PC1, y = PC2, color = group)
) +
  geom_point(
    size = point_size,
    alpha = point_alpha
  ) +
  stat_ellipse(
    aes(fill = group, color = group),
    geom = "polygon",
    level = 0.95,
    alpha = 0.2
  ) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  labs(
    title = pca_title,
    x = sprintf("Principal Component 1 (%.1f%%)", explained_var[1]),
    y = sprintf("Principal Component 2 (%.1f%%)", explained_var[2]),
    color = "Group",
    fill = "Group"
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      size = pca_title_size
    ),
    axis.title = element_text(
      size = pca_axis_title_size
    )
  )

pca_pdf_file <- paste0(
  output_folder,
  "/PCA_analysis_plot(K=",
  k_numb,
  ").pdf"
)

pdf(
  pca_pdf_file,
  width = pca_pdf_width,
  height = pca_pdf_height
)

print(pca_plot_pre)

dev.off()

# 24. 检查 t-SNE 参数
if (tsne_dims != 2) {
  stop("当前绘图代码只支持二维 t-SNE，请设置 tsne_dims <- 2。")
}

if (tsne_perplexity >= (nrow(df_filtered_pre) - 1) / 3) {
  stop("tsne_perplexity 过大。通常需要小于 (样本数 - 1) / 3。")
}

# 25. 运行 t-SNE
set.seed(seed)

tsne_res <- Rtsne::Rtsne(
  df_filtered_pre,
  dims = tsne_dims,
  perplexity = tsne_perplexity,
  verbose = TRUE,
  max_iter = tsne_max_iter
)

tsne_scores_pre <- as.data.frame(tsne_res$Y)

colnames(tsne_scores_pre) <- c("tSNE1", "tSNE2")

tsne_scores_pre$group <- AA

tsne_file_name <- paste0(
  output_folder,
  "/tSNE_scores(K=",
  k_numb,
  ").csv"
)

write.csv(
  tsne_scores_pre,
  file = tsne_file_name,
  row.names = TRUE
)

# 26. 绘制 t-SNE 图
tsne_plot_pre <- ggplot(
  tsne_scores_pre,
  aes(x = tSNE1, y = tSNE2, color = group)
) +
  geom_point(
    size = point_size,
    alpha = point_alpha
  ) +
  stat_ellipse(
    aes(fill = group, color = group),
    geom = "polygon",
    level = 0.95,
    alpha = 0.2
  ) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  labs(
    title = "t-SNE of Data",
    x = "tSNE 1",
    y = "tSNE 2",
    color = "Group",
    fill = "Group"
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      size = pca_title_size
    ),
    axis.title = element_text(
      size = pca_axis_title_size
    )
  )

tsne_pdf_file <- paste0(
  output_folder,
  "/tSNE_analysis_plot(K=",
  k_numb,
  ").pdf"
)

pdf(
  tsne_pdf_file,
  width = pca_pdf_width,
  height = pca_pdf_height
)

print(tsne_plot_pre)

dev.off()

# 27. 运行 UMAP
set.seed(seed)

umap_config <- umap::umap.defaults

umap_config$n_neighbors <- umap_n_neighbors
umap_config$min_dist <- umap_min_dist

umap_res <- umap::umap(
  df_filtered_pre,
  config = umap_config
)

umap_scores_pre <- as.data.frame(umap_res$layout)

colnames(umap_scores_pre) <- c("UMAP1", "UMAP2")

umap_scores_pre$group <- AA

umap_file_name <- paste0(
  output_folder,
  "/UMAP_scores(K=",
  k_numb,
  ").csv"
)

write.csv(
  umap_scores_pre,
  file = umap_file_name,
  row.names = TRUE
)

# 28. 绘制 UMAP 图
umap_plot_pre <- ggplot(
  umap_scores_pre,
  aes(x = UMAP1, y = UMAP2, color = group)
) +
  geom_point(
    size = point_size,
    alpha = point_alpha
  ) +
  stat_ellipse(
    aes(fill = group, color = group),
    geom = "polygon",
    level = 0.95,
    alpha = 0.2
  ) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  labs(
    title = "UMAP of Data",
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Group",
    fill = "Group"
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      size = pca_title_size
    ),
    axis.title = element_text(
      size = pca_axis_title_size
    )
  )

umap_pdf_file <- paste0(
  output_folder,
  "/UMAP_analysis_plot(K=",
  k_numb,
  ").pdf"
)

pdf(
  umap_pdf_file,
  width = pca_pdf_width,
  height = pca_pdf_height
)

print(umap_plot_pre)

dev.off()

# 29. 输出完成提示
message("一致性聚类分析完成。")
message("ConsensusClusterPlus 结果已保存到：", output_folder)
message("指定 K 的一致性矩阵热图已保存：", heatmapFile)
message("聚类分组文件已保存：", cluster_file_name)
message("PCA 图已保存：", pca_pdf_file)
message("t-SNE 图已保存：", tsne_pdf_file)
message("UMAP 图已保存：", umap_pdf_file)