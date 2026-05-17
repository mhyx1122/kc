library(ggplot2)
library(sva)
library(patchwork)

# 1. 随机生成颜色的函数
randomColor <- function() {
  paste0("#", paste0(sample(c(0:9, letters[1:6]), 6, replace = TRUE), collapse = ""))
}

# 2. 设置数据集数量
data_count <- 2

# 3. 设置每个批次数据文件路径
batch_files <- c(
  "第1个数据文件.csv",
  "第2个数据文件.csv"
)

# 4. 设置每个批次名称
dataset_names <- c("Dataset1", "Dataset2")

# 5. 设置每个批次样本个数
dataset_counts <- c(100, 100)

# 6. 设置分组信息文件路径
bulkQC_group_file <- "分组信息文件.csv"

# 7. 设置随机抽取样本数量进行箱线图可视化
bulkQC_boxplot <- 100

# 8. 设置箱线图参数
boxplot_width <- 15
boxplot_height <- 8
boxplot_color <- "lightblue"
boxplot_margins <- c(18, 4, 4, 2)

# 9. 设置 PCA 图参数
pca_plot_width <- 7
pca_plot_height <- 7
pca_colors <- c("red", "#E69F00", "#56B4E9", "#009E73")
group_colors <- c("red", "#E69F00", "green", "yellow")
pca_point_size <- 1
pca_alpha_value <- 0.2
pca_title_font_size <- 16
pca_axis_label_font_size <- 14

# 10. 如需生成 100 个随机颜色，可运行下面这行
generated_random_colors <- replicate(100, randomColor())
print(generated_random_colors)

# 11. 检查输入参数长度是否一致
if (length(batch_files) != data_count) {
  stop("batch_files 的数量与 data_count 不一致。")
}
if (length(dataset_names) != data_count) {
  stop("dataset_names 的数量与 data_count 不一致。")
}
if (length(dataset_counts) != data_count) {
  stop("dataset_counts 的数量与 data_count 不一致。")
}

# 12. 生成批次信息
bulkQC_batch <- factor(
  unlist(mapply(function(nm, ct) rep(nm, ct), dataset_names, dataset_counts, SIMPLIFY = FALSE)),
  levels = dataset_names
)

cat("生成的批次信息：\n")
for (i in seq_len(data_count)) {
  cat("数据集", i, ":", dataset_names[i], "- 样本数:", dataset_counts[i], "\n")
}

# 13. 读取多个数据集
all_data <- list()

for (i in seq_len(data_count)) {
  data <- read.csv(batch_files[i], stringsAsFactors = FALSE, check.names = FALSE)
  data <- data[!duplicated(data[, 1]), ]
  rownames(data) <- data[, 1]
  data <- data[, -1, drop = FALSE]
  all_data[[i]] <- data
}

# 14. 逐个合并数据集
merged_data <- all_data[[1]]

if (data_count >= 2) {
  for (i in 2:data_count) {
    merged_data <- merge(merged_data, all_data[[i]], by = "row.names", all = TRUE)
    rownames(merged_data) <- merged_data$Row.names
    merged_data <- merged_data[, -1, drop = FALSE]
  }
}

# 15. 读取分组信息
group_info <- read.csv(bulkQC_group_file, stringsAsFactors = FALSE, check.names = FALSE)

# 16. 设置 PCA 主题参数
pca_theme_plot_title <- element_text(hjust = 0.5, size = pca_title_font_size)
pca_theme_axis_text <- element_text(size = pca_axis_label_font_size)

# 17. 如果结果文件夹不存在，就创建
if (!dir.exists("合并数据去批次结果")) {
  dir.create("合并数据去批次结果")
}

# 18. 去除缺失值并保存合并后的数据
cleaned_data <- na.omit(merged_data)
write.csv(cleaned_data, "合并数据去批次结果/1.合并后的数据.csv", row.names = TRUE)
cat("3. 合并成功，文件已保存...\n")

# 19. 检查数据是否需要 log2 转换
data_type <- cleaned_data
ex <- data_type
qx <- as.numeric(quantile(ex, c(0.00, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))
LogC <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)

cat("4. 检查数据是否需要 log2 转换...\n")

if (LogC) {
  data_type <- log2(ex + 1)
  print("log2 transform finished")
} else {
  print("log2 transform not needed")
}

df <- data_type

# 20. 准备分组信息
group_info[, 1] <- gsub("-", ".", group_info[, 1])
group_info_reordered <- group_info[match(colnames(df), group_info$ID), ]
cat("5. 分组信息读取成功...\n")

# 21. 绘制去批次前箱线图
cat("6. 绘制箱线图看批次效应...\n")

bulkQC_boxplot_use <- min(bulkQC_boxplot, ncol(df))
selected_columns <- sort(sample(ncol(df), bulkQC_boxplot_use))
df_sampled <- df[, selected_columns, drop = FALSE]

pdf("合并数据去批次结果/2.去除批次前的箱线图.pdf", width = boxplot_width, height = boxplot_height)
par(mar = boxplot_margins)
boxplot(df_sampled, col = boxplot_color, las = 2, main = "Before Correction")
dev.off()

cat("7. 去除批次前的箱线图已绘制完毕并保存...\n")

# 22. 绘制去批次前 PCA 图（按批次给色）
df_t <- t(df)
var_nonzero_pre <- apply(df_t, 2, var) != 0
df_filtered_pre <- df_t[, var_nonzero_pre, drop = FALSE]

pca_res_pre <- prcomp(df_filtered_pre, scale. = TRUE)
explained_variance <- round(100 * pca_res_pre$sdev^2 / sum(pca_res_pre$sdev^2), 2)

pca_scores_pre <- as.data.frame(pca_res_pre$x[, 1:2])
pca_scores_pre$Batch <- bulkQC_batch

cat("8. 绘制 PCA 图的数据已准备...\n")

pca_colors_use <- rep_len(pca_colors, length(levels(bulkQC_batch)))

pca_plot_pre <- ggplot(pca_scores_pre, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point(size = pca_point_size) +
  stat_ellipse(aes(fill = Batch, color = Batch), geom = "polygon", level = 0.95, alpha = pca_alpha_value) +
  scale_color_manual(values = pca_colors_use) +
  scale_fill_manual(values = pca_colors_use) +
  theme_minimal() +
  labs(
    title = "PCA of Original Data",
    x = paste0("Principal Component 1 (", explained_variance[1], "%)"),
    y = paste0("Principal Component 2 (", explained_variance[2], "%)"),
    color = "Batch",
    fill = "Batch"
  ) +
  theme(
    plot.title = pca_theme_plot_title,
    axis.text = pca_theme_axis_text
  )

pdf("合并数据去批次结果/3.去除批次前的PCA图（按批次给色）.pdf", width = pca_plot_width, height = pca_plot_height)
print(pca_plot_pre)
dev.off()

cat("9. 去除批次前的 PCA 图（按批次给色）已绘制完毕并保存...\n")

# 23. 绘制去批次前 PCA 图（按分组给色）
group_labels <- factor(group_info_reordered$group)
group_colors_use <- rep_len(group_colors, length(levels(group_labels)))

pca_scores_pre$Group <- group_labels

pca_plot_pre_group <- ggplot(pca_scores_pre, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = pca_point_size) +
  stat_ellipse(aes(fill = Group, color = Group), geom = "polygon", level = 0.95, alpha = pca_alpha_value) +
  scale_color_manual(values = group_colors_use) +
  scale_fill_manual(values = group_colors_use) +
  theme_minimal() +
  labs(
    title = "PCA of Original Data by Group",
    x = paste0("Principal Component 1 (", explained_variance[1], "%)"),
    y = paste0("Principal Component 2 (", explained_variance[2], "%)"),
    color = "Group",
    fill = "Group"
  ) +
  theme(
    plot.title = pca_theme_plot_title,
    axis.text = pca_theme_axis_text
  )

pdf("合并数据去批次结果/4.去除批次前的PCA图（按分组给色）.pdf", width = pca_plot_width, height = pca_plot_height)
print(pca_plot_pre_group)
dev.off()

cat("10. 去除批次前的 PCA 图（按分组给色）已绘制完毕并保存...\n")

# 24. 开始去除批次效应
mod <- model.matrix(~ factor(group_info_reordered$group))
df_corrected <- ComBat(dat = as.matrix(df), batch = bulkQC_batch, mod = mod)

# 25. 绘制去批次后箱线图
selected_columns <- sort(sample(ncol(df_corrected), min(bulkQC_boxplot, ncol(df_corrected))))
df_sampled <- df_corrected[, selected_columns, drop = FALSE]

pdf("合并数据去批次结果/5.去除批次后的箱线图.pdf", width = boxplot_width, height = boxplot_height)
par(mar = boxplot_margins)
boxplot(df_sampled, col = boxplot_color, las = 2, main = "After Correction")
dev.off()

cat("11. 去除批次后的箱线图已绘制完毕并保存...\n")

# 26. 绘制去批次后 PCA 图（按批次给色）
df_corrected_t <- t(df_corrected)
var_nonzero <- apply(df_corrected_t, 2, var) != 0
df_filtered <- df_corrected_t[, var_nonzero, drop = FALSE]

pca_res <- prcomp(df_filtered, scale. = TRUE)
explained_variance_corrected <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 2)

pca_scores <- as.data.frame(pca_res$x[, 1:2])
pca_scores$Batch <- bulkQC_batch

pca_plot_aft <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point(size = pca_point_size) +
  stat_ellipse(aes(fill = Batch, color = Batch), geom = "polygon", level = 0.95, alpha = pca_alpha_value) +
  scale_color_manual(values = pca_colors_use) +
  scale_fill_manual(values = pca_colors_use) +
  theme_minimal() +
  labs(
    title = "PCA of Batch-Corrected Data",
    x = paste0("Principal Component 1 (", explained_variance_corrected[1], "%)"),
    y = paste0("Principal Component 2 (", explained_variance_corrected[2], "%)"),
    color = "Batch",
    fill = "Batch"
  ) +
  theme(
    plot.title = pca_theme_plot_title,
    axis.text = pca_theme_axis_text
  )

combined_pca_batch_plot <- pca_plot_pre + pca_plot_aft + plot_layout(ncol = 1)

pdf("合并数据去批次结果/6.去除批次后的PCA图（按批次给色）.pdf", width = pca_plot_width, height = pca_plot_height)
print(combined_pca_batch_plot)
dev.off()

cat("12. 去除批次后的 PCA 图（按批次给色）已绘制完毕并保存...\n")

# 27. 绘制去批次后 PCA 图（按分组给色）
pca_scores$Group <- group_labels

pca_plot_aft_group <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = pca_point_size) +
  stat_ellipse(aes(fill = Group, color = Group), geom = "polygon", level = 0.95, alpha = pca_alpha_value) +
  scale_color_manual(values = group_colors_use) +
  scale_fill_manual(values = group_colors_use) +
  theme_minimal() +
  labs(
    title = "PCA of Batch-Corrected Data by Group",
    x = paste0("Principal Component 1 (", explained_variance_corrected[1], "%)"),
    y = paste0("Principal Component 2 (", explained_variance_corrected[2], "%)"),
    color = "Group",
    fill = "Group"
  ) +
  theme(
    plot.title = pca_theme_plot_title,
    axis.text = pca_theme_axis_text
  )

combined_pca_group_plot <- pca_plot_pre_group + pca_plot_aft_group + plot_layout(ncol = 1)

pdf("合并数据去批次结果/7.去除批次后的PCA图（按分组给色）.pdf", width = pca_plot_width, height = pca_plot_height)
print(combined_pca_group_plot)
dev.off()

cat("13. 去除批次后的 PCA 图（按分组给色）已绘制完毕并保存...\n")

# 28. 保存校正后的数据
write.csv(df_corrected, "合并数据去批次结果/14.去除批次后的数据.csv", row.names = TRUE)
cat("14. 去除批次后的数据已保存，一键式分析结束...\n")
