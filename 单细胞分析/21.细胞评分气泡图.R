# 1. 加载必要 R 包

library(ggplot2)
library(dplyr)
library(tidyr)


# 2. 读取输入数据

# 2.1 输入文件参数
input_file <- "your_input_file.csv"

# 2.2 读取 CSV 文件
data <- read.csv(
  input_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

# 2.3 检查是否存在 cellType 列
if (!"cellType" %in% colnames(data)) {
  stop("输入数据中未找到 cellType 列，请检查 CSV 文件格式。")
}


# 3. 转换数据格式并计算每种细胞类型的平均评分

# 3.1 将宽格式数据转换成长格式
data_long <- data %>%
  gather(
    key = "score_type",
    value = "score",
    -cellType
  )

# 3.2 确保 score 为数值型
data_long$score <- as.numeric(data_long$score)

# 3.3 按 cellType 和 score_type 计算平均评分
data_long <- data_long %>%
  group_by(cellType, score_type) %>%
  summarise(
    mean_score = mean(score, na.rm = TRUE),
    .groups = "drop"
  )


# 4. 对不同评分方法进行量纲统一

# 4.1 定义 Min-Max 标准化函数
# 说明：
# 每一种 score_type 内部单独标准化
# 标准化后，最低值为 0，最高值为 1
# 如果某一种 score_type 内所有 cellType 的均值完全一样，则统一设为 0.5
min_max_scale <- function(x) {
  x_min <- min(x, na.rm = TRUE)
  x_max <- max(x, na.rm = TRUE)
  
  if (is.na(x_min) || is.na(x_max)) {
    return(rep(NA, length(x)))
  }
  
  if (x_max == x_min) {
    return(rep(0.5, length(x)))
  }
  
  return((x - x_min) / (x_max - x_min))
}

# 4.2 在每一种评分方法内部进行 Min-Max 标准化
data_long <- data_long %>%
  group_by(score_type) %>%
  mutate(
    mean_score_scaled = min_max_scale(mean_score)
  ) %>%
  ungroup()


# 5. 基于标准化后的评分计算气泡大小比例

# 5.1 每一种评分方法内部，计算标准化评分总和
data_long <- data_long %>%
  group_by(score_type) %>%
  mutate(
    total_scaled_score = sum(mean_score_scaled, na.rm = TRUE)
  ) %>%
  ungroup()

# 5.2 计算每个 cellType 在对应评分方法中的相对比例
data_long <- data_long %>%
  mutate(
    proportion = mean_score_scaled / total_scaled_score * 100
  )

# 5.3 如果某些异常情况导致 proportion 为 NA，则设为 0
data_long$proportion[is.na(data_long$proportion)] <- 0


# 6. 绘制气泡图

# 6.1 气泡图参数
point_color <- "#6F99ADFF"
border_color <- "#E64B35FF"

# 6.2 绘图
p <- ggplot(
  data_long,
  aes(
    x = score_type,
    y = cellType,
    size = proportion,
    color = mean_score_scaled
  )
) +
  geom_point(
    alpha = 0.6,
    shape = 16,
    stroke = 1
  ) +
  scale_size_continuous(
    name = "Relative proportion (%)",
    range = c(3, 12)
  ) +
  scale_color_gradient(
    low = point_color,
    high = border_color,
    name = "Scaled mean score"
  ) +
  labs(
    title = "Bubble Plot of Scaled Scores by Cell Type",
    x = "",
    y = ""
  ) +
  theme_minimal(base_size = 15) +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    plot.title = element_text(
      hjust = 0.5,
      size = 18,
      face = "bold"
    ),
    axis.title = element_text(size = 14),
    axis.text.x = element_text(
      size = 12,
      angle = 45,
      hjust = 1
    ),
    axis.text.y = element_text(size = 12)
  )

# 6.3 显示图片
print(p)


# 7. 保存气泡图

# 7.1 输出参数
output_dir <- "气泡图输出"
plot_width <- 8
plot_height <- 6
pdf_name <- "bubble_plot_scaled_output.pdf"

# 7.2 创建输出文件夹
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 7.3 设置 PDF 保存路径
pdf_path <- file.path(output_dir, pdf_name)

# 7.4 保存 PDF
pdf(
  file = pdf_path,
  width = plot_width,
  height = plot_height
)

print(p)

dev.off()

cat("气泡图已保存至：", pdf_path, "\n")


# 8. 保存整理后的绘图数据

# 8.1 输出标准化后的数据表
output_table <- file.path(output_dir, "bubble_plot_scaled_data.csv")

write.csv(
  data_long,
  file = output_table,
  row.names = FALSE
)

cat("绘图数据已保存至：", output_table, "\n")