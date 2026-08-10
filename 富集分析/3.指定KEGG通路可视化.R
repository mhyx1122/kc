library(ggplot2)
library(dplyr)

# 1. 设置输入文件路径
file_path <- "KEGG_result.csv"

# 2. 设置图像大小
Figwidth <- 7
Figheight <- 5

# 3. 设置 low 和 high 颜色
low_color <- "#B2182B"   # 小P值，显著，深红
high_color <- "#2166AC"  # 大P值，不显著，蓝色
# 4. 选择颜色映射依据
# 可选："pvalue" 或 "p.adjust"
color_choice <- "pvalue"

# 5. 设置保存文件夹
save_folder <- "指定通路可视化"

# 6. 读取数据
rt <- read.csv(file_path, header = TRUE, check.names = FALSE)

# 7. 计算 GeneRatio
rt$GeneRatio <- sapply(strsplit(rt$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))

# 8. 按照 GeneRatio 对 Description 排序
labels <- rt[order(rt$GeneRatio), "Description"]
rt$Description <- factor(rt$Description, levels = labels)

# 9. 设置颜色映射变量和颜色
color_var <- color_choice
low_col <- low_color
high_col <- high_color

# 10. 创建 ggplot 图形
p <- ggplot(rt, aes(x = GeneRatio, y = Description)) +
  geom_point(aes(size = Count, color = !!sym(color_var))) +
  scale_colour_gradient(low = high_col, high = low_col) +
  labs(color = color_var, size = "Count", x = "Gene ratio", y = "") +
  theme_bw() +
  theme(
    axis.text.x = element_text(color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# 11. 检查并创建保存文件夹
if (!dir.exists(save_folder)) {
  dir.create(save_folder, recursive = TRUE)
}

# 12. 设置输出文件路径
output_file <- file.path(save_folder, "KEGG指定图.pdf")

# 13. 保存图像为 PDF 文件
ggsave(output_file, plot = p, width = Figwidth, height = Figheight)

# 14. 输出完成提示
cat("数据处理完成。图像已保存为 PDF 文件：", output_file, "\n")
