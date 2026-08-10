library(ggplot2)
library(ggh4x)
library(stringr)
library(grid)

# 1. 设置输入文件路径
file_path <- "GO_result2.csv"

# 2. 设置输出文件夹
save_folder <- "指定通路可视化"

# 3. 设置图像大小
Figwidth <- 8
Figheight <- 6

# 4. 设置 GO 三大类别颜色
color_Biological <- "#B02DE2"
color_Molecular <- "#F07A48"
color_Cellular <- "#F0456A"

# 5. 读取数据
data <- read.csv(file_path, header = TRUE, check.names = FALSE)

# 6. 替换 ONTOLOGY 列中的简写
data$ONTOLOGY <- gsub("BP", "Biological Process", data$ONTOLOGY)
data$ONTOLOGY <- gsub("MF", "Molecular Function", data$ONTOLOGY)
data$ONTOLOGY <- gsub("CC", "Cellular Component", data$ONTOLOGY)

# 7. 设置每个分面的背景颜色
facet_background_colors <- c(
  "Biological Process" = color_Biological,
  "Molecular Function" = color_Molecular,
  "Cellular Component" = color_Cellular
)

# 8. 检查并创建保存文件夹
if (!dir.exists(save_folder)) {
  dir.create(save_folder, recursive = TRUE)
}

# 9. 创建 ggplot 图形
p <- ggplot(data, aes(x = Count, y = reorder(Description, Count), fill = ONTOLOGY)) +
  geom_bar(stat = "identity", color = "black", width = 0.8) +
  scale_fill_manual(values = facet_background_colors) +
  facet_grid2(
    rows = vars(ONTOLOGY),
    scales = "free_y",
    space = "free_y",
    strip = strip_themed(
      background_x = lapply(
        facet_background_colors,
        function(col) element_rect(fill = col, color = "black")
      )
    )
  ) +
  theme_classic() +
  labs(x = "Number of Genes", y = "", title = "") +
  theme(
    strip.text.y.left = element_text(angle = 90, hjust = 0, face = "bold", lineheight = 1.2),
    strip.placement = "outside",
    panel.spacing.y = unit(0.0, "lines"),
    legend.position = "none",
    axis.text.x = element_text(angle = 0, hjust = 1, colour = "black", size = 10),
    axis.text.y = element_text(colour = "black", size = 10)
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  geom_text(aes(label = Count), hjust = 1.1, size = 3.5) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 30))

# 10. 设置输出路径
output_file <- file.path(save_folder, "指定的GO图.pdf")

# 11. 保存 PDF
pdf(file = output_file, width = Figwidth, height = Figheight)
par(mar = c(6, 10, 3, 3))
print(p)
dev.off()

# 12. 输出完成提示
cat("数据处理完成。图像已保存为 PDF 文件：", output_file, "\n")