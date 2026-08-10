library(ggplot2)
library(cols4all)

# 1. 设置输入参数
input_file <- "your_data.csv"
color_palette <- "br_bg"
output_dir <- "定制图片"
plot_width <- 8
plot_height <- 6

# 2. 创建输出文件夹
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# 3. 读取CSV文件
dt <- read.csv(input_file,header = TRUE)

# 4. 计算 -log10pvalue
dt$"-log10pvalue" <- ifelse(
  dt$change == "up",
  -log10(dt$pvalue),
  -(-log10(dt$pvalue))
)

# 5. 设置 Description 的顺序
level_up <- dt$Description[dt$change == "up"]
level_down <- dt$Description[dt$change == "down"]
level <- c(level_up,level_down)
dt$Description <- factor(dt$Description,levels = rev(level))

# 6. 绘制棒棒糖图
p <- ggplot(dt,aes(x = `-log10pvalue`,y = Description)) +
  geom_col(
    aes(fill = `-log10pvalue`),
    width = 0.1
  ) +
  geom_point(
    aes(size = Count,color = `-log10pvalue`)
  ) +
  scale_size_continuous(
    range = c(2,7)
  ) +
  scale_color_continuous_c4a_div(
    color_palette,
    mid = 0,
    reverse = TRUE,
    labels = function(x) abs(x)
  ) +
  scale_fill_continuous_c4a_div(
    color_palette,
    mid = 0,
    reverse = TRUE,
    labels = function(x) abs(x)
  ) +
  scale_x_continuous(
    breaks = seq(0,100,by = 5),
    labels = function(x) abs(x)
  ) +
  ylab("") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 12),
    axis.title.x = element_text(size = 13),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12)
  )

# 7. 保存图形为PDF文件
output_file <- file.path(output_dir,"lollipop_plot.pdf")
ggsave(output_file,plot = p,width = plot_width,height = plot_height)

# 8. 输出完成提示
message("图像已保存到文件夹：",output_dir)