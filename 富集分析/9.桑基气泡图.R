library(ggplot2)
library(dplyr)
library(tidyr)
library(ggalluvial)
library(cowplot)
library(grid)

# 1. 设置输入参数
input_file <- "your_data.csv"
working_dir <- "定制图片"
plot_width <- 9
plot_height <- 8
colorlow <- "#e06738"
colorhigh <- "#3fc039"

# 2. 读取数据
data2 <- read.csv(input_file,header = TRUE)

# 3. 创建输出文件夹
if (!dir.exists(working_dir)) {
  dir.create(working_dir)
}

# 4. 准备原始数据
data <- data2

# 5. 转换数据为适合桑基图的格式
data_long <- data %>% tidyr::separate_rows(geneID,sep = "/")
data_long$geneID <- factor(data_long$geneID,levels = unique(data_long$geneID))
data_long$Description <- factor(data_long$Description,levels = unique(data$Description))

# 6. 计算Count列并按降序排序Description
data_long <- data_long %>%
  group_by(Description) %>%
  mutate(Count = n()) %>%
  ungroup()

# 7. 按Count值降序排列Description
data_long <- data_long %>% arrange(desc(Count))
data_long$Description <- factor(data_long$Description,levels = unique(data_long$Description))

# 8. 转换成长格式
data_l <- data_long[,c("geneID","Description")]
data_lodes <- to_lodes_form(data_l,key = "x",value = "stratum",id = "alluvium",axes = 1:2)

# 9. 绘制桑基图
sankey_plot <- ggplot(data_lodes,aes(x = x,stratum = stratum,alluvium = alluvium,fill = stratum,label = stratum)) +
  scale_x_discrete(expand = c(0,0)) +
  geom_flow(alpha = 0.3,width = 0,knot.pos = 0.1) +
  geom_stratum(width = 0.05,color = "white") +
  geom_text(stat = "stratum",aes(label = after_stat(stratum)),size = 3,hjust = 1,nudge_x = -0.03) +
  guides(fill = "none",color = "none") +
  expand_limits(x = c(0.79,1)) +
  theme_minimal() +
  labs(title = "",x = "",y = "") +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = unit(c(0,9,0,0),units = "cm")
  )

# 10. 绘制气泡图
data_long$Description <- factor(data_long$Description,levels = rev(unique(data_long$Description)))
data_long <- data_long %>% arrange(Description) %>% mutate(Description_num = cumsum(Count))
dot_plot <- ggplot(data_long,aes(x = Count,y = Description)) +
  geom_point(aes(size = Count,color = pvalue)) +
  scale_color_gradient(low = colorlow,high = colorhigh) +
  theme_classic() +
  labs(size = "Count",color = "Pvalue",y = "",x = "Count") +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

# 11. 拼图
final_plot <- ggdraw() +
  draw_plot(sankey_plot,0,0,1,1) +
  draw_plot(dot_plot,0.6,0.05,0.4,0.88)

# 12. 设置输出文件路径并保存PDF
pdf_path <- file.path(working_dir,"桑基图联合气泡图.pdf")
ggsave(pdf_path,plot = final_plot,device = "pdf",width = plot_width,height = plot_height)

# 13. 输出完成提示
message("数据处理完成。图像已保存为 PDF 文件：",pdf_path)