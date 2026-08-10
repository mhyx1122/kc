library(ggplot2)
library(dplyr)
library(tibble)
library(rlang)
library(gground)
library(ggprism)
library(grid)

# 1. 设置输入参数
input_file <- "your_pathway_data.csv"
p_value_col <- "p.adjust"
category_order <- c("BP", "CC", "MF", "KEGG")
color_palette1122 <- c("#F06292", "#81C784", "#FF8A65", "#42A5F5")
plot_width <- 10
plot_height <- 8
output_dir <- "定制图片"
color_palette <- color_palette1122

# 2. 读取数据
pathway_data <- read.csv(input_file,check.names = FALSE)

# 3. 创建输出文件夹
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# 4. 数据处理
pathway_data <- pathway_data %>%
  ungroup() %>%
  mutate(ONTOLOGY = factor(ONTOLOGY,levels = rev(category_order))) %>%
  group_by(ONTOLOGY) %>%
  dplyr::arrange(ONTOLOGY,-log10(!!sym(p_value_col))) %>%
  mutate(Description = factor(Description,levels = Description)) %>%
  tibble::rowid_to_column("index") %>%
  ungroup()

# 5. 设置分类标签和基因数量点图参数
width_value <- 0.5
xaxis_limit <- max(-log10(pathway_data[[p_value_col]])) + 1

# 6. 创建分类标签数据框
rect_data <- group_by(pathway_data,ONTOLOGY) %>%
  reframe(n = n()) %>%
  ungroup() %>%
  mutate(
    xmin = -3 * width_value,
    xmax = -2 * width_value,
    ymax = cumsum(n),
    ymin = lag(ymax,default = 0) + 0.6,
    ymax = ymax + 0.4
  )

# 7. 创建PDF文件路径
pdf_file_path <- file.path(output_dir,"KEGG和GO进阶版.pdf")

# 8. 输出PDF图形
pdf(pdf_file_path,width = plot_width,height = plot_height)

# 9. 作图
p <- pathway_data %>%
  ggplot(aes(-log10(!!sym(p_value_col)),y = index,fill = ONTOLOGY)) +
  geom_round_col(
    aes(y = Description),
    width = 0.6,
    alpha = 0.8
  ) +
  geom_text(
    aes(x = 0.05,label = Description),
    hjust = 0,
    size = 5
  ) +
  geom_text(
    aes(x = 0.1,label = geneID,colour = ONTOLOGY),
    hjust = 0,
    vjust = 2.6,
    size = 3.5,
    fontface = "italic",
    show.legend = FALSE
  ) +
  geom_point(
    aes(x = -width_value,size = Count),
    shape = 21
  ) +
  geom_text(
    aes(x = -width_value,label = Count)
  ) +
  scale_size_continuous(
    name = "Count",
    range = c(5,16)
  ) +
  geom_round_rect(
    aes(xmin = xmin,xmax = xmax,ymin = ymin,ymax = ymax,fill = ONTOLOGY),
    data = rect_data,
    radius = unit(2,"mm"),
    inherit.aes = FALSE
  ) +
  geom_text(
    aes(x = (xmin + xmax) / 2,y = (ymin + ymax) / 2,label = ONTOLOGY),
    data = rect_data,
    inherit.aes = FALSE,
    angle = 90
  ) +
  annotate(
    "segment",
    x = 0,
    xend = xaxis_limit,
    y = 0,
    yend = 0,
    linewidth = 1.5,
    color = "black"
  ) +
  labs(y = NULL) +
  scale_fill_manual(
    name = "Category",
    values = color_palette,
    breaks = category_order
  ) +
  scale_colour_manual(
    values = color_palette,
    breaks = category_order
  ) +
  scale_x_continuous(
    breaks = seq(0,xaxis_limit,2),
    expand = expansion(c(0,0))
  ) +
  theme_prism() +
  theme(
    axis.text.y = element_blank(),
    axis.line = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text()
  )

# 10. 打印并保存图形
print(p)
dev.off()

# 11. 输出完成提示
message("图像已生成并保存在文件夹: ", output_dir)