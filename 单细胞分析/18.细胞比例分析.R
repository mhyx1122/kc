suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(tidyr)
  library(ggalluvial)
  library(reshape2)
})

# 1. 检查全局环境中是否存在 seurat 对象
if (!exists("seurat", envir = .GlobalEnv)) {
  stop("全局环境中没有 seurat 对象，请先加载。")
}

srt <- get("seurat", envir = .GlobalEnv)

if (!"orig.ident" %in% colnames(srt@meta.data)) {
  stop("seurat@meta.data 中没有 orig.ident 列。")
}


# 2. 输出目录模块

# 2.1 输出目录参数
out_dir <- "细胞比例分析"

# 2.2 创建输出目录
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}


# 3. 样本分组设置模块

# 3.1 获取样本名称
sample_names <- unique(trimws(as.character(srt@meta.data$orig.ident)))

# 3.2 设置每个样本对应的分组
# 注意：sample_groups 的顺序必须和 sample_names 一一对应
# 先运行 sample_names 查看样本顺序，再修改 sample_groups
sample_groups <- sample_names

# 示例：
# sample_groups <- c("Normal", "Normal", "Disease", "Disease")

# 3.3 检查分组数量是否匹配
if (length(sample_groups) != length(sample_names)) {
  stop("sample_groups 的数量必须和 sample_names 的数量一致。")
}

if (any(trimws(sample_groups) == "")) {
  stop("sample_groups 中不能有空值。")
}

# 3.4 生成样本到分组的映射
group_mapping <- stats::setNames(
  trimws(as.character(sample_groups)),
  trimws(as.character(sample_names))
)

unique_groups <- unique(sample_groups)


# 4. 公共分析参数模块

# 4.1 差异分析方法
# 可选："t.test"、"wilcox.test"、"anova"
analysis_method <- "wilcox.test"

# 4.2 细胞类型颜色集合
celltype_palette_text <- "#E64B35FF,#4DBBD5FF,#00A087FF,#3C5488FF,#F39B7FFF,#8491B4FF,#91D1C2FF,#7E6148FF,#0073C2FF,#EFC000FF,#868686FF,#CD534CFF,#7AA6DCFF,#003C67FF,#8F7700FF,#3B3B3BFF,#A73030FF,#374E55FF,#DF8F44FF,#00A1D5FF,#B24745FF,#79AF97FF,#6A6599FF,#80796BFF,#BC3C29FF,#0072B5FF,#E18727FF,#20854EFF,#7876B1FF,#6F99ADFF,#FFDC91FF,#EE4C97FF"

celltype_palette <- unlist(strsplit(celltype_palette_text, ","))
celltype_palette <- trimws(celltype_palette)
celltype_palette <- celltype_palette[celltype_palette != ""]

# 4.3 分组颜色集合
group_palette_text <- "#E64B35FF,#4DBBD5FF,#00A087FF,#3C5488FF,#F39B7FFF,#8491B4FF,#91D1C2FF,#7E6148FF"

group_palette <- unlist(strsplit(group_palette_text, ","))
group_palette <- trimws(group_palette)
group_palette <- group_palette[group_palette != ""]


# 5. 为每个细胞映射分组模块

# 5.1 根据 orig.ident 映射 group
cell_group <- unname(group_mapping[trimws(as.character(srt$orig.ident))])

# 5.2 检查是否存在未匹配样本
if (any(is.na(cell_group))) {
  missing_samples <- unique(trimws(as.character(srt$orig.ident))[is.na(cell_group)])
  stop(paste0("以下样本未成功匹配到组别：", paste(missing_samples, collapse = ", ")))
}

# 5.3 检查长度是否一致
if (length(cell_group) != ncol(srt)) {
  stop("分组向量长度与 Seurat 对象细胞数不一致。")
}

# 5.4 写入 Seurat 对象
srt$group <- cell_group


# 6. 计算各样本细胞比例模块

# 6.1 按 Idents(srt) 和 orig.ident 计算每个样本内细胞比例
sample_prop <- prop.table(
  table(Idents(srt), srt$orig.ident),
  margin = 2
)

sample_prop <- as.data.frame(sample_prop)
colnames(sample_prop) <- c("CellName", "Sample", "Proportion")

# 6.2 添加样本所属分组
sample_prop$Group <- unname(group_mapping[trimws(as.character(sample_prop$Sample))])

# 6.3 保存各样本细胞比例表
name_table <- "各样本的细胞比例"

write.csv(
  sample_prop,
  file = file.path(out_dir, paste0(name_table, ".csv")),
  row.names = FALSE
)


# 7. 计算按分组汇总的比例模块

# 7.1 构建分组和细胞类型数据框
df_group <- data.frame(
  CellType = Idents(srt),
  Group = srt$group
)

# 7.2 按 Group 归一化：每个组别内各细胞类型比例和为 1
df_group_count <- df_group %>%
  dplyr::group_by(Group, CellType) %>%
  dplyr::summarise(Freq = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(Group) %>%
  dplyr::mutate(Proportion = Freq / sum(Freq)) %>%
  dplyr::ungroup()

# 7.3 按 CellType 归一化：每个细胞类型在不同组别中的比例和为 1
df_group_count1 <- df_group %>%
  dplyr::group_by(CellType, Group) %>%
  dplyr::summarise(Freq = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(CellType) %>%
  dplyr::mutate(Proportion = Freq / sum(Freq)) %>%
  dplyr::ungroup()

# 7.4 保存分组汇总比例表
write.csv(
  df_group_count,
  file = file.path(out_dir, "分组汇总细胞比例.csv"),
  row.names = FALSE
)

write.csv(
  df_group_count1,
  file = file.path(out_dir, "分组汇总细胞比例2.csv"),
  row.names = FALSE
)


# 8. 准备流向图数据模块

# 8.1 宽表
alluvial_wide <- df_group_count %>%
  dplyr::select(Group, CellType, Proportion) %>%
  tidyr::pivot_wider(
    names_from = CellType,
    values_from = Proportion,
    values_fill = 0
  ) %>%
  as.data.frame()

# 8.2 长表
alluvial_long <- reshape2::melt(
  alluvial_wide,
  id.vars = "Group"
)

colnames(alluvial_long) <- c("Group", "CellType", "value")


# 9. 单样本细胞比例堆叠图模块

# 9.1 图1参数
plot1_border_linewidth <- 0.5

plot1_width <- 12
plot1_height <- 8
name_plot1 <- "1.单样本的堆叠图"

# 9.2 准备颜色
cell_names_plot1 <- unique(as.character(sample_prop$CellName))
cell_pal_plot1 <- rep(celltype_palette, length.out = length(cell_names_plot1))
names(cell_pal_plot1) <- cell_names_plot1

# 9.3 生成图1
p_sample_stack <- ggplot(sample_prop) +
  geom_bar(
    aes(x = Sample, y = Proportion, fill = CellName),
    stat = "identity",
    width = 0.7,
    linewidth = 0.5,
    colour = "#222222"
  ) +
  theme_classic() +
  scale_fill_manual(values = cell_pal_plot1) +
  labs(x = "Sample", y = "Ratio") +
  coord_flip() +
  theme(
    panel.border = element_rect(
      fill = NA,
      color = "black",
      linewidth = plot1_border_linewidth,
      linetype = "solid"
    )
  )

# 9.4 保存图1
ggsave(
  filename = file.path(out_dir, paste0(name_plot1, ".pdf")),
  plot = p_sample_stack,
  width = plot1_width,
  height = plot1_height,
  device = "pdf"
)


# 10. 分组堆叠图模块

# 10.1 图2参数
plot2_width <- 12
plot2_height <- 8
name_plot2 <- "分组堆叠图"

# 10.2 准备颜色
group_names_plot2 <- unique(as.character(df_group_count1$Group))
group_pal_plot2 <- rep(group_palette, length.out = length(group_names_plot2))
names(group_pal_plot2) <- group_names_plot2

# 10.3 生成图2
p_group_stack <- ggplot(
  df_group_count1,
  aes(x = CellType, y = Proportion, fill = Group)
) +
  geom_bar(stat = "identity", width = 0.7, color = "black") +
  scale_fill_manual(values = group_pal_plot2) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 10.4 保存图2
ggsave(
  filename = file.path(out_dir, paste0(name_plot2, ".pdf")),
  plot = p_group_stack,
  width = plot2_width,
  height = plot2_height,
  device = "pdf"
)


# 11. 分组差异箱线图模块

# 11.1 图3参数
axis_text_size <- 9
legend_text_size <- 9
signif_label_y <- 0.6

plot3_width <- 12
plot3_height <- 8
name_plot3 <- "2.分组间的差异箱线图"

# 11.2 准备颜色
group_names_plot3 <- unique(as.character(sample_prop$Group))
group_pal_plot3 <- rep(group_palette, length.out = length(group_names_plot3))
names(group_pal_plot3) <- group_names_plot3

# 11.3 生成图3
p_box <- ggplot(
  sample_prop,
  aes(x = CellName, y = Proportion, fill = Group)
) +
  geom_boxplot(
    width = 0.3,
    position = position_dodge(0.5),
    outlier.colour = NA
  ) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_fill_manual(values = group_pal_plot3) +
  stat_compare_means(
    aes(group = Group),
    method = analysis_method,
    label = "p.signif",
    hide.ns = TRUE,
    label.y = signif_label_y
  ) +
  theme(
    axis.text.x = element_text(size = axis_text_size, colour = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = axis_text_size, angle = 0),
    axis.title.x = element_text(size = 5),
    axis.title.y = element_text(size = axis_text_size),
    legend.text = element_text(size = legend_text_size)
  ) +
  xlab("")

# 11.4 保存图3
ggsave(
  filename = file.path(out_dir, paste0(name_plot3, ".pdf")),
  plot = p_box,
  width = plot3_width,
  height = plot3_height,
  device = "pdf"
)


# 12. 细胞在不同组别中比例的流向图模块

# 12.1 图4参数
highlight_celltype <- unique(as.character(alluvial_long$CellType))[1]

curve_type <- "linear"
# 可选：
# curve_type <- "cubic"
# curve_type <- "sigmoid"

alluvial_border_color <- "darkgreen"
alluvial_border_size <- 0.5
show_axis_lines <- FALSE
show_axis_ticks <- TRUE

plot4_width <- 12
plot4_height <- 8
name_plot4 <- "3.细胞在不同组别中比例的流向图"

# 12.2 准备颜色
cell_names_plot4 <- unique(as.character(alluvial_long$CellType))
cell_pal_plot4 <- rep(celltype_palette, length.out = length(cell_names_plot4))
names(cell_pal_plot4) <- cell_names_plot4

# 12.3 检查高亮细胞类型
if (is.null(highlight_celltype) || highlight_celltype == "") {
  highlight_celltype <- cell_names_plot4[1]
}

# 12.4 生成图4
p_alluvial <- ggplot(
  alluvial_long,
  aes(
    x = Group,
    y = value,
    stratum = CellType,
    alluvium = CellType
  )
) +
  ggalluvial::geom_alluvium(
    aes(fill = CellType, alpha = CellType == highlight_celltype),
    width = 0.3,
    curve_type = curve_type,
    color = alluvial_border_color,
    linewidth = alluvial_border_size
  ) +
  ggalluvial::geom_stratum(
    aes(fill = CellType),
    width = 0.3
  ) +
  scale_fill_manual(values = cell_pal_plot4) +
  scale_alpha_manual(
    values = c("TRUE" = 0.6, "FALSE" = 0),
    guide = "none"
  ) +
  labs(
    x = "Group",
    y = "Percentage",
    fill = "Cell Type",
    title = "Cell Type Proportion Across Groups"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5),
    axis.line = if (show_axis_lines) {
      element_line(color = "black", linewidth = 0.5)
    } else {
      element_blank()
    },
    axis.ticks = if (show_axis_ticks) {
      element_line(color = "black", linewidth = 0.5)
    } else {
      element_blank()
    }
  )

# 12.5 保存图4
ggsave(
  filename = file.path(out_dir, paste0(name_plot4, ".pdf")),
  plot = p_alluvial,
  width = plot4_width,
  height = plot4_height,
  device = "pdf"
)


# 13. 参数记录模块

# 13.1 参数文件保存参数
name_params <- "cell_proportion_parameters"

# 13.2 生成参数记录文本
param_text <- paste0(
  "本次分析参数总结：\n",
  "- 样本分组：\n",
  paste(names(group_mapping), " -> ", unname(group_mapping), collapse = "\n"), "\n",
  "- 差异分析方法：", analysis_method, "\n",
  "- 细胞类型颜色集合：", celltype_palette_text, "\n",
  "- 分组颜色集合：", group_palette_text, "\n",
  "- 图1宽：", plot1_width, " 英寸\n",
  "- 图1高：", plot1_height, " 英寸\n",
  "- 图1边框线宽：", plot1_border_linewidth, "\n",
  "- 图2宽：", plot2_width, " 英寸\n",
  "- 图2高：", plot2_height, " 英寸\n",
  "- 图3宽：", plot3_width, " 英寸\n",
  "- 图3高：", plot3_height, " 英寸\n",
  "- 图3轴标签字体大小：", axis_text_size, "\n",
  "- 图3图例字体大小：", legend_text_size, "\n",
  "- 图3显著性标注位置：", signif_label_y, "\n",
  "- 图4高亮细胞类型：", highlight_celltype, "\n",
  "- 图4曲线类型：", curve_type, "\n",
  "- 图4流向带边框颜色：", alluvial_border_color, "\n",
  "- 图4流向带边框粗细：", alluvial_border_size, "\n",
  "- 图4显示轴线：", show_axis_lines, "\n",
  "- 图4显示刻度线：", show_axis_ticks, "\n",
  "- 样本数：", length(sample_names), "\n",
  "- 分组数：", length(unique_groups), "\n",
  "- 输出目录：", out_dir, "\n"
)

# 13.3 保存参数记录
writeLines(
  param_text,
  con = file.path(out_dir, paste0(name_params, ".txt"))
)


# 14. 写回 seurat 对象

seurat <- srt