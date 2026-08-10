library(ggplot2)
library(dplyr)
library(Seurat)
library(tidydr)
library(grid)

# 1. 检查全局环境中是否存在 seurat 对象
if (!exists("seurat", envir = .GlobalEnv)) {
  stop("全局环境中没有 seurat 对象")
}

srt <- get("seurat", envir = .GlobalEnv)


# 2. 输出目录模块

# 2.1 输出目录参数
out_dir_default <- "4.2带置信区间的UMAP和TSNE"
save_dir <- out_dir_default

# 2.2 创建输出目录
if (!dir.exists(save_dir)) {
  dir.create(save_dir, recursive = TRUE)
}


# 3. 绘图参数模块

# 3.1 基础绘图参数
line_size <- 0.01
bxalpha <- 0.4
point_size <- 0.5
line_type_1 <- 1
label_size <- 3

# 3.2 是否显示置信椭圆
# TRUE：显示95%置信椭圆
# FALSE：不显示95%置信椭圆
show_ellipse <- TRUE

# 3.3 分面参数
# 空字符串表示不分面
# "group" 表示按 seurat@meta.data$group 分面
# "orig.ident" 表示按 seurat@meta.data$orig.ident 分面
facet_column <- ""

# 3.4 分面排版参数
facet_ncol <- 2
facet_panel_spacing_x <- 1.5
facet_panel_spacing_y <- 1.5

# 3.5 tSNE 保存参数
w_tsne <- 12
h_tsne <- 6
name_tsne <- "tSNE_plot"

# 3.6 UMAP 保存参数
w_umap <- 12
h_umap <- 6
name_umap <- "UMAP_plot"

# 3.7 tSNE 颜色参数
palette_tsne_text <- "#E64B35FF,#4DBBD5FF,#00A087FF,#3C5488FF,#F39B7FFF,#8491B4FF,#91D1C2FF,#7E6148FF,#0073C2FF,#EFC000FF,#868686FF,#CD534CFF,#7AA6DCFF,#003C67FF,#8F7700FF,#3B3B3BFF,#A73030FF,#374E55FF,#DF8F44FF,#00A1D5FF,#B24745FF,#79AF97FF,#6A6599FF,#80796BFF,#BC3C29FF,#0072B5FF,#E18727FF,#20854EFF,#7876B1FF,#6F99ADFF,#FFDC91FF,#EE4C97FF"

palette_tsne <- unlist(strsplit(palette_tsne_text, ","))
palette_tsne <- trimws(palette_tsne)
palette_tsne <- palette_tsne[palette_tsne != ""]

# 3.8 UMAP 颜色参数
palette_umap_text <- "#E64B35FF,#4DBBD5FF,#00A087FF,#3C5488FF,#F39B7FFF,#8491B4FF,#91D1C2FF,#7E6148FF,#0073C2FF,#EFC000FF,#868686FF,#CD534CFF,#7AA6DCFF,#003C67FF,#8F7700FF,#3B3B3BFF,#A73030FF,#374E55FF,#DF8F44FF,#00A1D5FF,#B24745FF,#79AF97FF,#6A6599FF,#80796BFF,#BC3C29FF,#0072B5FF,#E18727FF,#20854EFF,#7876B1FF,#6F99ADFF,#FFDC91FF,#EE4C97FF"

palette_umap <- unlist(strsplit(palette_umap_text, ","))
palette_umap <- trimws(palette_umap)
palette_umap <- palette_umap[palette_umap != ""]


# 4. 检查 tSNE、UMAP 和 meta.data 列

# 4.1 检查 tSNE 是否存在
if (!("tsne" %in% names(srt@reductions))) {
  stop("seurat@reductions 中不存在：tsne（请先 RunTSNE）")
}

# 4.2 检查 UMAP 是否存在
if (!("umap" %in% names(srt@reductions))) {
  stop("seurat@reductions 中不存在：umap（请先 RunUMAP）")
}

# 4.3 检查 cellType 列是否存在
if (!("cellType" %in% colnames(srt@meta.data))) {
  stop("seurat@meta.data 中不存在列：cellType（请确保用于分组的列名为 cellType）")
}

# 4.4 检查分面列是否存在
facet_column <- trimws(facet_column)

if (facet_column != "" && !(facet_column %in% colnames(srt@meta.data))) {
  stop(paste0("seurat@meta.data 中不存在分面列：", facet_column))
}


# 5. 提取 tSNE、UMAP 坐标和 meta.data 模块

# 5.1 提取 tSNE 坐标
tSNE <- as.data.frame(srt@reductions$tsne@cell.embeddings)

# 5.2 提取 UMAP 坐标
UMAP <- as.data.frame(srt@reductions$umap@cell.embeddings)

# 5.3 提取 meta.data
meta <- srt@meta.data

# 5.4 统一 tSNE 和 UMAP 的前两列列名
if (ncol(tSNE) >= 2) {
  colnames(tSNE)[1:2] <- c("tSNE_1", "tSNE_2")
}

if (ncol(UMAP) >= 2) {
  colnames(UMAP)[1:2] <- c("umap_1", "umap_2")
}


# 6. 合并坐标和 meta.data 模块

# 6.1 合并 tSNE 坐标和 meta.data
tSNE_data <- merge(
  tSNE,
  meta,
  by = "row.names",
  drop = TRUE
)

rownames(tSNE_data) <- tSNE_data$Row.names
tSNE_data <- tSNE_data[, setdiff(colnames(tSNE_data), "Row.names"), drop = FALSE]

# 6.2 合并 UMAP 坐标和 meta.data
UMAP_data <- merge(
  UMAP,
  meta,
  by = "row.names",
  drop = TRUE
)

rownames(UMAP_data) <- UMAP_data$Row.names
UMAP_data <- UMAP_data[, setdiff(colnames(UMAP_data), "Row.names"), drop = FALSE]


# 7. tSNE 图模块

# 7.1 设置 tSNE 分组颜色
cell_levels_tsne <- sort(unique(as.character(tSNE_data$cellType)))
cols_tsne <- rep(palette_tsne, length.out = length(cell_levels_tsne))
names(cols_tsne) <- cell_levels_tsne

# 7.2 计算 tSNE 每个 cellType 的标签中心点
if (facet_column == "") {
  
  centers_tsne <- tSNE_data %>%
    group_by(cellType) %>%
    summarise(
      tSNE_1 = median(tSNE_1),
      tSNE_2 = median(tSNE_2),
      .groups = "drop"
    )
  
} else {
  
  centers_tsne <- tSNE_data %>%
    group_by(.data[[facet_column]], cellType) %>%
    summarise(
      tSNE_1 = median(tSNE_1),
      tSNE_2 = median(tSNE_2),
      .groups = "drop"
    )
}

# 7.3 生成 tSNE 基础图
p_tsne <- ggplot(tSNE_data, aes(x = tSNE_1, y = tSNE_2)) +
  geom_point(aes(color = cellType), size = point_size, alpha = 0.8) +
  theme_void() +
  theme_dr(
    xlength = 0.2,
    ylength = 0.2,
    arrow = grid::arrow(length = unit(0.1, "inches"), ends = "last", type = "closed")
  ) +
  theme(
    panel.grid = element_blank(),
    plot.margin = margin(10, 20, 10, 10),
    panel.spacing.x = unit(facet_panel_spacing_x, "cm"),
    panel.spacing.y = unit(facet_panel_spacing_y, "cm"),
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_blank()
  ) +
  geom_text(
    data = centers_tsne,
    aes(x = tSNE_1, y = tSNE_2, label = cellType),
    fontface = "bold",
    color = "black",
    size = label_size
  ) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_color_manual(values = cols_tsne) +
  scale_fill_manual(values = cols_tsne)

# 7.4 根据参数决定是否显示 tSNE 置信椭圆
if (show_ellipse) {
  p_tsne <- p_tsne +
    stat_ellipse(
      aes(color = cellType, fill = cellType),
      level = 0.95,
      linetype = as.numeric(line_type_1),
      linewidth = line_size,
      show.legend = FALSE,
      geom = "polygon",
      alpha = bxalpha
    )
}

# 7.5 根据参数决定是否分面
if (facet_column != "") {
  p_tsne <- p_tsne +
    facet_wrap(
      as.formula(paste("~", facet_column)),
      ncol = facet_ncol,
      scales = "free"
    )
}

# 7.6 保存 tSNE 图
ggsave(
  filename = file.path(save_dir, paste0(name_tsne, ".pdf")),
  plot = p_tsne,
  width = w_tsne,
  height = h_tsne,
  device = "pdf"
)


# 8. UMAP 图模块

# 8.1 设置 UMAP 分组颜色
cell_levels_umap <- sort(unique(as.character(UMAP_data$cellType)))
cols_umap <- rep(palette_umap, length.out = length(cell_levels_umap))
names(cols_umap) <- cell_levels_umap

# 8.2 计算 UMAP 每个 cellType 的标签中心点
if (facet_column == "") {
  
  centers_umap <- UMAP_data %>%
    group_by(cellType) %>%
    summarise(
      umap_1 = median(umap_1),
      umap_2 = median(umap_2),
      .groups = "drop"
    )
  
} else {
  
  centers_umap <- UMAP_data %>%
    group_by(.data[[facet_column]], cellType) %>%
    summarise(
      umap_1 = median(umap_1),
      umap_2 = median(umap_2),
      .groups = "drop"
    )
}

# 8.3 生成 UMAP 基础图
p_umap <- ggplot(UMAP_data, aes(x = umap_1, y = umap_2)) +
  geom_point(aes(color = cellType), size = point_size, alpha = 0.8) +
  theme_void() +
  theme_dr(
    xlength = 0.2,
    ylength = 0.2,
    arrow = grid::arrow(length = unit(0.1, "inches"), ends = "last", type = "closed")
  ) +
  theme(
    panel.grid = element_blank(),
    plot.margin = margin(10, 20, 10, 10),
    panel.spacing.x = unit(facet_panel_spacing_x, "cm"),
    panel.spacing.y = unit(facet_panel_spacing_y, "cm"),
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_blank()
  ) +
  geom_text(
    data = centers_umap,
    aes(x = umap_1, y = umap_2, label = cellType),
    fontface = "bold",
    color = "black",
    size = label_size
  ) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_color_manual(values = cols_umap) +
  scale_fill_manual(values = cols_umap)

# 8.4 根据参数决定是否显示 UMAP 置信椭圆
if (show_ellipse) {
  p_umap <- p_umap +
    stat_ellipse(
      aes(color = cellType, fill = cellType),
      level = 0.95,
      linetype = as.numeric(line_type_1),
      linewidth = line_size,
      show.legend = FALSE,
      geom = "polygon",
      alpha = bxalpha
    )
}

# 8.5 根据参数决定是否分面
if (facet_column != "") {
  p_umap <- p_umap +
    facet_wrap(
      as.formula(paste("~", facet_column)),
      ncol = facet_ncol,
      scales = "free"
    )
}

# 8.6 保存 UMAP 图
ggsave(
  filename = file.path(save_dir, paste0(name_umap, ".pdf")),
  plot = p_umap,
  width = w_umap,
  height = h_umap,
  device = "pdf"
)


# 9. 参数记录模块

# 9.1 参数文件保存参数
name_params <- "umap_tsne_ellipse_parameters"

# 9.2 生成参数记录文本
param_text <- paste0(
  "本次分析参数总结：\n",
  "- 线条粗细(line_size)：", line_size, "\n",
  "- 填充透明度(bxalpha)：", bxalpha, "\n",
  "- 点大小(point_size)：", point_size, "\n",
  "- 线条类型(line_type)：", line_type_1, "\n",
  "- 标签大小(label_size)：", label_size, "\n",
  "- 是否显示置信椭圆(show_ellipse)：", show_ellipse, "\n",
  "- 分面列(facet_column)：", ifelse(facet_column == "", "不分面", facet_column), "\n",
  "- 分面列数(facet_ncol)：", facet_ncol, "\n",
  "- 分面横向间距(facet_panel_spacing_x)：", facet_panel_spacing_x, "\n",
  "- 分面纵向间距(facet_panel_spacing_y)：", facet_panel_spacing_y, "\n",
  "- 输出文件夹：", save_dir, "\n",
  "- tSNE 颜色集合：", palette_tsne_text, "\n",
  "- UMAP 颜色集合：", palette_umap_text, "\n",
  "\n写作提示词（自行组装语言，或借助AI组装）：\n",
  "1.在已有UMAP/tSNE降维结果基础上，按cellType进行分组可视化。\n",
  "2.可选择是否使用95%置信椭圆（stat_ellipse, level=0.95）辅助展示各细胞类型在低维空间的分布范围与聚集趋势。\n",
  "3.可通过 facet_column 参数按指定 meta.data 列进行分面展示，例如 group 或 orig.ident。\n",
  "4.通过标签中心点（各组坐标中位数）标注cellType，提升可读性。"
)

# 9.3 保存参数记录
writeLines(
  param_text,
  con = file.path(save_dir, paste0(name_params, ".txt"))
)


# 10. 完成提示
cat("\nUMAP 和 tSNE 绘图完成。\n")
cat("结果保存目录：", save_dir, "\n")
cat("是否显示置信椭圆：", show_ellipse, "\n")
cat("分面列：", ifelse(facet_column == "", "不分面", facet_column), "\n")
cat("tSNE 图已保存：", file.path(save_dir, paste0(name_tsne, ".pdf")), "\n")
cat("UMAP 图已保存：", file.path(save_dir, paste0(name_umap, ".pdf")), "\n")
cat("参数文件已保存：", file.path(save_dir, paste0(name_params, ".txt")), "\n")