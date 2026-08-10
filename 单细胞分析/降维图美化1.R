suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(scop)
})

# 1. 检查全局环境中是否存在 seurat 对象
if (!exists("seurat", envir = .GlobalEnv)) {
  stop("全局环境中没有 seurat 对象")
}

srt <- get("seurat", envir = .GlobalEnv)


# 2. 输出目录模块

# 2.1 输出目录参数
out_dir <- "4.3降维图美化"

# 2.2 创建输出目录
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}


# 3. CellDimPlot 绘图参数模块

# 3.1 分组和分面参数
group_by <- "cellType"

split_by <- "None"
# 可选：
# split_by <- "group"
# split_by <- "orig.ident"

stat_by <- "None"
# 可选：
# stat_by <- "Phase"
# stat_by <- "group"

# 3.2 降维方法参数
reduction <- "umap"
# 可选：
# reduction <- "tsne"

# 3.3 点参数
pt_size <- 2
pt_alpha <- 0.8

# 3.4 置信椭圆 / 标记参数
add_mark <- FALSE
mark_type <- "ellipse"
# 可选：
# mark_type <- "hull"
# mark_type <- "ellipse"
# mark_type <- "rect"
# mark_type <- "circle"

mark_alpha <- 0.15

# 3.5 密度线参数
add_density <- FALSE
density_filled <- FALSE

density_filled_palcolor_text <- "#000000,#FF0000,#00FF00"
density_filled_palcolor <- unlist(strsplit(density_filled_palcolor_text, ","))
density_filled_palcolor <- trimws(density_filled_palcolor)
density_filled_palcolor <- density_filled_palcolor[density_filled_palcolor != ""]

# 3.6 标签参数
label <- FALSE
label_repel <- TRUE

# 3.7 颜色参数
palette_text <- "#FBB4AE,#B3CDE3,#CCEBC5,#DECBE4,#FED9A6,#FFFFCC,#E5D8BD,#FDDAEC,#F2F2F2"

palcolor <- unlist(strsplit(palette_text, ","))
palcolor <- trimws(palcolor)
palcolor <- palcolor[palcolor != ""]

# 3.8 保存图片参数
w_pdf <- 8
h_pdf <- 6
name_pdf <- "CellDimPlot"

# 3.9 参数文件名
name_params <- "CellDimPlot_parameters"


# 4. 参数整理和检查模块

# 4.1 将 None 转换为 NULL
if (identical(split_by, "None")) {
  split_by <- NULL
}

if (identical(stat_by, "None")) {
  stat_by <- NULL
}

# 4.2 检查 group.by 是否存在
if (!(group_by %in% colnames(srt@meta.data))) {
  stop(paste0("seurat@meta.data 中不存在列：", group_by))
}

# 4.3 检查 split.by 是否存在
if (!is.null(split_by) && !(split_by %in% colnames(srt@meta.data))) {
  stop(paste0("seurat@meta.data 中不存在 split.by 列：", split_by))
}

# 4.4 检查 stat.by 是否存在
if (!is.null(stat_by) && !(stat_by %in% colnames(srt@meta.data))) {
  stop(paste0("seurat@meta.data 中不存在 stat.by 列：", stat_by))
}

# 4.5 检查 reduction 是否存在
if (!(reduction %in% Reductions(srt))) {
  stop(paste0("seurat 中不存在 reduction：", reduction))
}


# 5. CellDimPlot 绘图模块

# 5.1 生成 CellDimPlot 图
p <- CellDimPlot(
  srt = srt,
  group.by = group_by,
  split.by = split_by,
  reduction = reduction,
  stat.by = stat_by,
  palcolor = palcolor,
  mark_type = mark_type,
  mark_alpha = mark_alpha,
  add_mark = add_mark,
  pt.size = pt_size,
  pt.alpha = pt_alpha,
  add_density = add_density,
  density_filled = density_filled,
  density_filled_palcolor = density_filled_palcolor,
  label = label,
  label_repel = label_repel
)

# 5.2 保存 CellDimPlot 图
ggsave(
  filename = file.path(out_dir, paste0(name_pdf, ".pdf")),
  plot = p,
  width = w_pdf,
  height = h_pdf,
  device = "pdf"
)


# 6. 参数记录模块

# 6.1 生成参数记录文本
param_text <- paste0(
  "本次分析参数总结：\n",
  "- group.by：", group_by, "\n",
  "- split.by：", ifelse(is.null(split_by), "NULL", split_by), "\n",
  "- stat.by：", ifelse(is.null(stat_by), "NULL", stat_by), "\n",
  "- reduction：", reduction, "\n",
  "- pt.size：", pt_size, "\n",
  "- pt.alpha：", pt_alpha, "\n",
  "- add_mark：", add_mark, "\n",
  "- mark_type：", mark_type, "\n",
  "- mark_alpha：", mark_alpha, "\n",
  "- add_density：", add_density, "\n",
  "- density_filled：", density_filled, "\n",
  "- density_filled_palcolor：", density_filled_palcolor_text, "\n",
  "- label：", label, "\n",
  "- label_repel：", label_repel, "\n",
  "- palcolor：", palette_text, "\n",
  "\n写作提示词（自行组装语言，或借助AI组装）：\n",
  "1.基于指定降维结果绘制细胞二维可视化图（CellDimPlot），并按指定分组信息进行着色。\n",
  "2.可选叠加置信椭圆/密度信息，以展示不同类别在降维空间中的分布特征。\n",
  "3.根据需要进行分面展示与标签标注，便于对不同组/样本进行对比。"
)

# 6.2 保存参数记录
writeLines(
  param_text,
  con = file.path(out_dir, paste0(name_params, ".txt"))
)