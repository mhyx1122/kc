suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(scop)
})

# 1. 读取Spatial_Data并运行差异分析

# 1.1 参数设置

out_dir <- "差异基因可视化图"

# 可选：cellType、group、cluster
group_by <- "cellType"

# 是否只保留上调基因
only_pos <- FALSE

# 差异分析参数
p_val_threshold <- 0.05
logfc_threshold <- 0.1
min_pct <- 0.1

# 火山图使用的统计值，可选：p_val、p_val_adj
p_value_choice <- "p_val_adj"

# 需要额外标注的基因，多个基因用英文逗号分隔
features_label_input <- ""

# 1.2 创建输出文件夹

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# 1.3 检查并读取Spatial_Data对象

if (!exists("Spatial_Data", envir = .GlobalEnv)) {
  stop("全局环境中没有 Spatial_Data 对象。")
}

srt <- get("Spatial_Data", envir = .GlobalEnv)

if (!inherits(srt, "Seurat")) {
  stop("Spatial_Data 不是有效的Seurat对象。")
}

# 1.4 检查分组字段

if (!group_by %in% colnames(srt@meta.data)) {
  stop(
    paste0(
      "Spatial_Data@meta.data 中不存在分组列：",
      group_by
    )
  )
}

# 1.5 运行FindAllMarkers差异分析

markers <- FindAllMarkers(
  srt,
  group.by = group_by,
  only.pos = only_pos,
  min.pct = min_pct,
  logfc.threshold = logfc_threshold
)

if (nrow(markers) == 0) {
  stop("FindAllMarkers()没有获得差异基因，请检查分组和筛选参数。")
}

# 1.6 将cluster列改为group1

de <- markers

colnames(de)[colnames(de) == "cluster"] <- "group1"

# 2. 整理火山图数据

# 2.1 根据选择设置P值类型

de2 <- de

if (p_value_choice == "p_val") {
  colnames(de2)[colnames(de2) == "p_val_adj"] <- "p_val_adj2"
  colnames(de2)[colnames(de2) == "p_val"] <- "p_val_adj"
  ylab_text <- "-log10(p-val)"
} else if (p_value_choice == "p_val_adj") {
  ylab_text <- "-log10(p-adjust)"
} else {
  stop("p_value_choice只能设置为p_val或p_val_adj。")
}

# 2.2 将差异分析结果写入Seurat对象的tools槽

if (is.null(srt@tools$DEtest_cellType)) {
  srt@tools$DEtest_cellType <- list()
}

srt@tools$DEtest_cellType$AllMarkers_wilcox <- de2

# 2.3 解析需要标注的基因

features_label <- unlist(strsplit(features_label_input, ","))
features_label <- trimws(features_label)
features_label <- features_label[features_label != ""]
features_label <- unique(features_label)

if (length(features_label) == 0) {
  features_label <- NULL
}

# 3. 绘制火山图

# 3.1 参数设置

img_width <- 12
img_height <- 9
plot_name <- "volcano_plot"

# 3.2 绘制火山图

volcano_plot <- VolcanoPlot(
  srt = srt,
  DE_threshold = paste0(
    "avg_log2FC > ", logfc_threshold,
    " & p_val_adj < ", p_val_threshold
  ),
  features_label = features_label,
  group_by = group_by,
  x_metric = "avg_log2FC",
  ylab = ylab_text
)

print(volcano_plot)

# 3.3 保存火山图

plot_file <- file.path(
  out_dir,
  paste0(plot_name, ".pdf")
)

ggsave(
  filename = plot_file,
  plot = volcano_plot,
  device = "pdf",
  width = img_width,
  height = img_height
)

# 4. 保存参数记录

# 4.1 参数设置

name_params <- "volcano_parameters"

# 4.2 生成参数记录

param_text <- paste0(
  "本次分析参数总结：\n",
  "- group_by：", group_by, "\n",
  "- only.pos：", only_pos, "\n",
  "- p值阈值：", p_val_threshold, "\n",
  "- LogFC阈值：", logfc_threshold, "\n",
  "- min.pct：", min_pct, "\n",
  "- 统计值类型：", p_value_choice, "\n",
  "- 特定基因标注：",
  if (is.null(features_label)) {
    "无"
  } else {
    paste(features_label, collapse = ", ")
  },
  "\n",
  "- 差异基因结果行数：", nrow(de), "\n",
  "- 保存宽高：", img_width, " × ", img_height, "英寸\n",
  "- 图片文件：", basename(plot_file), "\n",
  "- 输出目录：", out_dir, "\n",
  "写作提示词（自行组装语言，或借助AI组装）：\n",
  "1. 基于指定分组变量进行单细胞差异分析，筛选具有统计学意义的差异表达基因。\n",
  "2. 结合火山图展示基因表达倍数变化与显著性分布特征，识别上调和下调的关键候选基因。\n",
  "3. 如设置特定基因标注，可进一步突出关注基因在整体差异表达格局中的位置。"
)

# 4.3 保存参数记录

param_file <- file.path(
  out_dir,
  paste0(name_params, ".txt")
)

writeLines(
  param_text,
  con = param_file
)