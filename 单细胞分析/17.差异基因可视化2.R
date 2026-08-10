suppressPackageStartupMessages({
  library(ggplot2)
  library(Seurat)
  library(scop)
})

# 1. 检查全局环境中是否存在 seurat 对象
if (!exists("seurat", envir = .GlobalEnv)) {
  stop("全局环境中没有 seurat 对象。")
}

srt <- get("seurat", envir = .GlobalEnv)


# 2. 输出目录模块

# 2.1 输出目录参数
out_dir <- "7.差异基因可视化图"

# 2.2 创建输出目录
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}


# 3. 差异分析参数模块

# 3.1 差异比较依据
# 可选："cellType"、"group"、"cluster"
group_by <- "cellType"

# 3.2 FindAllMarkers 参数
only_pos <- FALSE
p_val_threshold <- 0.05
logfc_threshold <- 0.1
min_pct <- 0.1

# 3.3 火山图统计值选择
# 可选："p_val" 或 "p_val_adj"
p_value_choice <- "p_val_adj"

# 3.4 特定基因标注参数
features_label_text <- ""

features_label <- unlist(strsplit(features_label_text, ","))
features_label <- trimws(features_label)
features_label <- features_label[features_label != ""]

if (length(features_label) == 0) {
  features_label <- NULL
}

# 3.5 保存 PDF 参数
img_width <- 12
img_height <- 9
plot_name <- "volcano_plot"


# 4. 参数检查模块

# 4.1 检查分组列是否存在
if (!group_by %in% colnames(srt@meta.data)) {
  stop(paste0("seurat@meta.data 中不存在列：", group_by))
}


# 5. 差异分析模块

# 5.1 运行 FindAllMarkers
markers <- FindAllMarkers(
  srt,
  group.by = group_by,
  only.pos = only_pos,
  min.pct = min_pct,
  logfc.threshold = logfc_threshold
)

# 5.2 整理差异分析结果列名
de <- markers

colnames(de)[colnames(de) == "cluster"] <- "group1"


# 6. 火山图数据准备模块

# 6.1 根据 p_value_choice 设置火山图使用的统计值
de2 <- de

if (p_value_choice == "p_val") {
  colnames(de2)[colnames(de2) == "p_val_adj"] <- "p_val_adj2"
  colnames(de2)[colnames(de2) == "p_val"] <- "p_val_adj"
  ylab_text <- "-log10(p-val)"
} else if (p_value_choice == "p_val_adj") {
  ylab_text <- "-log10(p-adjust)"
} else {
  stop("p_value_choice 只能设置为 p_val 或 p_val_adj")
}

# 6.2 将差异结果写入 seurat@tools
tool_slot_name <- paste0("DEtest_", group_by)

if (is.null(srt@tools[[tool_slot_name]])) {
  srt@tools[[tool_slot_name]] <- list()
}

srt@tools[[tool_slot_name]]$AllMarkers_wilcox <- de2


# 7. 火山图绘制模块

# 7.1 生成火山图
p_volcano <- VolcanoPlot(
  srt = srt,
  DE_threshold = paste0(
    "avg_log2FC > ", logfc_threshold,
    " & ",
    "p_val_adj < ", p_val_threshold
  ),
  features_label = features_label,
  group_by = group_by,
  x_metric = "avg_log2FC",
  ylab = ylab_text
)

# 7.2 保存火山图
ggsave(
  filename = file.path(out_dir, paste0(plot_name, ".pdf")),
  plot = p_volcano,
  device = "pdf",
  width = img_width,
  height = img_height
)


# 8. 参数记录模块

# 8.1 参数文件保存参数
name_params <- "volcano_parameters"

# 8.2 生成参数记录文本
param_text <- paste0(
  "本次分析参数总结：\n",
  "- group_by：", group_by, "\n",
  "- only.pos：", only_pos, "\n",
  "- p值阈值：", p_val_threshold, "\n",
  "- LogFC阈值：", logfc_threshold, "\n",
  "- min.pct：", min_pct, "\n",
  "- 统计值类型：", p_value_choice, "\n",
  "- 特定基因标注：", if (is.null(features_label)) "无" else paste(features_label, collapse = ", "), "\n",
  "- 差异基因结果行数：", nrow(de), "\n",
  "写作提示词（自行组装语言，或借助AI组装）：\n",
  "1.基于指定分组变量进行单细胞差异分析，筛选具有统计学意义的差异表达基因。\n",
  "2.结合火山图展示基因表达倍数变化与显著性分布特征，识别上调和下调的关键候选基因。\n",
  "3.如设置特定基因标注，可进一步突出关注基因在整体差异表达格局中的位置。"
)

# 8.3 保存参数记录
writeLines(
  param_text,
  con = file.path(out_dir, paste0(name_params, ".txt"))
)


# 9. 写回 seurat 对象

seurat <- srt