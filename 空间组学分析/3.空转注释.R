suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(plyr)
})

# 1. 设置输出文件夹

out_dir <- "手动注释结果"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# 2. 读取空间转录组对象和Marker基因表

# 2.1 参数设置

spatial_object_name <- "Spatial_Data"
marker_table_name <- "g8L5w3N2x7R4k1P9"

# 2.2 读取Spatial_Data对象

if (!exists(spatial_object_name, envir = .GlobalEnv)) {
  stop("全局环境中不存在 Spatial_Data，请先完成前面的分析步骤。")
}

srt <- get(spatial_object_name, envir = .GlobalEnv)

if (!inherits(srt, "Seurat")) {
  stop("Spatial_Data 不是有效的 Seurat 对象。")
}

# 2.3 读取Marker基因表

if (!exists(marker_table_name, envir = .GlobalEnv)) {
  stop("全局环境中不存在 g8L5w3N2x7R4k1P9，请先完成前一步 marker 基因整理。")
}

cell_features <- get(marker_table_name, envir = .GlobalEnv)

# 2.4 将Marker基因表整理为基因列表

genes <- list()

for (col_name in names(cell_features)) {
  gene_vector <- as.character(cell_features[[col_name]])
  gene_vector <- gene_vector[!is.na(gene_vector)]
  gene_vector <- trimws(gene_vector)
  gene_vector <- gene_vector[gene_vector != ""]
  
  genes[[col_name]] <- gene_vector
}

# 3. 对cluster进行手动注释

# 3.1 参数设置

cluster_column <- "cluster"

# 按照cluster 0、cluster 1、cluster 2……的顺序填写注释名称。
# 注释名称数量必须等于最大cluster编号加1。
# 不同cluster可以填写相同的细胞类型名称。

annotations <- c(
  # "T cells",
  # "Macrophages",
  # "Fibroblasts",
  # "Endothelial cells"
)

# 3.2 检查cluster列

if (!cluster_column %in% colnames(srt@meta.data)) {
  stop(
    paste0(
      "Spatial_Data@meta.data 中不存在 ",
      cluster_column,
      " 列。"
    )
  )
}

# 3.3 按原始逻辑计算cluster数量

cluster_count <- max(
  as.numeric(as.character(srt@meta.data$cluster))
) + 1

# 3.4 检查手动注释信息

if (length(annotations) == 0) {
  stop(
    paste0(
      "请先在 annotations 中填写手动注释名称。",
      "当前对象共有 ",
      cluster_count,
      " 个cluster。"
    )
  )
}

annotations <- as.character(annotations)

if (length(annotations) != cluster_count) {
  stop(
    paste0(
      "手动注释名称数量与cluster数量不一致。",
      "当前有 ",
      cluster_count,
      " 个cluster，但填写了 ",
      length(annotations),
      " 个注释名称。"
    )
  )
}

if (any(is.na(annotations)) || any(trimws(annotations) == "")) {
  stop("annotations 中存在空值，请为每个cluster填写注释名称。")
}

# 3.5 删除已有cellType并恢复cluster身份

if ("cellType" %in% colnames(srt@meta.data)) {
  srt@meta.data$cellType <- NULL
}

Idents(srt) <- cluster_column

Idents(srt) <- factor(
  Idents(srt),
  levels = sort(as.numeric(levels(Idents(srt))))
)

# 3.6 将cluster映射为cellType

seuratidens <- plyr::mapvalues(
  Idents(srt),
  from = levels(Idents(srt)),
  to = annotations
)

Idents(srt) <- seuratidens
srt$cellType <- Idents(srt)
Idents(srt) <- srt$cellType

# 3.7 更新全局环境中的Spatial_Data

assign(
  spatial_object_name,
  srt,
  envir = .GlobalEnv
)

# 3.8 生成手动注释结果摘要

celltype_table <- table(srt$cellType)

annotation_summary <- paste0(
  "手动注释已完成。\n",
  "cellType数量：", length(celltype_table), "\n",
  "各cellType的spot/cell数量：\n",
  paste(
    paste0(
      " - ",
      names(celltype_table),
      ": ",
      as.integer(celltype_table)
    ),
    collapse = "\n"
  )
)

# 3.9 记录手动注释参数

param_step1_text <- paste0(
  "Step 1：cluster手动注释\n",
  "运行流程说明：\n",
  "1. 从全局环境读取 Spatial_Data 对象。\n",
  "2. 从全局环境读取 marker 基因表 g8L5w3N2x7R4k1P9。\n",
  "3. 根据最大 cluster 编号加1计算 cluster 数量。\n",
  "4. 根据填写的注释名称，使用 plyr::mapvalues() 将 cluster 映射为 cellType。\n",
  "5. 将新的 cellType 写入 Seurat 对象 meta.data，并更新全局 Spatial_Data。\n",
  "6. 不同 cluster 可以被注释为相同的 cellType。\n\n",
  "本次运行参数：\n",
  "- cluster列：", cluster_column, "\n",
  "- cluster数量：", cluster_count, "\n",
  "- 注释信息：\n",
  paste0(
    "  cluster ",
    0:(length(annotations) - 1),
    " -> ",
    annotations,
    collapse = "\n"
  ),
  "\n"
)

# 4. 绘制第一次注释气泡图

# 4.1 参数设置

dotplot_width <- 30
dotplot_height <- 15
dot_scale <- 12
legend_framewidth <- 2
font_size <- 20
name_dotplot <- "1.第一次注释气泡结果图"

# 4.2 检查do_DotPlot函数

if (!exists("do_DotPlot", mode = "function")) {
  stop(
    paste0(
      "当前R环境中没有找到 do_DotPlot() 函数。",
      "请先加载包含 do_DotPlot() 的R包或函数代码。"
    )
  )
}

# 4.3 绘制第一次注释气泡图

dotplot_plot <- do_DotPlot(
  sample = srt,
  features = genes,
  dot.scale = dot_scale,
  legend.framewidth = legend_framewidth,
  font.size = font_size
)

print(dotplot_plot)

# 4.4 保存第一次注释气泡图

ggsave(
  filename = file.path(out_dir, paste0(name_dotplot, ".pdf")),
  plot = dotplot_plot,
  width = dotplot_width,
  height = dotplot_height,
  device = "pdf"
)

# 4.5 记录第一次注释气泡图参数

param_step2_text <- paste0(
  "Step 2：第一次注释气泡图\n",
  "运行流程说明：\n",
  "1. 将全局 marker 基因表整理为带名称的基因列表。\n",
  "2. 调用 do_DotPlot() 绘制第一次注释气泡图。\n",
  "3. 使用 dot.scale、legend.framewidth 和 font.size 控制图形样式。\n\n",
  "本次运行参数：\n",
  "- dot.scale：", dot_scale, "\n",
  "- legend.framewidth：", legend_framewidth, "\n",
  "- font.size：", font_size, "\n",
  "- 保存宽高：", dotplot_width, " × ", dotplot_height, " 英寸\n",
  "- 文件名：", name_dotplot, ".pdf\n"
)

# 5. 绘制分类UMAP图

# 5.1 参数设置

umap_width <- 7.5
umap_height <- 5.5
umap_label_size <- 3.5
umap_pt_size <- 0.5

umap_colors <- c(
  "#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF",
  "#91D1C2FF", "#7E6148FF", "#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF",
  "#7AA6DCFF", "#003C67FF", "#8F7700FF", "#3B3B3BFF", "#A73030FF", "#374E55FF",
  "#DF8F44FF", "#00A1D5FF", "#B24745FF", "#79AF97FF", "#6A6599FF", "#80796BFF",
  "#BC3C29FF", "#0072B5FF", "#E18727FF", "#20854EFF", "#7876B1FF", "#6F99ADFF",
  "#FFDC91FF", "#EE4C97FF"
)

name_umap <- "2.分类UMAP"

# 5.2 检查UMAP降维结果

if (!"umap" %in% Reductions(srt)) {
  stop("Spatial_Data 中不存在 umap 降维结果，无法绘制分类UMAP。")
}

# 5.3 检查UMAP颜色数量

if (length(umap_colors) < length(levels(Idents(srt)))) {
  stop("UMAP颜色数量不足以覆盖所有cellType。")
}

# 5.4 绘制分类UMAP

umap_plot <- DimPlot(
  srt,
  reduction = "umap",
  label = TRUE,
  cols = umap_colors,
  label.size = umap_label_size,
  pt.size = umap_pt_size,
  raster = FALSE
) +
  theme_classic() +
  theme(
    panel.border = element_rect(
      fill = NA,
      color = "black",
      size = 0.5,
      linetype = "solid"
    ),
    legend.position = "right"
  )

print(umap_plot)

# 5.5 保存分类UMAP

ggsave(
  filename = file.path(out_dir, paste0(name_umap, ".pdf")),
  plot = umap_plot,
  width = umap_width,
  height = umap_height,
  device = "pdf"
)

# 5.6 记录分类UMAP参数

param_step3_text <- paste0(
  "Step 3：分类UMAP\n",
  "运行流程说明：\n",
  "1. 使用 Seurat::DimPlot(reduction = umap) 绘制注释后的分类UMAP。\n",
  "2. 使用 label.size、pt.size 和颜色集合控制图形样式。\n",
  "3. 图中每个类别对应手动注释后的 cellType。\n\n",
  "本次运行参数：\n",
  "- label.size：", umap_label_size, "\n",
  "- pt.size：", umap_pt_size, "\n",
  "- colors：", paste(umap_colors, collapse = ","), "\n",
  "- 保存宽高：", umap_width, " × ", umap_height, " 英寸\n",
  "- 文件名：", name_umap, ".pdf\n"
)

# 6. 绘制分类tSNE图

# 6.1 参数设置

tsne_width <- 7.5
tsne_height <- 5.5
tsne_label_size <- 3.5
tsne_pt_size <- 0.5

tsne_colors <- c(
  "#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF",
  "#91D1C2FF", "#7E6148FF", "#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF",
  "#7AA6DCFF", "#003C67FF", "#8F7700FF", "#3B3B3BFF", "#A73030FF", "#374E55FF",
  "#DF8F44FF", "#00A1D5FF", "#B24745FF", "#79AF97FF", "#6A6599FF", "#80796BFF",
  "#BC3C29FF", "#0072B5FF", "#E18727FF", "#20854EFF", "#7876B1FF", "#6F99ADFF",
  "#FFDC91FF", "#EE4C97FF"
)

name_tsne <- "3.分类TSNE"

# 6.2 检查tSNE降维结果

if (!"tsne" %in% Reductions(srt)) {
  stop("Spatial_Data 中不存在 tsne 降维结果，无法绘制分类tSNE。")
}

# 6.3 检查tSNE颜色数量

if (length(tsne_colors) < length(levels(Idents(srt)))) {
  stop("tSNE颜色数量不足以覆盖所有cellType。")
}

# 6.4 绘制分类tSNE

tsne_plot <- DimPlot(
  srt,
  reduction = "tsne",
  label = TRUE,
  cols = tsne_colors,
  label.size = tsne_label_size,
  pt.size = tsne_pt_size,
  raster = FALSE
) +
  theme_classic() +
  theme(
    panel.border = element_rect(
      fill = NA,
      color = "black",
      size = 0.5,
      linetype = "solid"
    ),
    legend.position = "right"
  )

print(tsne_plot)

# 6.5 保存分类tSNE

ggsave(
  filename = file.path(out_dir, paste0(name_tsne, ".pdf")),
  plot = tsne_plot,
  width = tsne_width,
  height = tsne_height,
  device = "pdf"
)

# 6.6 记录分类tSNE参数

param_step4_text <- paste0(
  "Step 4：分类tSNE\n",
  "运行流程说明：\n",
  "1. 使用 Seurat::DimPlot(reduction = tsne) 绘制注释后的分类tSNE。\n",
  "2. 使用 label.size、pt.size 和颜色集合控制图形样式。\n",
  "3. 图中每个类别对应手动注释后的 cellType。\n\n",
  "本次运行参数：\n",
  "- label.size：", tsne_label_size, "\n",
  "- pt.size：", tsne_pt_size, "\n",
  "- colors：", paste(tsne_colors, collapse = ","), "\n",
  "- 保存宽高：", tsne_width, " × ", tsne_height, " 英寸\n",
  "- 文件名：", name_tsne, ".pdf\n"
)

# 7. 绘制空间定位图

# 7.1 参数设置

spatial_width <- 7.5
spatial_height <- 5.5
spatial_label_size <- 5
spatial_pt_size <- 1

spatial_colors <- c(
  "#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF",
  "#91D1C2FF", "#7E6148FF", "#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF",
  "#7AA6DCFF", "#003C67FF", "#8F7700FF", "#3B3B3BFF", "#A73030FF", "#374E55FF",
  "#DF8F44FF", "#00A1D5FF", "#B24745FF", "#79AF97FF", "#6A6599FF", "#80796BFF",
  "#BC3C29FF", "#0072B5FF", "#E18727FF", "#20854EFF", "#7876B1FF", "#6F99ADFF",
  "#FFDC91FF", "#EE4C97FF"
)

name_spatial <- "4.空间定位图"

# 7.2 检查空间图像信息

if (length(Images(srt)) == 0) {
  stop("Spatial_Data 中不存在空间图像信息，无法绘制SpatialPlot。")
}

# 7.3 整理空间定位图颜色

num_clusters <- length(unique(srt$cellType))
pal_use <- spatial_colors[1:num_clusters]

if (length(pal_use) < num_clusters) {
  stop("空间定位图颜色数量不足以覆盖所有唯一的cellType值。")
}

cluster_colors <- setNames(
  pal_use,
  levels(srt$cellType)
)

# 7.4 绘制空间定位图

spatial_plot <- suppressWarnings(
  SpatialPlot(
    srt,
    label = FALSE,
    label.size = spatial_label_size,
    cols = cluster_colors,
    pt.size.factor = spatial_pt_size
  )
)

print(spatial_plot)

# 7.5 保存空间定位图

ggsave(
  filename = file.path(out_dir, paste0(name_spatial, ".pdf")),
  plot = spatial_plot,
  width = spatial_width,
  height = spatial_height,
  device = "pdf"
)

# 7.6 记录空间定位图参数

param_step5_text <- paste0(
  "Step 5：空间定位图\n",
  "运行流程说明：\n",
  "1. 依据 cellType 设置空间类别颜色。\n",
  "2. 使用 Seurat::SpatialPlot() 绘制注释后的空间定位图。\n",
  "3. 自动检查颜色数量是否足够覆盖所有 cellType。\n\n",
  "本次运行参数：\n",
  "- label.size：", spatial_label_size, "\n",
  "- pt.size.factor：", spatial_pt_size, "\n",
  "- colors：", paste(spatial_colors, collapse = ","), "\n",
  "- 保存宽高：", spatial_width, " × ", spatial_height, " 英寸\n",
  "- 文件名：", name_spatial, ".pdf\n"
)

# 8. 生成二次注释可选项

# 8.1 按cellType拆分Seurat对象

aff <- SplitObject(
  srt,
  split.by = "cellType"
)

# 8.2 获取可以用于二次注释的cellType名称

subcluster_choices <- names(aff)

# 9. 保存全部参数记录

# 9.1 参数设置

name_params <- "Cell_RemarkSpatial_parameters"

# 9.2 合并各步骤参数记录

parameter_summary <- paste(
  param_step1_text,
  param_step2_text,
  param_step3_text,
  param_step4_text,
  param_step5_text,
  sep = "\n\n------------------------------\n\n"
)

# 9.3 保存参数记录

writeLines(
  parameter_summary,
  con = file.path(out_dir, paste0(name_params, ".txt"))
)