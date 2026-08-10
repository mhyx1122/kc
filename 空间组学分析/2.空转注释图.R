suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(SCpubr)
})

# 1. 设置公共输出文件夹

out_dir <- "手动注释"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# 2. 读取Marker基因表并整理基因分组

# 2.1 参数设置

cell_features_file <- "Cell Features.csv"

# 2.2 读取Marker基因表

if (!file.exists(cell_features_file)) {
  stop(paste0("未找到Marker基因文件：", cell_features_file))
}

cell_features <- read.csv(
  cell_features_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (ncol(cell_features) == 0) {
  stop("Marker基因表中没有可用的列。")
}

# 2.3 将每一列整理为一个Marker基因分组

marker_groups <- list()

for (group_name in names(cell_features)) {
  gene_vector <- as.character(cell_features[[group_name]])
  gene_vector <- gene_vector[!is.na(gene_vector)]
  gene_vector <- trimws(gene_vector)
  gene_vector <- gene_vector[gene_vector != ""]
  
  marker_groups[[group_name]] <- gene_vector
}

if (length(marker_groups) == 0) {
  stop("没有读取到有效的Marker基因分组。")
}

# 2.4 将Marker基因分组保存到全局环境

assign(
  "g8L5w3N2x7R4k1P9",
  marker_groups,
  envir = .GlobalEnv
)

# 2.5 读取全局环境中的Spatial_Data对象

if (!exists("Spatial_Data", envir = .GlobalEnv)) {
  stop("全局环境中不存在 Spatial_Data，请先完成前面的分析步骤。")
}

srt <- get("Spatial_Data", envir = .GlobalEnv)

if (!inherits(srt, "Seurat")) {
  stop("Spatial_Data 不是有效的Seurat对象。")
}

assign("Spatial_Data", srt, envir = .GlobalEnv)

# 2.6 生成Marker基因载入摘要

marker_load_summary <- paste0(
  "Marker基因组已载入。\n",
  "分组数量：", length(marker_groups), "\n",
  "各组基因数：\n",
  paste(
    paste0(
      " - ",
      names(marker_groups),
      ": ",
      lengths(marker_groups)
    ),
    collapse = "\n"
  )
)

# 2.7 记录Marker基因载入参数

param_step1_text <- paste0(
  "Step 1：载入Marker基因表并初始化\n",
  "运行流程说明：\n",
  "1. 使用read.csv()读取Marker基因CSV文件。\n",
  "2. 对CSV各列进行整理，每一列作为一个Marker基因分组。\n",
  "3. 从全局环境中读取Spatial_Data对象。\n",
  "4. 将Marker基因分组保存为全局变量g8L5w3N2x7R4k1P9。\n\n",
  "本次运行参数：\n",
  "- 文件路径：", cell_features_file, "\n",
  "- Marker分组数：", length(marker_groups), "\n",
  "- 各组名称：", paste(names(marker_groups), collapse = ", "), "\n"
)

# 3. 绘制第一次注释气泡图

# 3.1 参数设置

dot_scale <- 20
legend_framewidth <- 2
font_size <- 20

dotplot_colors <- c("#FFFFFF", "#FEE08B", "#D73027")

w_dotplot <- 30
h_dotplot <- 15
name_dotplot <- "1.第一次注释气泡图"

# 3.2 检查气泡图颜色参数

if (length(dotplot_colors) < 3) {
  stop("气泡图颜色至少需要3个，例如：#FFFFFF、#FEE08B、#D73027。")
}

# 3.3 检查do_DotPlot函数

if (!exists("do_DotPlot", mode = "function")) {
  stop("当前环境中未找到do_DotPlot()函数，请检查SCpubr包是否正确加载。")
}

# 3.4 绘制气泡图

dotplot_plot <- do_DotPlot(
  sample = srt,
  features = marker_groups,
  dot.scale = dot_scale,
  legend.framewidth = legend_framewidth,
  font.size = font_size
)

print(dotplot_plot)

# 3.5 保存气泡图

ggsave(
  filename = file.path(out_dir, paste0(name_dotplot, ".pdf")),
  plot = dotplot_plot,
  width = w_dotplot,
  height = h_dotplot,
  device = "pdf"
)

# 3.6 记录气泡图参数

param_step2_text <- paste0(
  "Step 2：第一次注释气泡图\n",
  "运行流程说明：\n",
  "1. 从全局Spatial_Data读取Seurat对象。\n",
  "2. 调用do_DotPlot()按照Marker基因分组绘制注释气泡图。\n",
  "3. 使用dot.scale、legend.framewidth和font.size控制图形样式。\n\n",
  "本次运行参数：\n",
  "- dot.scale：", dot_scale, "\n",
  "- legend.framewidth：", legend_framewidth, "\n",
  "- font.size：", font_size, "\n",
  "- colors：", paste(dotplot_colors, collapse = ","), "\n",
  "- 保存宽高：", w_dotplot, " × ", h_dotplot, "英寸\n",
  "- 文件名：", name_dotplot, ".pdf\n"
)

# 4. 绘制UMAP基因表达分布图

# 4.1 参数设置

marker_group_umap <- names(marker_groups)[1]

feature_pt_size_umap <- 0.1
feature_ncol_umap <- 2
feature_colors_umap <- c("lightgrey", "#FF0000")

w_feature_umap <- 12
h_feature_umap <- 12

name_feature_umap <- "2.基因分布UMAP"
name_all_umap_prefix <- "2.基因分布UMAP"

# 4.2 检查UMAP降维结果

if (!"umap" %in% Reductions(srt)) {
  stop("Spatial_Data中不存在UMAP降维结果。")
}

# 4.3 检查选择的Marker基因组

if (!marker_group_umap %in% names(marker_groups)) {
  stop(
    paste0(
      "Marker基因组不存在：",
      marker_group_umap,
      "。可选分组为：",
      paste(names(marker_groups), collapse = ", ")
    )
  )
}

feature_genes_umap <- marker_groups[[marker_group_umap]]

if (length(feature_genes_umap) == 0) {
  stop("当前选择的UMAP Marker基因组中没有有效基因。")
}

if (length(feature_colors_umap) < 2) {
  feature_colors_umap <- c("lightgrey", "#FF0000")
}

# 4.4 绘制当前Marker基因组的UMAP基因分布图

feature_umap_plot <- FeaturePlot(
  srt,
  features = feature_genes_umap,
  reduction = "umap",
  min.cutoff = NA,
  max.cutoff = NA,
  ncol = feature_ncol_umap,
  pt.size = feature_pt_size_umap,
  cols = feature_colors_umap[1:2]
)

print(feature_umap_plot)

# 4.5 保存当前Marker基因组的UMAP基因分布图

ggsave(
  filename = file.path(
    out_dir,
    paste0(name_feature_umap, "_", marker_group_umap, ".pdf")
  ),
  plot = feature_umap_plot,
  width = w_feature_umap,
  height = h_feature_umap,
  device = "pdf"
)

# 4.6 批量导出全部Marker基因组的UMAP基因分布图

for (group_name in names(marker_groups)) {
  current_genes <- marker_groups[[group_name]]
  
  current_umap_plot <- FeaturePlot(
    srt,
    features = current_genes,
    reduction = "umap",
    min.cutoff = NA,
    max.cutoff = NA,
    ncol = feature_ncol_umap,
    pt.size = feature_pt_size_umap,
    cols = feature_colors_umap[1:2]
  )
  
  ggsave(
    filename = file.path(
      out_dir,
      paste0(name_all_umap_prefix, "_", group_name, ".pdf")
    ),
    plot = current_umap_plot,
    width = w_feature_umap,
    height = h_feature_umap,
    device = "pdf"
  )
}

# 4.7 记录UMAP基因分布图参数

param_step3_text <- paste0(
  "Step 3：UMAP基因分布图\n",
  "运行流程说明：\n",
  "1. 从全局Spatial_Data读取Seurat对象。\n",
  "2. 根据选择的Marker基因组提取对应基因列表。\n",
  "3. 调用Seurat::FeaturePlot(reduction = umap)绘制基因表达分布图。\n",
  "4. 导出当前Marker基因组图形，并批量导出全部Marker组的UMAP图。\n\n",
  "本次运行参数：\n",
  "- marker_group：", marker_group_umap, "\n",
  "- pt.size：", feature_pt_size_umap, "\n",
  "- ncol：", feature_ncol_umap, "\n",
  "- colors：", paste(feature_colors_umap, collapse = ","), "\n",
  "- 保存宽高：", w_feature_umap, " × ", h_feature_umap, "英寸\n",
  "- 当前图文件名前缀：", name_feature_umap, "\n",
  "- 批量导出文件名前缀：", name_all_umap_prefix, "\n"
)

# 5. 绘制tSNE基因表达分布图

# 5.1 参数设置

marker_group_tsne <- names(marker_groups)[1]

feature_pt_size_tsne <- 0.1
feature_ncol_tsne <- 2
feature_colors_tsne <- c("lightgrey", "#FF0000")

w_feature_tsne <- 12
h_feature_tsne <- 12

name_feature_tsne <- "3.基因分布TSNE"
name_all_tsne_prefix <- "3.基因分布TSNE"

# 5.2 检查tSNE降维结果

if (!"tsne" %in% Reductions(srt)) {
  stop("Spatial_Data中不存在tSNE降维结果。")
}

# 5.3 检查选择的Marker基因组

if (!marker_group_tsne %in% names(marker_groups)) {
  stop(
    paste0(
      "Marker基因组不存在：",
      marker_group_tsne,
      "。可选分组为：",
      paste(names(marker_groups), collapse = ", ")
    )
  )
}

feature_genes_tsne <- marker_groups[[marker_group_tsne]]

if (length(feature_genes_tsne) == 0) {
  stop("当前选择的tSNE Marker基因组中没有有效基因。")
}

if (length(feature_colors_tsne) < 2) {
  feature_colors_tsne <- c("lightgrey", "#FF0000")
}

# 5.4 绘制当前Marker基因组的tSNE基因分布图

feature_tsne_plot <- FeaturePlot(
  srt,
  features = feature_genes_tsne,
  reduction = "tsne",
  min.cutoff = NA,
  max.cutoff = NA,
  ncol = feature_ncol_tsne,
  pt.size = feature_pt_size_tsne,
  cols = feature_colors_tsne[1:2]
)

print(feature_tsne_plot)

# 5.5 保存当前Marker基因组的tSNE基因分布图

ggsave(
  filename = file.path(
    out_dir,
    paste0(name_feature_tsne, "_", marker_group_tsne, ".pdf")
  ),
  plot = feature_tsne_plot,
  width = w_feature_tsne,
  height = h_feature_tsne,
  device = "pdf"
)

# 5.6 批量导出全部Marker基因组的tSNE基因分布图

for (group_name in names(marker_groups)) {
  current_genes <- marker_groups[[group_name]]
  
  current_tsne_plot <- FeaturePlot(
    srt,
    features = current_genes,
    reduction = "tsne",
    min.cutoff = NA,
    max.cutoff = NA,
    ncol = feature_ncol_tsne,
    pt.size = feature_pt_size_tsne,
    cols = feature_colors_tsne[1:2]
  )
  
  ggsave(
    filename = file.path(
      out_dir,
      paste0(name_all_tsne_prefix, "_", group_name, ".pdf")
    ),
    plot = current_tsne_plot,
    width = w_feature_tsne,
    height = h_feature_tsne,
    device = "pdf"
  )
}

# 5.7 记录tSNE基因分布图参数

param_step4_text <- paste0(
  "Step 4：tSNE基因分布图\n",
  "运行流程说明：\n",
  "1. 从全局Spatial_Data读取Seurat对象。\n",
  "2. 根据选择的Marker基因组提取对应基因列表。\n",
  "3. 调用Seurat::FeaturePlot(reduction = tsne)绘制基因表达分布图。\n",
  "4. 导出当前Marker基因组图形，并批量导出全部Marker组的tSNE图。\n\n",
  "本次运行参数：\n",
  "- marker_group：", marker_group_tsne, "\n",
  "- pt.size：", feature_pt_size_tsne, "\n",
  "- ncol：", feature_ncol_tsne, "\n",
  "- colors：", paste(feature_colors_tsne, collapse = ","), "\n",
  "- 保存宽高：", w_feature_tsne, " × ", h_feature_tsne, "英寸\n",
  "- 当前图文件名前缀：", name_feature_tsne, "\n",
  "- 批量导出文件名前缀：", name_all_tsne_prefix, "\n"
)

# 6. 保存全部参数记录

# 6.1 参数设置

name_params <- "MarkerGene_annotation_parameters"

# 6.2 合并各步骤参数

parameter_summary <- paste(
  param_step1_text,
  param_step2_text,
  param_step3_text,
  param_step4_text,
  sep = "\n\n------------------------------\n\n"
)

# 6.3 保存参数记录

writeLines(
  parameter_summary,
  con = file.path(out_dir, paste0(name_params, ".txt"))
)