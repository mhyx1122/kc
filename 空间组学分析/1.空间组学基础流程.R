suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(clustree)
  library(patchwork)
})

# 1. 读取10X空间转录组数据并绘制预览图

# 1.1 参数设置

data_dir <- "10X_Spatial"
filename <- "matrix.h5"
slice_name <- "cgxr410"
sample_name <- "SAMPLE1"
step1_feature <- "nCount_Spatial"

step1_output_dir <- "1.空转预览图"
step1_plot_name <- "SpatialFeaturePlot"
step1_plot_width <- 8
step1_plot_height <- 6

# 1.2 创建输出文件夹

if (!dir.exists(step1_output_dir)) {
  dir.create(step1_output_dir, recursive = TRUE)
}

# 1.3 读取空间转录组数据

Spatial_Data <- Load10X_Spatial(
  data.dir = data_dir,
  filename = filename,
  slice = slice_name
)

# 1.4 写入样本名称

Spatial_Data$orig.ident <- sample_name

# 1.5 绘制空间特征预览图

step1_plot <- suppressWarnings(
  SpatialFeaturePlot(
    Spatial_Data,
    features = step1_feature
  )
)

print(step1_plot)

ggsave(
  filename = file.path(step1_output_dir, paste0(step1_plot_name, ".pdf")),
  plot = step1_plot,
  width = step1_plot_width,
  height = step1_plot_height,
  device = "pdf"
)

# 1.6 记录本步骤参数

param_step1_text <- paste0(
  "Step 1：读取空间数据\n",
  "运行流程说明：\n",
  "1. 使用 Seurat::Load10X_Spatial() 读取 10X 空间转录组数据。\n",
  "2. 写入 orig.ident 样本标签。\n",
  "3. 将对象保存为 Spatial_Data。\n",
  "4. 使用 Seurat::SpatialFeaturePlot() 绘制空间特征图。\n\n",
  "本次运行参数：\n",
  "- data.dir: ", data_dir, "\n",
  "- filename: ", filename, "\n",
  "- slice: ", slice_name, "\n",
  "- sample_name: ", sample_name, "\n",
  "- feature: ", step1_feature, "\n"
)

# 2. 使用SCTransform对Spatial assay进行标准化

# 2.1 参数设置

sct_variable_features_n <- 3000
sct_vst_flavor <- "v2"
sct_vars_to_regress <- character(0)
sct_return_only_var_genes <- FALSE

step2_summary_name <- "SCTransform_summary"

# 2.2 整理需要回归的变量

vars_regress <- sct_vars_to_regress

if (length(vars_regress) == 0) {
  vars_regress <- NULL
}

# 2.3 运行SCTransform

Spatial_Data <- SCTransform(
  Spatial_Data,
  assay = "Spatial",
  variable.features.n = sct_variable_features_n,
  vst.flavor = sct_vst_flavor,
  vars.to.regress = vars_regress,
  return.only.var.genes = sct_return_only_var_genes,
  verbose = FALSE
)

# 2.4 保存SCTransform运行摘要

step2_summary_text <- paste0(
  "SCTransform 已完成。\n",
  "当前默认 assay: ", DefaultAssay(Spatial_Data), "\n",
  "对象中 assays: ", paste(Assays(Spatial_Data), collapse = ", "), "\n",
  "当前细胞数: ", ncol(Spatial_Data), "\n",
  "当前基因数: ", nrow(Spatial_Data)
)

writeLines(
  step2_summary_text,
  con = file.path(step1_output_dir, paste0(step2_summary_name, ".txt"))
)

# 2.5 记录本步骤参数

param_step2_text <- paste0(
  "Step 2：SCTransform标准化\n",
  "运行流程说明：\n",
  "1. 使用 Seurat::SCTransform() 对 Spatial assay 进行标准化。\n",
  "2. 将标准化后的对象继续保存为 Spatial_Data。\n\n",
  "本次运行参数：\n",
  "- assay: Spatial\n",
  "- variable.features.n: ", sct_variable_features_n, "\n",
  "- vst.flavor: ", sct_vst_flavor, "\n",
  "- vars.to.regress: ",
  ifelse(
    length(sct_vars_to_regress) == 0,
    "NULL",
    paste(sct_vars_to_regress, collapse = ", ")
  ),
  "\n",
  "- return.only.var.genes: ", sct_return_only_var_genes, "\n"
)

# 3. 主成分分析与主成分选择

# 3.1 参数设置

pc_reduction <- "pca"
num_replicate <- 0
dims_start <- 1
dims_end <- 20
cells_to_use <- 1000
elbow_dims <- 30

step3_output_dir <- "2.确定主成分的个数"

pc_heatmap_name <- "1.PC热图"
pc_heatmap_width <- 7.5
pc_heatmap_height <- 9

jack_plot_name <- "2.jackstrawplot"
jack_plot_width <- 7.5
jack_plot_height <- 5.5

elbow_plot_name <- "3.ElbowPlot"
elbow_plot_width <- 5
elbow_plot_height <- 4

# 3.2 创建输出文件夹

if (!dir.exists(step3_output_dir)) {
  dir.create(step3_output_dir, recursive = TRUE)
}

# 3.3 运行PCA

Spatial_Data <- RunPCA(
  Spatial_Data,
  assay = "SCT",
  verbose = FALSE
)

# 3.4 检查需要使用的降维结果

if (!pc_reduction %in% Reductions(Spatial_Data)) {
  stop(
    paste0(
      "对象中不存在 reduction：",
      pc_reduction,
      "。原代码没有运行 Harmony，若使用 harmony，需先自行生成 harmony reduction。"
    )
  )
}

# 3.5 绘制主成分热图

pc_heatmap <- suppressWarnings(
  DimHeatmap(
    Spatial_Data,
    nfeatures = 10,
    dims = dims_start:dims_end,
    cells = cells_to_use,
    balanced = TRUE,
    reduction = pc_reduction,
    fast = FALSE
  )
)

print(pc_heatmap)

ggsave(
  filename = file.path(step3_output_dir, paste0(pc_heatmap_name, ".pdf")),
  plot = pc_heatmap,
  width = pc_heatmap_width,
  height = pc_heatmap_height,
  device = "pdf"
)

# 3.6 绘制ElbowPlot

elbow_plot <- suppressWarnings(
  ElbowPlot(
    Spatial_Data,
    ndims = elbow_dims,
    reduction = pc_reduction
  )
)

print(elbow_plot)

ggsave(
  filename = file.path(step3_output_dir, paste0(elbow_plot_name, ".pdf")),
  plot = elbow_plot,
  width = elbow_plot_width,
  height = elbow_plot_height,
  device = "pdf"
)

# 3.7 根据参数决定是否运行JackStraw

if (num_replicate > 0) {
  Spatial_Data <- JackStraw(
    Spatial_Data,
    num.replicate = num_replicate
  )
  
  Spatial_Data <- ScoreJackStraw(
    Spatial_Data,
    dims = dims_start:dims_end
  )
  
  jack_plot <- suppressWarnings(
    JackStrawPlot(
      Spatial_Data,
      dims = dims_start:dims_end
    )
  )
} else {
  jack_plot <- ggplot() +
    theme_void() +
    ggtitle("num.replicate = 0，未运行 JackStraw")
}

print(jack_plot)

ggsave(
  filename = file.path(step3_output_dir, paste0(jack_plot_name, ".pdf")),
  plot = jack_plot,
  width = jack_plot_width,
  height = jack_plot_height,
  device = "pdf"
)

# 3.8 记录本步骤参数

param_step3_text <- paste0(
  "Step 3：主成分选择\n",
  "运行流程说明：\n",
  "1. 使用 Seurat::RunPCA() 计算主成分。\n",
  "2. 使用 Seurat::DimHeatmap() 绘制主成分热图。\n",
  "3. 若 num.replicate > 0，则使用 JackStraw()、ScoreJackStraw()、JackStrawPlot() 评估主成分显著性。\n",
  "4. 使用 ElbowPlot() 绘制拐点图。\n",
  "5. 所有处理结果继续保存于 Spatial_Data。\n\n",
  "本次运行参数：\n",
  "- reduction: ", pc_reduction, "\n",
  "- num.replicate: ", num_replicate, "\n",
  "- dims: ", dims_start, ":", dims_end, "\n",
  "- cells_to_use: ", cells_to_use, "\n",
  "- elbow_dims: ", elbow_dims, "\n"
)

# 4. 多分辨率聚类并绘制clustree

# 4.1 参数设置

resolution_1 <- 0.01
resolution_2 <- 0.05
resolution_3 <- 0.1
resolution_4 <- 0.2
resolution_5 <- 0.3
resolution_6 <- 0.4
resolution_7 <- 0.5
resolution_8 <- 0.6
resolution_9 <- 0.7
resolution_10 <- 0.8

clustree_reduction <- "pca"
clustree_seuratPC <- 20

step4_output_dir <- "3.选择分辨率"
clustree_plot_name <- "挑选分辨率"
clustree_plot_width <- 12
clustree_plot_height <- 10

# 4.2 创建输出文件夹

if (!dir.exists(step4_output_dir)) {
  dir.create(step4_output_dir, recursive = TRUE)
}

# 4.3 检查需要使用的降维结果

if (!clustree_reduction %in% Reductions(Spatial_Data)) {
  stop(paste0("对象中不存在 reduction：", clustree_reduction))
}

# 4.4 整理分辨率集合

res_vec <- c(
  resolution_1,
  resolution_2,
  resolution_3,
  resolution_4,
  resolution_5,
  resolution_6,
  resolution_7,
  resolution_8,
  resolution_9,
  resolution_10
)

# 4.5 构建邻接图

Spatial_Data <- FindNeighbors(
  Spatial_Data,
  dims = 1:clustree_seuratPC,
  reduction = clustree_reduction
)

# 4.6 在多个分辨率下运行聚类

for (r in res_vec) {
  Spatial_Data <- FindClusters(
    Spatial_Data,
    graph.name = "SCT_snn",
    resolution = r,
    algorithm = 1
  )
}

# 4.7 绘制clustree

clustree_plot <- clustree(
  Spatial_Data@meta.data,
  prefix = "SCT_snn_res."
)

print(clustree_plot)

ggsave(
  filename = file.path(step4_output_dir, paste0(clustree_plot_name, ".pdf")),
  plot = clustree_plot,
  width = clustree_plot_width,
  height = clustree_plot_height,
  device = "pdf"
)

# 4.8 记录本步骤参数

param_step4_text <- paste0(
  "Step 4：分辨率选择\n",
  "运行流程说明：\n",
  "1. 使用 Seurat::FindNeighbors() 构建邻接图。\n",
  "2. 在多个 resolution 下重复运行 Seurat::FindClusters()。\n",
  "3. 使用 clustree::clustree() 可视化不同分辨率下的聚类层级关系。\n",
  "4. 所有处理结果继续保存于 Spatial_Data。\n\n",
  "本次运行参数：\n",
  "- reduction: ", clustree_reduction, "\n",
  "- seuratPC: ", clustree_seuratPC, "\n",
  "- 分辨率集合: ", paste(res_vec, collapse = ", "), "\n"
)

# 5. 正式聚类、Marker分析、聚类热图、UMAP和tSNE

# 5.1 参数设置

cluster_resolution <- 0.1
only_pos <- TRUE
min_pct <- 0.25
logfc_threshold <- 0.25
top_genes_per_cluster <- 5
cluster_reduction <- "pca"
cluster_seuratPC <- 20

step5_output_dir <- "4.分群聚类"

heat_group_colors <- c(
  "#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF",
  "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#7E6148FF",
  "#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF",
  "#7AA6DCFF", "#003C67FF", "#8F7700FF", "#3B3B3BFF",
  "#A73030FF", "#374E55FF", "#DF8F44FF", "#00A1D5FF",
  "#B24745FF", "#79AF97FF", "#6A6599FF", "#80796BFF",
  "#BC3C29FF", "#0072B5FF", "#E18727FF", "#20854EFF",
  "#7876B1FF", "#6F99ADFF", "#FFDC91FF", "#EE4C97FF"
)

heat_low <- "#0099CC"
heat_mid <- "#FFFFFF"
heat_high <- "#CC0033"
heat_name <- "Z-score"

cluster_heatmap_name <- "1.聚类热图"
cluster_heatmap_width <- 22
cluster_heatmap_height <- 16

pt_size_umap <- 0.5
palette_umap <- heat_group_colors
umap_plot_name <- "2.聚类后UMAP(按细胞群给色)"
umap_plot_width <- 6.5
umap_plot_height <- 5.5

pt_size_tsne <- 0.5
palette_tsne <- heat_group_colors
tsne_plot_name <- "3.聚类后TSNE(按细胞群给色)"
tsne_plot_width <- 6.5
tsne_plot_height <- 5.5

# 5.2 创建输出文件夹

if (!dir.exists(step5_output_dir)) {
  dir.create(step5_output_dir, recursive = TRUE)
}

# 5.3 检查需要使用的降维结果

if (!cluster_reduction %in% Reductions(Spatial_Data)) {
  stop(paste0("对象中不存在 reduction：", cluster_reduction))
}

# 5.4 按选定分辨率完成正式聚类

Spatial_Data <- FindClusters(
  Spatial_Data,
  resolution = cluster_resolution
)

Spatial_Data@meta.data$cluster <- as.character(Idents(Spatial_Data))

# 5.5 识别各cluster的Marker基因

markers_all <- FindAllMarkers(
  Spatial_Data,
  only.pos = only_pos,
  min.pct = min_pct,
  logfc.threshold = logfc_threshold
)

if (nrow(markers_all) == 0) {
  stop("FindAllMarkers 未返回结果，请调整参数。")
}

write.csv(
  markers_all,
  file = file.path(step5_output_dir, "每个聚类的marker基因.csv"),
  row.names = FALSE
)

# 5.6 按cluster拆分并保存每个cluster的前100个Marker基因

split_tables <- markers_all %>%
  group_by(cluster) %>%
  arrange(cluster, desc(avg_log2FC)) %>%
  group_split()

all_genes <- data.frame()

for (i in seq_along(split_tables)) {
  cluster_data <- split_tables[[i]]
  top100 <- cluster_data[1:min(100, nrow(cluster_data)), ]
  all_genes <- rbind(all_genes, top100)
  
  write.csv(
    top100,
    file.path(
      step5_output_dir,
      paste0("cluster_", i - 1, "_sorted.csv")
    ),
    row.names = FALSE
  )
}

write.csv(
  all_genes,
  file.path(step5_output_dir, "所有聚类_top100合并表.csv"),
  row.names = FALSE
)

# 5.7 提取每个cluster的Top Marker基因

top_markers <- markers_all %>%
  group_by(cluster) %>%
  top_n(
    n = top_genes_per_cluster,
    wt = avg_log2FC
  )

# 5.8 运行UMAP和tSNE

Spatial_Data <- RunUMAP(
  Spatial_Data,
  dims = 1:cluster_seuratPC,
  reduction = cluster_reduction,
  verbose = FALSE
)

Spatial_Data <- RunTSNE(
  Spatial_Data,
  dims = 1:cluster_seuratPC,
  reduction = cluster_reduction,
  check_duplicates = FALSE,
  verbose = FALSE
)

# 5.9 绘制聚类热图

cluster_heatmap <- suppressWarnings(
  DoHeatmap(
    Spatial_Data,
    features = top_markers$gene,
    group.colors = heat_group_colors
  )
) +
  scale_fill_gradient2(
    low = heat_low,
    mid = heat_mid,
    high = heat_high,
    name = heat_name
  )

print(cluster_heatmap)

ggsave(
  filename = file.path(
    step5_output_dir,
    paste0(cluster_heatmap_name, ".pdf")
  ),
  plot = cluster_heatmap,
  width = cluster_heatmap_width,
  height = cluster_heatmap_height,
  device = "pdf"
)

# 5.10 绘制UMAP图

umap_plot <- suppressWarnings(
  DimPlot(
    Spatial_Data,
    reduction = "umap",
    label = TRUE,
    label.size = 3.5,
    pt.size = pt_size_umap,
    cols = palette_umap,
    raster = FALSE
  ) +
    theme_classic() +
    theme(
      panel.border = element_rect(
        fill = NA,
        color = "black",
        size = 0.5
      ),
      legend.position = "right"
    )
)

print(umap_plot)

ggsave(
  filename = file.path(
    step5_output_dir,
    paste0(umap_plot_name, ".pdf")
  ),
  plot = umap_plot,
  width = umap_plot_width,
  height = umap_plot_height,
  device = "pdf"
)

# 5.11 绘制tSNE图

tsne_plot <- suppressWarnings(
  DimPlot(
    Spatial_Data,
    reduction = "tsne",
    label = TRUE,
    label.size = 3.5,
    pt.size = pt_size_tsne,
    cols = palette_tsne,
    raster = FALSE
  ) +
    theme_classic() +
    theme(
      panel.border = element_rect(
        fill = NA,
        color = "black",
        size = 0.5
      ),
      legend.position = "right"
    )
)

print(tsne_plot)

ggsave(
  filename = file.path(
    step5_output_dir,
    paste0(tsne_plot_name, ".pdf")
  ),
  plot = tsne_plot,
  width = tsne_plot_width,
  height = tsne_plot_height,
  device = "pdf"
)

# 5.12 记录本步骤参数

param_step5_text <- paste0(
  "Step 5：聚类与 marker 分析\n",
  "运行流程说明：\n",
  "1. 使用 Seurat::FindClusters() 完成正式聚类。\n",
  "2. 使用 Seurat::FindAllMarkers() 识别各 cluster 的差异基因。\n",
  "3. 写出 marker 基因表及各 cluster 的 top100 表格。\n",
  "4. 提取每个 cluster 的 top marker 用于 DoHeatmap() 绘图。\n",
  "5. 使用 RunUMAP() 与 RunTSNE() 进行降维。\n",
  "6. 所有处理结果继续保存于 Spatial_Data。\n\n",
  "本次运行参数：\n",
  "- cluster resolution: ", cluster_resolution, "\n",
  "- only.pos: ", only_pos, "\n",
  "- min.pct: ", min_pct, "\n",
  "- logfc.threshold: ", logfc_threshold, "\n",
  "- top_genes_per_cluster: ", top_genes_per_cluster, "\n",
  "- reduction: ", cluster_reduction, "\n",
  "- cluster_seuratPC: ", cluster_seuratPC, "\n"
)

# 6. 绘制空间定位图和空间图与UMAP组合图

# 6.1 参数设置

spatial_label_size <- 5

palette_spatial <- c(
  "#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF",
  "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#7E6148FF",
  "#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF",
  "#7AA6DCFF", "#003C67FF", "#8F7700FF", "#3B3B3BFF",
  "#A73030FF", "#374E55FF", "#DF8F44FF", "#00A1D5FF",
  "#B24745FF", "#79AF97FF", "#6A6599FF", "#80796BFF",
  "#BC3C29FF", "#0072B5FF", "#E18727FF", "#20854EFF",
  "#7876B1FF", "#6F99ADFF", "#FFDC91FF", "#EE4C97FF"
)

step6_output_dir <- "5.空转定位图"

spatial_plot_name <- "空转给色定位图"
spatial_plot_width <- 6.5
spatial_plot_height <- 5.5

combined_plot_name <- "组合图"
combined_plot_width <- 13
combined_plot_height <- 5.5

# 6.2 创建输出文件夹

if (!dir.exists(step6_output_dir)) {
  dir.create(step6_output_dir, recursive = TRUE)
}

# 6.3 检查UMAP结果是否存在

if (!"umap" %in% Reductions(Spatial_Data)) {
  stop("对象中不存在 umap，请先完成第5步。")
}

# 6.4 根据cluster数量整理空间图颜色

num_clusters <- length(levels(Idents(Spatial_Data)))

if (length(palette_spatial) < num_clusters) {
  stop("颜色数量不足以覆盖所有 cluster。")
}

cluster_colors <- setNames(
  palette_spatial[1:num_clusters],
  levels(Idents(Spatial_Data))
)

# 6.5 绘制空间定位图

spatial_plot <- suppressWarnings(
  SpatialPlot(
    Spatial_Data,
    label = TRUE,
    label.size = spatial_label_size,
    cols = cluster_colors
  )
)

print(spatial_plot)

ggsave(
  filename = file.path(
    step6_output_dir,
    paste0(spatial_plot_name, ".pdf")
  ),
  plot = spatial_plot,
  width = spatial_plot_width,
  height = spatial_plot_height,
  device = "pdf"
)

# 6.6 绘制用于组合图的UMAP

combined_umap_plot <- suppressWarnings(
  DimPlot(
    Spatial_Data,
    reduction = "umap",
    label = TRUE,
    label.size = 3.5,
    pt.size = pt_size_umap,
    cols = palette_spatial,
    raster = FALSE
  ) +
    theme_classic() +
    theme(
      panel.border = element_rect(
        fill = NA,
        color = "black",
        size = 0.5
      ),
      legend.position = "right"
    )
)

# 6.7 拼接空间定位图与UMAP图

combined_plot <- spatial_plot + combined_umap_plot

print(combined_plot)

ggsave(
  filename = file.path(
    step6_output_dir,
    paste0(combined_plot_name, ".pdf")
  ),
  plot = combined_plot,
  width = combined_plot_width,
  height = combined_plot_height,
  device = "pdf"
)

# 6.8 记录本步骤参数

param_step6_text <- paste0(
  "Step 6：空间定位与组合图\n",
  "运行流程说明：\n",
  "1. 使用 Seurat::SpatialPlot() 绘制 cluster 在组织切片中的空间定位图。\n",
  "2. 使用 DimPlot() 绘制 UMAP，并通过 patchwork 拼接为空间图与 UMAP 的组合图。\n",
  "3. 所有分析结果继续保存在 Spatial_Data 中。\n\n",
  "本次运行参数：\n",
  "- spatial label.size: ", spatial_label_size, "\n"
)

# 7. 保存全部步骤的参数记录

# 7.1 参数设置

parameter_file_name <- "Spatial_analysis_parameters"

# 7.2 合并并保存参数记录

parameter_summary <- paste(
  param_step1_text,
  param_step2_text,
  param_step3_text,
  param_step4_text,
  param_step5_text,
  param_step6_text,
  sep = "\n\n------------------------------\n\n"
)

writeLines(
  parameter_summary,
  con = paste0(parameter_file_name, ".txt")
)