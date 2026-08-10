suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(monocle)
  library(Biobase)
})

# 1. 读取Spatial_Data对象

# 1.1 参数设置

out_dir <- "Monocle2结果"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# 1.2 检查并读取Spatial_Data

if (!exists("Spatial_Data", envir = .GlobalEnv)) {
  stop("全局环境中没有 Spatial_Data 对象。")
}

seurat_obj <- get("Spatial_Data", envir = .GlobalEnv)

if (!inherits(seurat_obj, "Seurat")) {
  stop("Spatial_Data 不是有效的Seurat对象。")
}

if (!"seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
  stop("Spatial_Data@meta.data 中不存在 seurat_clusters 列，无法绘制cluster轨迹图。")
}

# 2. 构建Monocle2对象并筛选轨迹构建基因

# 2.1 参数设置

monocle2cores <- 1
mono_mean_expression <- 0.1

# 2.2 提取counts表达矩阵

expression_data <- as(
  GetAssayData(
    object = seurat_obj,
    layer = "counts"
  ),
  "sparseMatrix"
)

# 2.3 构建细胞注释信息

pd <- new(
  "AnnotatedDataFrame",
  data = seurat_obj@meta.data
)

# 2.4 构建基因注释信息

feature_data <- data.frame(
  gene_short_name = rownames(expression_data),
  row.names = rownames(expression_data)
)

fd <- new(
  "AnnotatedDataFrame",
  data = feature_data
)

# 2.5 构建Monocle2 CellDataSet对象

mycds <- newCellDataSet(
  expression_data,
  phenoData = pd,
  featureData = fd,
  expressionFamily = negbinomial.size()
)

# 2.6 估计size factor和基因离散度

mycds <- estimateSizeFactors(mycds)

mycds <- estimateDispersions(
  mycds,
  cores = monocle2cores,
  relative_expr = TRUE
)

# 2.7 根据表达量和离散度筛选轨迹构建基因

disp_table <- dispersionTable(mycds)

disp.genes <- subset(
  disp_table,
  mean_expression >= mono_mean_expression &
    dispersion_empirical >= dispersion_fit
)$gene_id

if (length(disp.genes) == 0) {
  stop(
    paste0(
      "没有筛选到轨迹构建基因，请适当降低 mono_mean_expression。目前设置为：",
      mono_mean_expression
    )
  )
}

mycds <- setOrderingFilter(
  mycds,
  disp.genes
)

# 2.8 绘制轨迹构建基因图

ordering_genes_plot <- plot_ordering_genes(mycds)

print(ordering_genes_plot)

# 3. 保存轨迹构建基因图

# 3.1 参数设置

w_ordering_genes <- 8
h_ordering_genes <- 8
name_ordering_genes <- "0.轨迹构建基因"

# 3.2 保存PDF

ordering_genes_file <- file.path(
  out_dir,
  paste0(name_ordering_genes, ".pdf")
)

pdf(
  file = ordering_genes_file,
  width = w_ordering_genes,
  height = h_ordering_genes
)

print(ordering_genes_plot)

dev.off()

# 4. 使用DDRTree进行降维

mycds <- reduceDimension(
  mycds,
  max_components = 2,
  method = "DDRTree"
)

# 5. 第一遍拟时序排序及轨迹图

# 5.1 参数设置

w_first_pseudotime <- 8
h_first_pseudotime <- 8
name_first_pseudotime <- "1.第一遍轨迹图_Pseudotime"

w_first_cluster <- 8
h_first_cluster <- 8
name_first_cluster <- "2.第一遍轨迹图_seurat_clusters"

w_first_state <- 8
h_first_state <- 8
name_first_state <- "3.第一遍轨迹图_State"

# 5.2 第一遍排序

mycds_first <- orderCells(mycds)

# 5.3 绘制第一遍Pseudotime轨迹图

first_pseudotime_plot <- plot_cell_trajectory(
  mycds_first,
  color_by = "Pseudotime"
)

print(first_pseudotime_plot)

# 5.4 保存第一遍Pseudotime轨迹图

first_pseudotime_file <- file.path(
  out_dir,
  paste0(name_first_pseudotime, ".pdf")
)

pdf(
  file = first_pseudotime_file,
  width = w_first_pseudotime,
  height = h_first_pseudotime
)

print(first_pseudotime_plot)

dev.off()

# 5.5 绘制第一遍seurat_clusters轨迹图

first_cluster_plot <- plot_cell_trajectory(
  mycds_first,
  color_by = "seurat_clusters"
)

print(first_cluster_plot)

# 5.6 保存第一遍seurat_clusters轨迹图

first_cluster_file <- file.path(
  out_dir,
  paste0(name_first_cluster, ".pdf")
)

pdf(
  file = first_cluster_file,
  width = w_first_cluster,
  height = h_first_cluster
)

print(first_cluster_plot)

dev.off()

# 5.7 绘制第一遍State轨迹图

first_state_plot <- plot_cell_trajectory(
  mycds_first,
  color_by = "State"
)

print(first_state_plot)

# 5.8 保存第一遍State轨迹图

first_state_file <- file.path(
  out_dir,
  paste0(name_first_state, ".pdf")
)

pdf(
  file = first_state_file,
  width = w_first_state,
  height = h_first_state
)

print(first_state_plot)

dev.off()

# 6. 指定起点进行第二遍拟时序排序

# 6.1 参数设置

root_state <- 1

w_second_pseudotime <- 8
h_second_pseudotime <- 8
name_second_pseudotime <- "4.第二遍轨迹图_Pseudotime"

w_second_cluster <- 8
h_second_cluster <- 8
name_second_cluster <- "5.第二遍轨迹图_seurat_clusters"

w_second_state <- 8
h_second_state <- 8
name_second_state <- "6.第二遍轨迹图_State"

# 6.2 查看并检查第一遍产生的State

available_states <- sort(
  unique(
    pData(mycds_first)$State
  )
)

if (!root_state %in% available_states) {
  stop(
    paste0(
      "root_state = ",
      root_state,
      " 不在当前State中。可选值为：",
      paste(available_states, collapse = ", ")
    )
  )
}

# 6.3 按指定起点重新排序

mycds_second <- orderCells(
  mycds_first,
  root_state = root_state
)

# 6.4 绘制第二遍Pseudotime轨迹图

second_pseudotime_plot <- plot_cell_trajectory(
  mycds_second,
  color_by = "Pseudotime"
)

print(second_pseudotime_plot)

# 6.5 保存第二遍Pseudotime轨迹图

second_pseudotime_file <- file.path(
  out_dir,
  paste0(name_second_pseudotime, ".pdf")
)

pdf(
  file = second_pseudotime_file,
  width = w_second_pseudotime,
  height = h_second_pseudotime
)

print(second_pseudotime_plot)

dev.off()

# 6.6 绘制第二遍seurat_clusters轨迹图

second_cluster_plot <- plot_cell_trajectory(
  mycds_second,
  color_by = "seurat_clusters"
)

print(second_cluster_plot)

# 6.7 保存第二遍seurat_clusters轨迹图

second_cluster_file <- file.path(
  out_dir,
  paste0(name_second_cluster, ".pdf")
)

pdf(
  file = second_cluster_file,
  width = w_second_cluster,
  height = h_second_cluster
)

print(second_cluster_plot)

dev.off()

# 6.8 绘制第二遍State轨迹图

second_state_plot <- plot_cell_trajectory(
  mycds_second,
  color_by = "State"
)

print(second_state_plot)

# 6.9 保存第二遍State轨迹图

second_state_file <- file.path(
  out_dir,
  paste0(name_second_state, ".pdf")
)

pdf(
  file = second_state_file,
  width = w_second_state,
  height = h_second_state
)

print(second_state_plot)

dev.off()

# 6.10 将第二遍结果作为最终mycds对象

mycds <- mycds_second

# 7. 绘制空间拟时序分布图

# 7.1 参数设置

w_spatial_pseudotime <- 8
h_spatial_pseudotime <- 8
name_spatial_pseudotime <- "2.空转时序图"

# 7.2 提取第二遍排序后的Pseudotime

spot_Pseudotime <- pData(mycds_second)[
  ,
  "Pseudotime",
  drop = FALSE
]

# 7.3 匹配Monocle2细胞与Spatial_Data中的spot

common_cells <- intersect(
  rownames(spot_Pseudotime),
  colnames(seurat_obj)
)

if (length(common_cells) == 0) {
  stop("Monocle2对象和Spatial_Data之间没有匹配的spot/cell名称。")
}

metadata_matched <- spot_Pseudotime[
  common_cells,
  ,
  drop = FALSE
]

# 7.4 将Pseudotime添加到空间对象

Spatial_Data_with_pseudotime <- AddMetaData(
  seurat_obj,
  metadata = metadata_matched
)

# 7.5 绘制空间拟时序图

spatial_pseudotime_plot <- SpatialFeaturePlot(
  Spatial_Data_with_pseudotime,
  features = "Pseudotime"
)

print(spatial_pseudotime_plot)

# 7.6 保存空间拟时序图

spatial_pseudotime_file <- file.path(
  out_dir,
  paste0(name_spatial_pseudotime, ".pdf")
)

pdf(
  file = spatial_pseudotime_file,
  width = w_spatial_pseudotime,
  height = h_spatial_pseudotime
)

print(spatial_pseudotime_plot)

dev.off()

# 8. 保存参数记录

# 8.1 参数设置

name_params <- "pseudotime_parameters"

# 8.2 生成参数记录

param_text <- paste0(
  "本次拟时序分析参数总结：\n",
  "- monocle2cores：", monocle2cores, "\n",
  "- mono_mean_expression：", mono_mean_expression, "\n",
  "- 轨迹构建基因数量：", length(disp.genes), "\n",
  "- 降维方法：DDRTree\n",
  "- max_components：2\n",
  "- 第一遍：orderCells(mycds)\n",
  "- 第一遍State：", paste(available_states, collapse = ", "), "\n",
  "- 第二遍：orderCells(mycds, root_state = ", root_state, ")\n",
  "- 最终结果对象：mycds\n",
  "- 第一遍结果对象：mycds_first\n",
  "- 第二遍结果对象：mycds_second\n",
  "- 添加Pseudotime的空间对象：Spatial_Data_with_pseudotime\n",
  "- 输出目录：", out_dir, "\n"
)

# 8.3 保存参数记录

param_file <- file.path(
  out_dir,
  paste0(name_params, ".txt")
)

writeLines(
  param_text,
  con = param_file
)