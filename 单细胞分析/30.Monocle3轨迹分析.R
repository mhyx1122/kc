# 1. 加载必要的 R 包

suppressPackageStartupMessages({
  library(Seurat)
  library(monocle3)
  library(ggplot2)
  library(patchwork)
})


# 2. 检查 Seurat 对象并创建输出文件夹

out_dir <- "monocle3拟时序分析"

if (!exists("seurat", envir = .GlobalEnv)) {
  stop("全局环境中没有找到 seurat 对象。")
}

seurat_obj <- get("seurat", envir = .GlobalEnv)

if (!inherits(seurat_obj, "Seurat")) {
  stop("全局环境中的 seurat 不是 Seurat 对象。")
}

if (!"RNA" %in% names(seurat_obj@assays)) {
  stop("seurat 对象中不存在 RNA assay。")
}

if (!"umap" %in% names(seurat_obj@reductions)) {
  stop("seurat 对象中不存在 umap 降维结果。")
}

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}


# 3. 使用 Seurat 对象构建 Monocle3 cell_data_set 对象

color_cells_by <- "cellType"

if (!color_cells_by %in% colnames(seurat_obj@meta.data)) {
  stop(
    paste0(
      "seurat@meta.data 中不存在列：",
      color_cells_by
    )
  )
}

message("正在提取 Seurat 对象中的 RNA counts 矩阵")

count_matrix <- GetAssayData(
  object = seurat_obj,
  assay = "RNA",
  layer = "counts"
)

cell_metadata <- seurat_obj@meta.data

if (!all(colnames(count_matrix) %in% rownames(cell_metadata))) {
  stop("表达矩阵中的部分细胞不存在于 seurat@meta.data 中。")
}

cell_metadata <- cell_metadata[
  colnames(count_matrix),
  ,
  drop = FALSE
]

gene_metadata <- data.frame(
  gene_short_name = rownames(count_matrix),
  row.names = rownames(count_matrix),
  stringsAsFactors = FALSE
)

message("正在构建 Monocle3 cell_data_set 对象")

cds <- new_cell_data_set(
  expression_data = count_matrix,
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata
)


# 4. 对 Monocle3 对象进行预处理

num_dim <- 20

message("正在运行 preprocess_cds()")

cds <- preprocess_cds(
  cds,
  num_dim = num_dim
)


# 5. 绘制并保存 PCA 方差解释图

pc_var_width <- 8
pc_var_height <- 6
pc_var_file <- "1.PCA方差解释图.pdf"

p_pc <- plot_pc_variance_explained(cds)

ggsave(
  filename = file.path(out_dir, pc_var_file),
  plot = p_pc,
  width = pc_var_width,
  height = pc_var_height,
  device = "pdf"
)

message(
  "PCA 方差解释图已保存至：",
  file.path(out_dir, pc_var_file)
)


# 6. 运行 Monocle3 UMAP 降维

cores <- 8
color_seed <- 1234

set.seed(color_seed)

message("正在运行 Monocle3 UMAP 降维")

cds <- reduce_dimension(
  cds,
  reduction_method = "UMAP",
  cores = cores,
  umap.fast_sgd = TRUE
)


# 7. 绘制并保存 Monocle3 原始 UMAP

umap_monocle_width <- 8
umap_monocle_height <- 6
umap_monocle_file <- "2.cds_umap.pdf"

p_umap_monocle <- plot_cells(
  cds,
  reduction_method = "UMAP",
  color_cells_by = color_cells_by
) +
  ggtitle("cds.umap")

ggsave(
  filename = file.path(out_dir, umap_monocle_file),
  plot = p_umap_monocle,
  width = umap_monocle_width,
  height = umap_monocle_height,
  device = "pdf"
)

message(
  "Monocle3 原始 UMAP 已保存至：",
  file.path(out_dir, umap_monocle_file)
)


# 8. 使用 Seurat UMAP 坐标替换 Monocle3 UMAP 坐标

message("正在使用 Seurat UMAP 坐标替换 Monocle3 UMAP 坐标")

cds_umap <- cds@int_colData$reducedDims$UMAP

seurat_umap <- Embeddings(
  object = seurat_obj,
  reduction = "umap"
)

if (!all(rownames(cds_umap) %in% rownames(seurat_umap))) {
  missing_cells <- setdiff(
    rownames(cds_umap),
    rownames(seurat_umap)
  )
  
  stop(
    paste0(
      "以下细胞不存在于 Seurat UMAP 中：",
      paste(head(missing_cells, 20), collapse = ", ")
    )
  )
}

seurat_umap <- seurat_umap[
  rownames(cds_umap),
  ,
  drop = FALSE
]

cds@int_colData$reducedDims$UMAP <- seurat_umap


# 9. 绘制并保存替换后的 Seurat UMAP

umap_integrated_width <- 8
umap_integrated_height <- 6
umap_integrated_file <- "3.int_umap.pdf"

p_umap_integrated <- plot_cells(
  cds,
  reduction_method = "UMAP",
  color_cells_by = color_cells_by
) +
  ggtitle("int.umap")

ggsave(
  filename = file.path(out_dir, umap_integrated_file),
  plot = p_umap_integrated,
  width = umap_integrated_width,
  height = umap_integrated_height,
  device = "pdf"
)

message(
  "替换后的 Seurat UMAP 已保存至：",
  file.path(out_dir, umap_integrated_file)
)


# 10. 绘制并保存两种 UMAP 对比图

umap_compare_width <- 10
umap_compare_height <- 5
umap_compare_file <- "UMAP_VS.pdf"

p_umap_compare <- p_umap_monocle | p_umap_integrated

ggsave(
  filename = file.path(out_dir, umap_compare_file),
  plot = p_umap_compare,
  width = umap_compare_width,
  height = umap_compare_height,
  device = "pdf"
)

message(
  "UMAP 对比图已保存至：",
  file.path(out_dir, umap_compare_file)
)


# 11. 对细胞进行聚类并学习轨迹主图

Setminimal_branch <- 10

message("正在进行 Louvain 聚类")

cds <- cluster_cells(
  cds,
  cluster_method = "louvain"
)

message("正在运行 learn_graph() 学习细胞轨迹")

cds <- learn_graph(
  cds,
  learn_graph_control = list(
    minimal_branch_len = Setminimal_branch
  )
)


# 12. 绘制并保存 learn_graph 轨迹图

trajectory_width <- 8
trajectory_height <- 6
trajectory_file <- "Trajectory.pdf"

p_trajectory <- plot_cells(
  cds,
  color_cells_by = color_cells_by,
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE
)

ggsave(
  filename = file.path(out_dir, trajectory_file),
  plot = p_trajectory,
  width = trajectory_width,
  height = trajectory_height,
  device = "pdf"
)

message(
  "learn_graph 轨迹图已保存至：",
  file.path(out_dir, trajectory_file)
)


# 13. 绘制并保存带叶节点和分支点标注的轨迹图

graph_label_size <- 1.5

trajectory_labeled_width <- 8
trajectory_labeled_height <- 6
trajectory_labeled_file <- "embryo.time.bin.pdf"

p_trajectory_labeled <- plot_cells(
  cds,
  color_cells_by = color_cells_by,
  label_cell_groups = FALSE,
  label_leaves = TRUE,
  label_branch_points = TRUE,
  graph_label_size = graph_label_size
)

ggsave(
  filename = file.path(out_dir, trajectory_labeled_file),
  plot = p_trajectory_labeled,
  width = trajectory_labeled_width,
  height = trajectory_labeled_height,
  device = "pdf"
)

message(
  "带叶节点和分支点标注的轨迹图已保存至：",
  file.path(out_dir, trajectory_labeled_file)
)


# 14. 保存本次 Monocle3 分析参数

params_file <- "monocle3_run_parameters.txt"

param_text <- paste0(
  "本次 Monocle3 轨迹学习参数总结：\n",
  "- cores：", cores, "\n",
  "- num_dim：", num_dim, "\n",
  "- minimal_branch_len：", Setminimal_branch, "\n",
  "- color_cells_by：", color_cells_by, "\n",
  "- graph_label_size：", graph_label_size, "\n",
  "- color_seed：", color_seed, "\n",
  "- 输出文件夹：", out_dir, "\n",
  "- 结果对象：cds\n",
  "\n",
  "本步骤实际运行的核心代码与流程包括：\n",
  "1. 使用 Seurat 对象中的 RNA counts 矩阵、细胞注释信息和基因注释信息构建 monocle3 的 cell_data_set 对象。\n",
  "2. 运行 preprocess_cds(cds, num_dim = ", num_dim, ") 对数据进行预处理和降维准备。\n",
  "3. 绘制 plot_pc_variance_explained(cds) 观察主成分解释方差情况。\n",
  "4. 运行 reduce_dimension(cds, reduction_method = \"UMAP\") 生成 monocle3 自身的 UMAP 结果。\n",
  "5. 提取 Seurat 对象中的 UMAP 坐标，并替换 monocle3 内部的 UMAP 坐标。\n",
  "6. 运行 cluster_cells(cds, cluster_method = \"louvain\") 对细胞进行聚类。\n",
  "7. 运行 learn_graph(cds, learn_graph_control = list(minimal_branch_len = ",
  Setminimal_branch,
  ")) 学习细胞轨迹主图结构。\n",
  "8. 绘制普通轨迹图以及带叶节点和分支点标注的轨迹图。\n",
  "\n",
  "本步骤的目的在于：\n",
  "1. 在 monocle3 框架下完成对象构建、降维、聚类和轨迹学习。\n",
  "2. 比较 monocle3 原始 UMAP 与 Seurat UMAP 的可视化差异。\n",
  "3. 获得可用于后续 order_cells() 排序的轨迹主图。\n"
)

writeLines(
  text = param_text,
  con = file.path(out_dir, params_file)
)

message(
  "Monocle3 分析参数已保存至：",
  file.path(out_dir, params_file)
)