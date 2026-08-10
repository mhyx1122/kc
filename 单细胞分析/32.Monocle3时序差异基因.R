# 1. 加载必要的 R 包

suppressPackageStartupMessages({
  library(Seurat)
  library(monocle3)
  library(ggplot2)
  library(dplyr)
})


# 2. 检查分析所需对象

out_dir1 <- "monocle3拟时序分析"
out_dir2 <- "monocle3拟时序分析"

if (!exists("cds", envir = .GlobalEnv)) {
  stop("全局环境中没有 cds 对象，请先运行前面的 Monocle3 分析代码。")
}

if (!exists("seurat", envir = .GlobalEnv)) {
  stop("全局环境中没有 seurat 对象。")
}

cds_obj <- get("cds", envir = .GlobalEnv)
seurat_obj <- get("seurat", envir = .GlobalEnv)

if (!inherits(seurat_obj, "Seurat")) {
  stop("全局环境中的 seurat 不是 Seurat 对象。")
}

if (is.null(cds_obj@principal_graph_aux$UMAP$pseudotime)) {
  stop(
    paste0(
      "cds 中不存在 UMAP pseudotime，",
      "请先完成 order_cells(cds, reduction_method = \"UMAP\")。"
    )
  )
}

if (!dir.exists(out_dir1)) {
  dir.create(out_dir1, recursive = TRUE)
}

if (!dir.exists(out_dir2)) {
  dir.create(out_dir2, recursive = TRUE)
}


# 3. 绘制并保存 Monocle3 拟时序图

pseudotime_width <- 8
pseudotime_height <- 6
pseudotime_file <- "Trajectory_Pseudotime.pdf"

p_pseudotime <- plot_cells(
  cds_obj,
  color_cells_by = "pseudotime",
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE
)

ggsave(
  filename = file.path(out_dir1, pseudotime_file),
  plot = p_pseudotime,
  width = pseudotime_width,
  height = pseudotime_height,
  device = "pdf"
)

message(
  "拟时序图已保存至：",
  file.path(out_dir1, pseudotime_file)
)


# 4. 使用 graph_test 筛选轨迹相关基因

cores <- 8

min_expr_tab2 <- 0.5
porp_tab2 <- 30
top_gene_n_tab2 <- 10
color_cells_by_tab2 <- "cellType"

gene_jitter_width <- 8
gene_jitter_height <- 6
gene_jitter_file <- "Genes_Jitterplot.pdf"

trajectory_csv_file <- "Trajectory_genes.csv"

if (!color_cells_by_tab2 %in% colnames(colData(cds_obj))) {
  stop(
    paste0(
      "cds 的 colData 中不存在列：",
      color_cells_by_tab2
    )
  )
}

corescg112 <- cores

message("正在运行 graph_test() 筛选轨迹相关基因")

trajectory_genes <- graph_test(
  cds_obj,
  neighbor_graph = "principal_graph",
  cores = corescg112
)

trajectory_genes <- trajectory_genes[, c(5, 2, 3, 4, 1, 6)] %>%
  dplyr::filter(q_value < 1e-3)


# 5. 根据有效拟时序细胞中的表达比例过滤轨迹相关基因

pseudotime_values <- cds_obj@principal_graph_aux$UMAP$pseudotime

valid_pseudotime_cells <- rownames(cds_obj@colData)[
  !is.na(pseudotime_values) &
    is.finite(pseudotime_values)
]

expr_matrix <- GetAssayData(
  seurat_obj,
  assay = "RNA",
  layer = "counts"
)[
  ,
  valid_pseudotime_cells,
  drop = FALSE
]

genes_to_filter <- trajectory_genes$gene_short_name

genes_present <- genes_to_filter[
  genes_to_filter %in% rownames(expr_matrix)
]

expr_subset <- expr_matrix[
  genes_present,
  ,
  drop = FALSE
]

num_cells_expressed <- Matrix::rowSums(
  expr_subset > 0
)

percent_expressed <- (
  num_cells_expressed / ncol(expr_subset)
) * 100

genes_filtered <- names(percent_expressed)[
  percent_expressed >= porp_tab2
]

trajectory_genes <- trajectory_genes[
  trajectory_genes$gene_short_name %in% genes_filtered,
]


# 6. 提取 Top 轨迹相关基因

top_trajectory_genes <- trajectory_genes %>%
  dplyr::top_n(
    n = top_gene_n_tab2,
    wt = morans_I
  ) %>%
  dplyr::pull(gene_short_name) %>%
  as.character()


# 7. 绘制并保存 Top 轨迹相关基因拟时序表达图

p_gene_jitter <- plot_genes_in_pseudotime(
  cds_obj[top_trajectory_genes, ],
  color_cells_by = color_cells_by_tab2,
  min_expr = min_expr_tab2,
  ncol = 2
)

ggsave(
  filename = file.path(out_dir1, gene_jitter_file),
  plot = p_gene_jitter,
  width = gene_jitter_width,
  height = gene_jitter_height,
  device = "pdf"
)

message(
  "Top 轨迹相关基因拟时序表达图已保存至：",
  file.path(out_dir1, gene_jitter_file)
)


# 8. 保存轨迹相关基因结果

write.csv(
  trajectory_genes,
  file = file.path(out_dir1, trajectory_csv_file),
  row.names = FALSE
)

message(
  "轨迹相关基因结果已保存至：",
  file.path(out_dir1, trajectory_csv_file)
)


# 9. 使用 Seurat FeaturePlot 展示指定基因

genes_plot_text_tab3 <- "AIF1,C1QA,C1QC"
umap_color_text_tab3 <- "grey,#FF0000"

pt_size_tab3 <- 0.5
color_seed <- 1234

feature_seurat_width <- 12
feature_seurat_height <- 12
feature_seurat_file <- "seurat展示指定基因分布.pdf"

genes_plot_tab3 <- trimws(
  unlist(
    strsplit(
      genes_plot_text_tab3,
      split = ","
    )
  )
)

genes_plot_tab3 <- genes_plot_tab3[
  genes_plot_tab3 != ""
]

umap_color_tab3 <- trimws(
  unlist(
    strsplit(
      umap_color_text_tab3,
      split = ","
    )
  )
)

umap_color_tab3 <- umap_color_tab3[
  umap_color_tab3 != ""
]

set.seed(color_seed)

p_feature_seurat <- FeaturePlot(
  seurat_obj,
  features = genes_plot_tab3,
  reduction = "umap",
  cols = umap_color_tab3,
  min.cutoff = 0,
  max.cutoff = NA,
  ncol = 2,
  pt.size = pt_size_tab3
)

ggsave(
  filename = file.path(out_dir2, feature_seurat_file),
  plot = p_feature_seurat,
  width = feature_seurat_width,
  height = feature_seurat_height,
  device = "pdf"
)

message(
  "Seurat 指定基因表达分布图已保存至：",
  file.path(out_dir2, feature_seurat_file)
)


# 10. 使用 Monocle3 cds 对象展示指定基因

genes_plot_text_tab4 <- "AIF1,C1QA,C1QC"

feature_cds_width <- 12
feature_cds_height <- 12
feature_cds_file <- "cds展示指定基因分布.pdf"

genes_plot_tab4 <- trimws(
  unlist(
    strsplit(
      genes_plot_text_tab4,
      split = ","
    )
  )
)

genes_plot_tab4 <- genes_plot_tab4[
  genes_plot_tab4 != ""
]

p_feature_cds <- plot_cells(
  cds_obj,
  genes = genes_plot_tab4,
  label_cell_groups = FALSE,
  show_trajectory_graph = FALSE
)

ggsave(
  filename = file.path(out_dir2, feature_cds_file),
  plot = p_feature_cds,
  width = feature_cds_width,
  height = feature_cds_height,
  device = "pdf"
)

message(
  "Monocle3 cds 指定基因表达分布图已保存至：",
  file.path(out_dir2, feature_cds_file)
)


# 11. 绘制指定基因拟时序表达图

genes_plot_text_tab5 <- "AIF1,C1QA,C1QC"

min_expr_tab5 <- 0.5
cell_size_tab5 <- 1.5

gene_pseudotime_width <- 8
gene_pseudotime_height <- 8
gene_pseudotime_file <- "拟时序差异基因分布.pdf"

if (!"cellType" %in% colnames(colData(cds_obj))) {
  stop("cds 的 colData 中不存在列：cellType。")
}

genes_plot_tab5 <- trimws(
  unlist(
    strsplit(
      genes_plot_text_tab5,
      split = ","
    )
  )
)

genes_plot_tab5 <- genes_plot_tab5[
  genes_plot_tab5 != ""
]

genes_cds <- cds_obj[
  rowData(cds_obj)$gene_short_name %in% genes_plot_tab5,
]

p_gene_pseudotime <- plot_genes_in_pseudotime(
  genes_cds,
  color_cells_by = "cellType",
  min_expr = min_expr_tab5,
  cell_size = cell_size_tab5
)

ggsave(
  filename = file.path(out_dir2, gene_pseudotime_file),
  plot = p_gene_pseudotime,
  width = gene_pseudotime_width,
  height = gene_pseudotime_height,
  device = "pdf"
)

message(
  "指定基因拟时序表达图已保存至：",
  file.path(out_dir2, gene_pseudotime_file)
)


# 12. 保存本次 Monocle3 后续分析参数

color_cells_by <- "cellType"
params_file <- "monocle3_followup_parameters.txt"

param_text <- paste0(
  "本次 Monocle3 后续分析参数总结：\n",
  "- 固定输出文件夹1：", out_dir1, "\n",
  "- 固定输出文件夹2：", out_dir2, "\n",
  "- cores：", cores, "\n",
  "- color_cells_by：", color_cells_by, "\n",
  "- color_seed：", color_seed, "\n",
  "- graph_test min_expr：", min_expr_tab2, "\n",
  "- graph_test porp：", porp_tab2, "\n",
  "- graph_test top_gene_n：", top_gene_n_tab2, "\n",
  "- graph_test color_cells_by：", color_cells_by_tab2, "\n",
  "- Seurat指定基因：",
  paste(genes_plot_tab3, collapse = ", "),
  "\n",
  "- FeaturePlot颜色：",
  paste(umap_color_tab3, collapse = ", "),
  "\n",
  "- FeaturePlot点大小：", pt_size_tab3, "\n",
  "- cds指定基因：",
  paste(genes_plot_tab4, collapse = ", "),
  "\n",
  "- 拟时序表达指定基因：",
  paste(genes_plot_tab5, collapse = ", "),
  "\n",
  "- 指定基因拟时序 min_expr：", min_expr_tab5, "\n",
  "- 指定基因拟时序 cell_size：", cell_size_tab5, "\n",
  "- Trajectory_genes.csv：",
  file.path(out_dir1, trajectory_csv_file),
  "\n"
)

writeLines(
  text = param_text,
  con = file.path(out_dir1, params_file)
)

message(
  "Monocle3 后续分析参数已保存至：",
  file.path(out_dir1, params_file)
)