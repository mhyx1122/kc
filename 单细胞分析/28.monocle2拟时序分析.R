# 1. 加载必要的 R 包

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(Biobase)
  library(monocle)
  library(ggplot2)
})


# 2. 检查 Seurat 对象并创建输出文件夹

out_dir <- "Monocle2结果"

if (!exists("seurat", envir = .GlobalEnv)) {
  stop("全局环境中没有找到 seurat 对象，请先加载 Seurat 对象。")
}

seurat_obj <- get("seurat", envir = .GlobalEnv)

if (!inherits(seurat_obj, "Seurat")) {
  stop("全局环境中的 seurat 不是 Seurat 对象。")
}

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}


# 3. 从 Seurat 对象构建 Monocle2 对象

monocle2cores <- 1

mono_mean_expression <- 0.1

message("正在提取 Seurat 对象中的原始表达矩阵")

expression_data <- as(
  GetAssayData(
    object = seurat_obj,
    layer = "counts"
  ),
  "sparseMatrix"
)

cell_metadata <- seurat_obj@meta.data

if (!identical(colnames(expression_data), rownames(cell_metadata))) {
  cell_metadata <- cell_metadata[colnames(expression_data), , drop = FALSE]
}

pd <- new(
  "AnnotatedDataFrame",
  data = cell_metadata
)

gene_metadata <- data.frame(
  gene_short_name = rownames(expression_data),
  row.names = rownames(expression_data),
  stringsAsFactors = FALSE
)

fd <- new(
  "AnnotatedDataFrame",
  data = gene_metadata
)

message("正在构建 Monocle2 CellDataSet 对象")

mycds <- newCellDataSet(
  expression_data,
  phenoData = pd,
  featureData = fd,
  expressionFamily = negbinomial.size()
)

message("正在估计细胞大小因子")

mycds <- estimateSizeFactors(mycds)

message("正在估计基因离散度")

mycds <- estimateDispersions(
  mycds,
  cores = monocle2cores,
  relative_expr = TRUE
)


# 4. 筛选轨迹构建基因

disp_table <- dispersionTable(mycds)

disp.genes <- subset(
  disp_table,
  mean_expression >= mono_mean_expression &
    dispersion_empirical >= 1 * dispersion_fit
)$gene_id

if (length(disp.genes) == 0) {
  stop(
    "没有筛选到符合条件的轨迹构建基因，请适当降低 mono_mean_expression。"
  )
}

message(
  "筛选得到的轨迹构建基因数量：",
  length(disp.genes)
)

mycds <- setOrderingFilter(
  mycds,
  disp.genes
)


# 5. 绘制并保存轨迹构建基因图

ordering_genes_width <- 8

ordering_genes_height <- 8

ordering_genes_file <- "0.轨迹构建基因.pdf"

ordering_plot <- plot_ordering_genes(mycds)

pdf(
  file = file.path(out_dir, ordering_genes_file),
  width = ordering_genes_width,
  height = ordering_genes_height
)

print(ordering_plot)

dev.off()

message(
  "轨迹构建基因图已保存至：",
  file.path(out_dir, ordering_genes_file)
)


# 6. DDRTree 降维并进行第一遍细胞排序

message("正在进行 DDRTree 降维")

mycds <- reduceDimension(
  mycds,
  max_components = 2,
  method = "DDRTree"
)

message("正在进行第一遍细胞排序")

mycds_first <- orderCells(mycds)

mycds <- mycds_first


# 7. 绘制并保存第一遍 Pseudotime 轨迹图

first_pseudotime_width <- 8

first_pseudotime_height <- 8

first_pseudotime_file <- "1.第一遍轨迹图_Pseudotime.pdf"

p_first_pseudotime <- tryCatch(
  {
    plot_cell_trajectory(
      mycds_first,
      color_by = "Pseudotime"
    )
  },
  error = function(e) {
    ggplot() +
      theme_void() +
      ggtitle(
        paste0(
          "无法绘制 Pseudotime：",
          e$message
        )
      )
  }
)

pdf(
  file = file.path(out_dir, first_pseudotime_file),
  width = first_pseudotime_width,
  height = first_pseudotime_height
)

print(p_first_pseudotime)

dev.off()

message(
  "第一遍 Pseudotime 轨迹图已保存至：",
  file.path(out_dir, first_pseudotime_file)
)


# 8. 绘制并保存第一遍 seurat_clusters 轨迹图

first_cluster_width <- 8

first_cluster_height <- 8

first_cluster_file <- "2.第一遍轨迹图_seurat_clusters.pdf"

p_first_cluster <- tryCatch(
  {
    plot_cell_trajectory(
      mycds_first,
      color_by = "seurat_clusters"
    )
  },
  error = function(e) {
    ggplot() +
      theme_void() +
      ggtitle(
        paste0(
          "无法绘制 seurat_clusters：",
          e$message
        )
      )
  }
)

pdf(
  file = file.path(out_dir, first_cluster_file),
  width = first_cluster_width,
  height = first_cluster_height
)

print(p_first_cluster)

dev.off()

message(
  "第一遍 seurat_clusters 轨迹图已保存至：",
  file.path(out_dir, first_cluster_file)
)


# 9. 绘制并保存第一遍 State 轨迹图

first_state_width <- 8

first_state_height <- 8

first_state_file <- "3.第一遍轨迹图_State.pdf"

p_first_state <- tryCatch(
  {
    plot_cell_trajectory(
      mycds_first,
      color_by = "State"
    )
  },
  error = function(e) {
    ggplot() +
      theme_void() +
      ggtitle(
        paste0(
          "无法绘制 State：",
          e$message
        )
      )
  }
)

pdf(
  file = file.path(out_dir, first_state_file),
  width = first_state_width,
  height = first_state_height
)

print(p_first_state)

dev.off()

message(
  "第一遍 State 轨迹图已保存至：",
  file.path(out_dir, first_state_file)
)


# 10. 查看第一遍轨迹中的 State 编号

available_states <- sort(
  unique(
    pData(mycds_first)$State
  )
)

cat(
  "第一遍轨迹中可选择的 State 编号：",
  paste(available_states, collapse = ", "),
  "\n"
)


# 11. 指定轨迹起点并进行第二遍细胞排序

# 请根据第一遍 State 轨迹图设置起始 State
root_state <- 1

if (!(root_state %in% available_states)) {
  stop(
    paste0(
      "root_state = ",
      root_state,
      " 不在当前 State 中。可选值为：",
      paste(available_states, collapse = ", ")
    )
  )
}

message(
  "正在以 State ",
  root_state,
  " 作为起点重新计算拟时序"
)

mycds_second <- orderCells(
  mycds_first,
  root_state = root_state
)

mycds <- mycds_second


# 12. 绘制并保存第二遍 Pseudotime 轨迹图

second_pseudotime_width <- 8

second_pseudotime_height <- 8

second_pseudotime_file <- "4.第二遍轨迹图_Pseudotime.pdf"

p_second_pseudotime <- tryCatch(
  {
    plot_cell_trajectory(
      mycds_second,
      color_by = "Pseudotime"
    )
  },
  error = function(e) {
    ggplot() +
      theme_void() +
      ggtitle(
        paste0(
          "无法绘制 Pseudotime：",
          e$message
        )
      )
  }
)

pdf(
  file = file.path(out_dir, second_pseudotime_file),
  width = second_pseudotime_width,
  height = second_pseudotime_height
)

print(p_second_pseudotime)

dev.off()

message(
  "第二遍 Pseudotime 轨迹图已保存至：",
  file.path(out_dir, second_pseudotime_file)
)


# 13. 绘制并保存第二遍 seurat_clusters 轨迹图

second_cluster_width <- 8

second_cluster_height <- 8

second_cluster_file <- "5.第二遍轨迹图_seurat_clusters.pdf"

p_second_cluster <- tryCatch(
  {
    plot_cell_trajectory(
      mycds_second,
      color_by = "seurat_clusters"
    )
  },
  error = function(e) {
    ggplot() +
      theme_void() +
      ggtitle(
        paste0(
          "无法绘制 seurat_clusters：",
          e$message
        )
      )
  }
)

pdf(
  file = file.path(out_dir, second_cluster_file),
  width = second_cluster_width,
  height = second_cluster_height
)

print(p_second_cluster)

dev.off()

message(
  "第二遍 seurat_clusters 轨迹图已保存至：",
  file.path(out_dir, second_cluster_file)
)


# 14. 绘制并保存第二遍 State 轨迹图

second_state_width <- 8

second_state_height <- 8

second_state_file <- "6.第二遍轨迹图_State.pdf"

p_second_state <- tryCatch(
  {
    plot_cell_trajectory(
      mycds_second,
      color_by = "State"
    )
  },
  error = function(e) {
    ggplot() +
      theme_void() +
      ggtitle(
        paste0(
          "无法绘制 State：",
          e$message
        )
      )
  }
)

pdf(
  file = file.path(out_dir, second_state_file),
  width = second_state_width,
  height = second_state_height
)

print(p_second_state)

dev.off()

message(
  "第二遍 State 轨迹图已保存至：",
  file.path(out_dir, second_state_file)
)


# 15. 保存本次拟时序分析参数

params_file <- "pseudotime_parameters.txt"

param_text <- paste0(
  "本次拟时序分析参数总结：\n",
  "- monocle2cores：", monocle2cores, "\n",
  "- mono_mean_expression：", mono_mean_expression, "\n",
  "- 轨迹构建基因数量：", length(disp.genes), "\n",
  "- 降维方法：DDRTree\n",
  "- 第一遍：orderCells(mycds)\n",
  "- 第二遍：orderCells(mycds_first, root_state = ",
  root_state,
  ")\n",
  "- 第一遍结果对象：mycds_first\n",
  "- 第二遍结果对象：mycds_second\n",
  "- 最终结果对象：mycds\n",
  "- 输出文件夹：", out_dir, "\n"
)

writeLines(
  text = param_text,
  con = file.path(out_dir, params_file)
)

cat(param_text)

message(
  "参数文件已保存至：",
  file.path(out_dir, params_file)
)

