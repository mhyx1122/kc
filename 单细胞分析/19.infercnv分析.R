# 1. 加载必要 R 包

suppressPackageStartupMessages({
  library(Seurat)
  library(infercnv)
  library(Matrix)
})


# 2. 检查全局环境中的 seurat 对象

if (!exists("seurat", envir = .GlobalEnv)) {
  stop("全局环境中没有 seurat 对象，请先加载 seurat 对象。")
}

srt <- get("seurat", envir = .GlobalEnv)

if (!inherits(srt, "Seurat")) {
  stop("全局环境中的 seurat 不是 Seurat 对象。")
}

if (!"cellType" %in% colnames(srt@meta.data)) {
  stop("seurat@meta.data 中没有 cellType 列，请先添加细胞类型注释。")
}


# 3. 设置 InferCNV 分析参数

sample_size <- 500

unique_cell_types <- unique(srt@meta.data$cellType)
unique_cell_types <- unique_cell_types[!is.na(unique_cell_types)]

ref_group_names <- unique_cell_types[1]

cutoff_value <- 0.1

output_directory <- "cnv_analysis/"

cluster_by_groups <- TRUE

hclust_method <- "ward.D2"

plot_steps <- FALSE

write_expr_matrix <- TRUE

color_palette410 <- c("#8DD3C7", "white", "#BC80BD")

gene_order_file <- "gene_order_file.txt"


# 4. 创建输出目录

if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
}


# 5. 检查 gene_order_file 文件

if (!file.exists(gene_order_file)) {
  stop(
    "没有找到 gene_order_file 文件。请准备 inferCNV 需要的基因位置信息文件，并修改 gene_order_file 路径。\n",
    "gene_order_file 通常包含 4 列：gene、chr、start、end。"
  )
}


# 6. 抽样细胞

set.seed(123)

all_cells <- colnames(srt)

if (sample_size >= length(all_cells)) {
  sampled_cells <- all_cells
} else {
  sampled_cells <- sample(all_cells, sample_size)
}

srt_sub <- subset(srt, cells = sampled_cells)


# 7. 提取表达矩阵

default_assay <- DefaultAssay(srt_sub)

expr_matrix <- tryCatch(
  {
    GetAssayData(
      object = srt_sub,
      assay = default_assay,
      layer = "counts"
    )
  },
  error = function(e) {
    GetAssayData(
      object = srt_sub,
      assay = default_assay,
      slot = "counts"
    )
  }
)

if (nrow(expr_matrix) == 0 || ncol(expr_matrix) == 0) {
  stop("提取到的表达矩阵为空，请检查 Seurat 对象中的 counts 数据。")
}


# 8. 生成细胞注释文件

annotation_df <- data.frame(
  cell = colnames(expr_matrix),
  group = srt_sub@meta.data[colnames(expr_matrix), "cellType"],
  stringsAsFactors = FALSE
)

annotation_df <- annotation_df[!is.na(annotation_df$group), ]

expr_matrix <- expr_matrix[, annotation_df$cell, drop = FALSE]

annotation_file <- file.path(output_directory, "cell_annotations.txt")

write.table(
  annotation_df,
  file = annotation_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)


# 9. 写出表达矩阵文件

expr_matrix_file <- file.path(output_directory, "raw_counts_matrix.txt")

write.table(
  as.matrix(expr_matrix),
  file = expr_matrix_file,
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)


# 10. 检查参考组是否存在

ref_group_names <- intersect(ref_group_names, unique(annotation_df$group))

if (length(ref_group_names) == 0) {
  stop("设置的 ref_group_names 在 annotation_df$group 中不存在，请重新设置参考组。")
}


# 11. 创建 inferCNV 对象

infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = expr_matrix_file,
  annotations_file = annotation_file,
  delim = "\t",
  gene_order_file = gene_order_file,
  ref_group_names = ref_group_names
)


# 12. 运行 inferCNV 分析

infercnv_obj <- infercnv::run(
  infercnv_obj = infercnv_obj,
  cutoff = cutoff_value,
  out_dir = output_directory,
  cluster_by_groups = cluster_by_groups,
  hclust_method = hclust_method,
  plot_steps = plot_steps,
  write_expr_matrix = write_expr_matrix,
  denoise = TRUE,
  HMM = FALSE
)


# 13. 保存 inferCNV 对象

saveRDS(
  infercnv_obj,
  file = file.path(output_directory, "infercnv_obj.rds")
)