suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

# 1. APP参数设置

# 选择 MAD 值排名前 n 个基因
n_genesmad <- 3000

# 保存文件夹名称
folder_path <- "转换为bulk的结果"

# 使用的 assay 和 layer
assay_use <- "RNA"
layer_use <- "data"

# 输出文件名
output_all_file <- "single_for_bulk_data.csv"
output_top_file <- "single_for_bulk_data(Top).csv"


# 2. 检查 seurat 对象

if (!exists("seurat", envir = .GlobalEnv)) {
  stop("全局环境中没有找到 seurat 对象，请先加载 seurat 对象。")
}

seurat_obj <- get("seurat", envir = .GlobalEnv)

if (!inherits(seurat_obj, "Seurat")) {
  stop("全局环境中的 seurat 不是 Seurat 对象。")
}

if (!"cellType" %in% colnames(seurat_obj@meta.data)) {
  stop("seurat@meta.data 中没有 cellType 列，请先完成细胞类型注释。")
}


# 3. 创建输出文件夹

if (!dir.exists(folder_path)) {
  dir.create(folder_path, recursive = TRUE)
}


# 4. 提取单细胞表达矩阵

expr_matrix <- as.matrix(
  GetAssayData(
    seurat_obj,
    assay = assay_use,
    layer = layer_use
  )
)
cell_types <- seurat_obj@meta.data[colnames(expr_matrix), "cellType"]
cell_types2 <- seurat_obj@meta.data$cellType
cell_ids <- colnames(expr_matrix)

if (length(cell_types) != length(cell_ids)) {
  stop("cellType 数量与表达矩阵细胞数量不一致，请检查 seurat 对象。")
}

if (any(is.na(cell_types))) {
  stop("cellType 中存在 NA，请先处理未注释细胞。")
}


# 5. 修改细胞列名：cellType_cellID

new_cell_ids <- paste(cell_types, cell_ids, sep = "_")
colnames(expr_matrix) <- new_cell_ids


# 6. 按 cellType 对细胞排序

sort_index <- order(cell_types)
expr_matrix_sorted <- expr_matrix[, sort_index, drop = FALSE]


# 7. 保存排序后的完整表达矩阵

write.csv(
  expr_matrix_sorted,
  file = file.path(folder_path, output_all_file),
  row.names = TRUE
)


# 8. 根据 MAD 值筛选前 n 个高变基因
data_mad <- apply(expr_matrix_sorted, 1, mad)

n_genesmad <- min(n_genesmad, nrow(expr_matrix_sorted))

top_gene_index <- order(data_mad, decreasing = TRUE)[1:n_genesmad]

datExpr0 <- expr_matrix_sorted[top_gene_index, , drop = FALSE]

# 9. 保存 MAD 前 n 个基因表达矩阵
write.csv(
  datExpr0,
  file = file.path(folder_path, output_top_file),
  row.names = TRUE
)


# 10. 完成提示

cat("处理完成！\n")
cat("完整矩阵已保存：", file.path(folder_path, output_all_file), "\n")
cat("Top MAD 基因矩阵已保存：", file.path(folder_path, output_top_file), "\n")
cat("筛选基因数量：", n_genesmad, "\n")
cat("细胞数量：", ncol(expr_matrix_sorted), "\n")
cat("基因数量：", nrow(expr_matrix_sorted), "\n")