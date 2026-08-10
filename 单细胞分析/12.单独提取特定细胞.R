suppressPackageStartupMessages({
  library(Seurat)
})

# 1. 检查全局环境中是否存在 seurat 对象
if (!exists("seurat", envir = .GlobalEnv)) {
  stop("全局环境中没有 seurat 对象")
}

# 2. 检查 cellType 列是否存在
if (!"cellType" %in% colnames(seurat@meta.data)) {
  stop("seurat@meta.data 中没有 cellType 列")
}


# 3. 查看当前 Seurat 对象中的细胞类型

cell_types <- unique(seurat@meta.data$cellType)


# 4. 设置要提取的细胞类型

# 4.1 默认示例：提取第一个细胞类型
selected_celltypes <- cell_types[1]

# 4.2 如果要提取多个细胞类型，可以改成类似下面这样
# selected_celltypes <- c("T cells", "B cells", "Macrophages")


# 5. 检查选择的细胞类型是否存在

if (!all(selected_celltypes %in% unique(seurat$cellType))) {
  stop("某些细胞类型在 Seurat 对象中不存在，请检查 selected_celltypes 参数。")
}


# 6. 提取指定细胞类型

seurat <- subset(
  seurat,
  subset = cellType %in% selected_celltypes
)