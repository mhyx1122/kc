suppressPackageStartupMessages({
  library(Seurat)
})

# 1. 检查全局环境中是否存在 seurat 对象
if (!exists("seurat", envir = .GlobalEnv)) {
  stop("全局环境中没有 seurat 对象")
}

# 2. 检查 group 列是否存在
if (!"group" %in% colnames(seurat@meta.data)) {
  stop("seurat@meta.data 中没有 group 列")
}


# 3. 查看当前 Seurat 对象中的分组

group_types <- unique(seurat@meta.data$group)


# 4. 设置要提取的分组

# 4.1 默认示例：提取第一个分组
selected_groups <- group_types[1]

# 4.2 如果要提取多个分组，可以改成类似下面这样
# selected_groups <- c("Control", "Treat")


# 5. 检查选择的分组是否存在

if (!all(selected_groups %in% unique(seurat$group))) {
  stop("某些分组在 Seurat 对象中不存在，请检查 selected_groups 参数。")
}


# 6. 提取指定分组

seurat <- subset(
  seurat,
  subset = group %in% selected_groups
)