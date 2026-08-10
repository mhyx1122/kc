library(Seurat)

# 1. 获取当前 Seurat 对象中的样本名称
sample_names <- unique(seurat@meta.data$orig.ident)

# 2. 查看样本名称
sample_names

# 3. 样本分组设置模块
# 这里的分组顺序必须和 sample_names 的顺序一一对应
# 例如 sample_names 显示为：
# "sample1" "sample2" "sample3" "sample4" "sample5" "sample6" "sample7" "sample8"
# 那么 sample_groups 也要写 8 个组别

sample_groups <- c(
  "normal",
  "normal",
  "normal",
  "normal",
  "lesion",
  "lesion",
  "lesion",
  "lesion"
)

# 4. 创建样本名和组别的对应关系
group_mapping <- setNames(sample_groups, sample_names)

# 5. 将分组信息添加到 Seurat 对象的 meta.data 中
seurat@meta.data$group <- group_mapping[seurat@meta.data$orig.ident]

# 6. 查看每个样本对应的分组
group_mapping

# 7. 查看分组后的细胞数量
table(seurat@meta.data$group)

# 8. 查看每个样本在不同分组中的细胞数量
table(seurat@meta.data$orig.ident, seurat@meta.data$group)