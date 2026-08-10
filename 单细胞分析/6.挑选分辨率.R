library(Seurat)
library(clustree)

# 1. 检查全局环境中是否存在 seurat 对象
if (!exists("seurat", envir = .GlobalEnv)) {
  stop("全局环境中没有 seurat 对象")
}

srt0 <- get("seurat", envir = .GlobalEnv)


# 2. 多分辨率聚类参数模块

# 2.1 降维方法参数
# 可选值需要存在于 Reductions(seurat) 中
# 常见可选值："pca"、"harmony"、"integrated.cca"、"integrated.rpca"、"integrated.mnn"、"integrated.Join"
reduction <- "pca"

# 2.2 使用的 PC 个数
pcs_end <- 20
dims_use <- 1:pcs_end

# 2.3 resolutions 参数
res_vec_text <- "0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6"

res_vec <- unlist(strsplit(res_vec_text, ","))
res_vec <- trimws(res_vec)
res_vec <- res_vec[res_vec != ""]
res_vec <- suppressWarnings(as.numeric(res_vec))

# 2.4 graph.name 参数
graph_name <- "RNA_snn"

# 2.5 FindClusters algorithm 参数
algorithm <- 1


# 3. 参数检查模块

# 3.1 检查 reduction 是否存在
if (!reduction %in% Reductions(srt0)) {
  stop(paste0("seurat 中不存在 reduction：", reduction))
}

# 3.2 检查 resolutions 是否解析成功
if (length(res_vec) == 0 || any(is.na(res_vec))) {
  stop("resolutions 向量解析失败：请用逗号分隔的数字，例如 0.2,0.4,0.6")
}


# 4. 输出目录模块

# 4.1 输出目录参数
out_dir <- "2.2选择分辨率"

# 4.2 创建输出目录
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}


# 5. 构建邻近图模块

# 5.1 基于指定 reduction 和 dims 运行 FindNeighbors
srt <- srt0

srt <- FindNeighbors(
  srt,
  dims = dims_use,
  reduction = reduction
)


# 6. 多分辨率聚类模块

# 6.1 循环不同 resolution 运行 FindClusters
for (i in seq_along(res_vec)) {
  
  cat("正在运行 resolution = ", res_vec[[i]], "  ", i, "/", length(res_vec), "\n", sep = "")
  
  srt <- FindClusters(
    srt,
    graph.name = graph_name,
    resolution = res_vec[[i]],
    algorithm = algorithm
  )
}


# 7. 提取各分辨率聚类结果模块

# 7.1 提取 meta.data
md <- srt@meta.data

# 7.2 查找当前 graph.name 对应的 resolution 聚类列
res_cols <- grep(
  paste0(graph_name, "_res."),
  colnames(md),
  value = TRUE,
  fixed = TRUE
)

# 7.3 兼容常见默认前缀 RNA_snn_res. 或 SCT_snn_res.
if (length(res_cols) == 0) {
  res_cols <- grep("RNA_snn_res.", colnames(md), value = TRUE, fixed = TRUE)
  res_cols <- c(
    res_cols,
    grep("SCT_snn_res.", colnames(md), value = TRUE, fixed = TRUE)
  )
}

# 7.4 检查是否找到聚类列
if (length(res_cols) == 0) {
  stop("未在 meta.data 中找到 *_res.* 聚类列，检查 graph.name 是否正确。")
}


# 8. 统计各分辨率簇大小模块

# 8.1 统计每个 resolution 下各 cluster 的细胞数
cluster_tables <- apply(
  md[, res_cols, drop = FALSE],
  2,
  table
)

# 8.2 输出各分辨率簇大小表
cat("\n各分辨率簇大小表：\n")
print(cluster_tables)


# 9. clustree 图模块

# 9.1 clustree 图保存参数
w_clustree <- 12
h_clustree <- 8
name_clustree <- "挑选分辨率"

# 9.2 生成 clustree 图
p_tree <- clustree(
  md,
  prefix = paste0(graph_name, "_res.")
)

# 9.3 显示 clustree 图
print(p_tree)

# 9.4 保存 clustree 图
pdf(
  file = file.path(out_dir, paste0(name_clustree, ".pdf")),
  width = w_clustree,
  height = h_clustree
)

print(p_tree)

dev.off()


# 10. 将聚类后的对象写回 seurat

# 10.1 写回当前 R 环境中的 seurat 对象
seurat <- srt


# 11. 参数记录模块

# 11.1 参数文件保存参数
name_params <- "resolution_choose_parameters"

# 11.2 生成参数记录文本
param_text <- paste0(
  "本次分析参数总结：\n",
  "- reduction：", reduction, "\n",
  "- dims：1 ~ ", pcs_end, "\n",
  "- graph.name：", graph_name, "\n",
  "- algorithm：", algorithm, "\n",
  "- resolutions：", paste(res_vec, collapse = ", "), "\n",
  "写作提示词（自行组装语言，或借助AI组装）：\n",
  "1.基于指定降维空间构建 KNN 图（FindNeighbors），并在多个分辨率下进行聚类（FindClusters）。\n",
  "2.通过 clustree 可视化不同分辨率下簇的分裂/合并关系，用于辅助确定最终分辨率。\n",
  "3.综合考虑簇数量变化幅度、簇稳定性及生物学可解释性，选择最终分辨率用于后续聚类与下游分析。"
)

# 11.3 保存参数记录
writeLines(
  param_text,
  con = file.path(out_dir, paste0(name_params, ".txt"))
)


# 12. 完成提示
cat("\n2.2 多分辨率聚类与 clustree 分析完成。\n")
cat("结果保存目录：", out_dir, "\n")
cat("clustree 图已保存：", file.path(out_dir, paste0(name_clustree, ".pdf")), "\n")
cat("参数文件已保存：", file.path(out_dir, paste0(name_params, ".txt")), "\n")
cat("当前 R 环境中的 seurat 对象已更新。\n")