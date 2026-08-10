library(Seurat)
library(dplyr)
library(ggplot2)

# 1. 检查全局环境中是否存在 seurat 对象
if (!exists("seurat", envir = .GlobalEnv)) {
  stop("全局环境中没有 seurat 对象")
}

srt <- get("seurat", envir = .GlobalEnv)


# 2. 输出目录模块

# 2.1 输出目录参数
out_dir <- "2.3分群聚类"

# 2.2 创建输出目录
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}


# 3. 聚类参数模块

# 3.1 聚类分辨率
resolution <- 0.1

# 3.2 降维方法
# 可选值需要存在于 Reductions(seurat) 中
# 常见可选值："pca"、"harmony"、"integrated.cca"、"integrated.rpca"、"integrated.mnn"、"integrated.Join"
reduction <- "pca"

# 3.3 PCA 维度
seuratPC <- 20

# 3.4 RunUMAP 参数
umap_n_neighbors <- 30
umap_min_dist <- 0.3
umap_spread <- 1

# 3.5 检查 reduction 是否存在
if (!reduction %in% Reductions(srt)) {
  stop(paste0("seurat 中不存在 reduction：", reduction))
}


# 4. FindAllMarkers 参数模块

# 4.1 差异基因参数
only_pos <- TRUE
min_pct <- 0.01
logfc_threshold <- 0.1

# 4.2 每个 cluster 选择 Top marker 数量
top_genes_per_cluster <- 10


# 5. FindClusters 聚类模块

# 5.1 按指定 resolution 进行聚类
srt <- FindClusters(
  srt,
  resolution = resolution
)

# 5.2 保存 cluster 信息到 meta.data
srt@meta.data$cluster <- as.character(Idents(srt))


# 6. FindAllMarkers 差异基因模块

# 6.1 计算每个 cluster 的 marker 基因
markers_all <- FindAllMarkers(
  srt,
  only.pos = only_pos,
  min.pct = min_pct,
  logfc.threshold = logfc_threshold
)

# 6.2 保存所有 marker 基因
write.csv(
  markers_all,
  file = file.path(out_dir, "每个聚类的marker基因.csv"),
  row.names = FALSE
)


# 7. 每个 cluster 导出 Top100 marker 模块

# 7.1 按 cluster 拆分 marker 表
split_tables <- markers_all %>%
  group_by(cluster) %>%
  arrange(cluster, desc(avg_log2FC)) %>%
  group_split()

# 7.2 初始化合并表
all_genes <- data.frame()

# 7.3 每个 cluster 保存前100个 marker
for (i in seq_along(split_tables)) {
  
  cluster_data <- split_tables[[i]]
  top_genes <- cluster_data[1:min(100, nrow(cluster_data)), ]
  
  all_genes <- rbind(all_genes, top_genes)
  
  file_name <- file.path(
    out_dir,
    paste0("cluster_", i - 1, "_sorted.csv")
  )
  
  write.csv(
    top_genes,
    file = file_name,
    row.names = FALSE
  )
}

# 7.4 保存所有 cluster 的 Top100 合并表
write.csv(
  all_genes,
  file = file.path(out_dir, "所有聚类_top100合并表.csv"),
  row.names = FALSE
)


# 8. 选择 Top markers 模块

# 8.1 每个 cluster 选择 Top marker
top_markers <- markers_all %>%
  group_by(cluster) %>%
  top_n(n = top_genes_per_cluster, wt = avg_log2FC)

# 8.2 查看 Top markers
cat("Top markers 数量：", nrow(top_markers), "\n")
print(head(top_markers, 50))


# 9. RunUMAP 和 RunTSNE 模块

# 9.1 运行 UMAP
srt <- RunUMAP(
  srt,
  dims = 1:seuratPC,
  reduction = reduction,
  n.neighbors = umap_n_neighbors,
  min.dist = umap_min_dist,
  spread = umap_spread,
  verbose = FALSE
)

# 9.2 运行 tSNE
srt <- RunTSNE(
  srt,
  dims = 1:seuratPC,
  reduction = reduction,
  check_duplicates = FALSE,
  verbose = FALSE
)


# 10. 热图模块

# 10.1 热图参数
heat_group_colors_text <- "#E64B35FF,#4DBBD5FF,#00A087FF,#3C5488FF,#F39B7FFF,#8491B4FF,#91D1C2FF,#7E6148FF,#0073C2FF,#EFC000FF,#868686FF,#CD534CFF,#7AA6DCFF,#003C67FF,#8F7700FF,#3B3B3BFF,#A73030FF,#374E55FF,#DF8F44FF,#00A1D5FF,#B24745FF,#79AF97FF,#6A6599FF,#80796BFF,#BC3C29FF,#0072B5FF,#E18727FF,#20854EFF,#7876B1FF,#6F99ADFF,#FFDC91FF,#EE4C97FF"

heat_group_colors <- unlist(strsplit(heat_group_colors_text, ","))
heat_group_colors <- trimws(heat_group_colors)
heat_group_colors <- heat_group_colors[heat_group_colors != ""]

heat_low <- "#2166AC"
heat_mid <- "#FFFFFF"
heat_high <- "#B2182B"
heat_name <- "Expression"

w_heatmap <- 12
h_heatmap <- 8
name_heatmap <- "1.聚类热图"

# 10.2 生成聚类热图
p_heatmap <- suppressWarnings(
  DoHeatmap(
    srt,
    features = top_markers$gene,
    group.colors = heat_group_colors
  )
)

p_heatmap <- p_heatmap +
  scale_fill_gradient2(
    low = heat_low,
    mid = heat_mid,
    high = heat_high,
    name = heat_name
  )

# 10.3 保存聚类热图
pdf(
  file = file.path(out_dir, paste0(name_heatmap, ".pdf")),
  width = w_heatmap,
  height = h_heatmap
)

print(p_heatmap)

dev.off()


# 11. UMAP 图模块

# 11.1 UMAP 图参数
pt.size_umap <- 0.5

palette_umap_text <- "#E64B35FF,#4DBBD5FF,#00A087FF,#3C5488FF,#F39B7FFF,#8491B4FF,#91D1C2FF,#7E6148FF,#0073C2FF,#EFC000FF,#868686FF,#CD534CFF,#7AA6DCFF,#003C67FF,#8F7700FF,#3B3B3BFF,#A73030FF,#374E55FF,#DF8F44FF,#00A1D5FF,#B24745FF,#79AF97FF,#6A6599FF,#80796BFF,#BC3C29FF,#0072B5FF,#E18727FF,#20854EFF,#7876B1FF,#6F99ADFF,#FFDC91FF,#EE4C97FF"

palette_umap <- unlist(strsplit(palette_umap_text, ","))
palette_umap <- trimws(palette_umap)
palette_umap <- palette_umap[palette_umap != ""]

w_umap <- 10
h_umap <- 7
name_umap <- "2.聚类后UMAP"

# 11.2 生成 cluster UMAP
p_umap_cluster <- suppressWarnings(
  DimPlot(
    srt,
    reduction = "umap",
    group.by = NULL,
    label = TRUE,
    label.size = 3.5,
    pt.size = pt.size_umap,
    cols = palette_umap,
    raster = FALSE
  ) +
    theme_classic() +
    theme(
      panel.border = element_rect(fill = NA, color = "black", size = 0.5),
      legend.position = "right"
    )
)

# 11.3 生成 group UMAP
if ("group" %in% colnames(srt[[]])) {
  
  p_umap_group <- suppressWarnings(
    DimPlot(
      srt,
      reduction = "umap",
      group.by = "group",
      label = FALSE,
      label.size = 3.5,
      pt.size = pt.size_umap,
      cols = palette_umap,
      raster = FALSE
    ) +
      theme_classic() +
      theme(
        panel.border = element_rect(fill = NA, color = "black", size = 0.5),
        legend.position = "right"
      )
  )
  
} else {
  
  p_umap_group <- ggplot() +
    theme_void() +
    ggtitle("meta.data 中不存在列：group（此图跳过）")
}

# 11.4 生成 orig.ident UMAP
if ("orig.ident" %in% colnames(srt[[]])) {
  
  p_umap_orig <- suppressWarnings(
    DimPlot(
      srt,
      reduction = "umap",
      group.by = "orig.ident",
      label = FALSE,
      label.size = 3.5,
      pt.size = pt.size_umap,
      cols = palette_umap,
      raster = FALSE
    ) +
      theme_classic() +
      theme(
        panel.border = element_rect(fill = NA, color = "black", size = 0.5),
        legend.position = "right"
      )
  )
  
} else {
  
  p_umap_orig <- ggplot() +
    theme_void() +
    ggtitle("meta.data 中不存在列：orig.ident（此图跳过）")
}

# 11.5 保存 UMAP 三图到一个 PDF
pdf(
  file = file.path(out_dir, paste0(name_umap, "_3plots.pdf")),
  width = w_umap,
  height = h_umap
)

print(p_umap_cluster)
print(p_umap_group)
print(p_umap_orig)

dev.off()


# 12. tSNE 图模块

# 12.1 tSNE 图参数
pt.size_tsne <- 0.5

palette_tsne_text <- "#E64B35FF,#4DBBD5FF,#00A087FF,#3C5488FF,#F39B7FFF,#8491B4FF,#91D1C2FF,#7E6148FF,#0073C2FF,#EFC000FF,#868686FF,#CD534CFF,#7AA6DCFF,#003C67FF,#8F7700FF,#3B3B3BFF,#A73030FF,#374E55FF,#DF8F44FF,#00A1D5FF,#B24745FF,#79AF97FF,#6A6599FF,#80796BFF,#BC3C29FF,#0072B5FF,#E18727FF,#20854EFF,#7876B1FF,#6F99ADFF,#FFDC91FF,#EE4C97FF"

palette_tsne <- unlist(strsplit(palette_tsne_text, ","))
palette_tsne <- trimws(palette_tsne)
palette_tsne <- palette_tsne[palette_tsne != ""]

w_tsne <- 10
h_tsne <- 7
name_tsne <- "3.聚类后TSNE"

# 12.2 生成 cluster tSNE
p_tsne_cluster <- suppressWarnings(
  DimPlot(
    srt,
    reduction = "tsne",
    group.by = NULL,
    label = TRUE,
    label.size = 3.5,
    pt.size = pt.size_tsne,
    cols = palette_tsne,
    raster = FALSE
  ) +
    theme_classic() +
    theme(
      panel.border = element_rect(fill = NA, color = "black", size = 0.5),
      legend.position = "right"
    )
)

# 12.3 生成 group tSNE
if ("group" %in% colnames(srt[[]])) {
  
  p_tsne_group <- suppressWarnings(
    DimPlot(
      srt,
      reduction = "tsne",
      group.by = "group",
      label = FALSE,
      label.size = 3.5,
      pt.size = pt.size_tsne,
      cols = palette_tsne,
      raster = FALSE
    ) +
      theme_classic() +
      theme(
        panel.border = element_rect(fill = NA, color = "black", size = 0.5),
        legend.position = "right"
      )
  )
  
} else {
  
  p_tsne_group <- ggplot() +
    theme_void() +
    ggtitle("meta.data 中不存在列：group（此图跳过）")
}

# 12.4 生成 orig.ident tSNE
if ("orig.ident" %in% colnames(srt[[]])) {
  
  p_tsne_orig <- suppressWarnings(
    DimPlot(
      srt,
      reduction = "tsne",
      group.by = "orig.ident",
      label = FALSE,
      label.size = 3.5,
      pt.size = pt.size_tsne,
      cols = palette_tsne,
      raster = FALSE
    ) +
      theme_classic() +
      theme(
        panel.border = element_rect(fill = NA, color = "black", size = 0.5),
        legend.position = "right"
      )
  )
  
} else {
  
  p_tsne_orig <- ggplot() +
    theme_void() +
    ggtitle("meta.data 中不存在列：orig.ident（此图跳过）")
}

# 12.5 保存 tSNE 三图到一个 PDF
pdf(
  file = file.path(out_dir, paste0(name_tsne, "_3plots.pdf")),
  width = w_tsne,
  height = h_tsne
)

print(p_tsne_cluster)
print(p_tsne_group)
print(p_tsne_orig)

dev.off()


# 13. 参数记录模块

# 13.1 参数文件保存参数
name_params <- "clustering_parameters"

# 13.2 生成参数记录文本
param_text <- paste0(
  "本次分析参数总结：\n",
  "- resolution：", resolution, "\n",
  "- RunUMAP/RunTSNE reduction：", reduction, "\n",
  "- PCA dims：1 ~ ", seuratPC, "\n",
  "- RunUMAP n.neighbors：", umap_n_neighbors, "\n",
  "- RunUMAP min.dist：", umap_min_dist, "\n",
  "- RunUMAP spread：", umap_spread, "\n",
  "- FindAllMarkers only.pos：", only_pos, "\n",
  "- FindAllMarkers min.pct：", min_pct, "\n",
  "- FindAllMarkers logfc.threshold：", logfc_threshold, "\n",
  "- 每cluster Top markers：", top_genes_per_cluster, "\n",
  "写作提示词（自行组装语言，或借助AI组装）：\n",
  "1.基于指定分辨率进行细胞聚类（FindClusters），并对各聚类进行差异基因筛选（FindAllMarkers）。\n",
  "2.结合每个聚类的Top marker基因热图与降维可视化结果（UMAP/tSNE）评估聚类合理性。\n",
  "3.后续可在此基础上进行细胞类型注释与下游功能分析。"
)

# 13.3 保存参数记录
writeLines(
  param_text,
  con = file.path(out_dir, paste0(name_params, ".txt"))
)


# 14. 将结果写回 seurat

# 14.1 更新当前 R 环境中的 seurat 对象
seurat <- srt


# 15. 完成提示
cat("\n2.3 分群聚类分析完成。\n")
cat("结果保存目录：", out_dir, "\n")
cat("marker 总表已保存：", file.path(out_dir, "每个聚类的marker基因.csv"), "\n")
cat("聚类热图已保存：", file.path(out_dir, paste0(name_heatmap, ".pdf")), "\n")
cat("UMAP 三图已保存：", file.path(out_dir, paste0(name_umap, "_3plots.pdf")), "\n")
cat("tSNE 三图已保存：", file.path(out_dir, paste0(name_tsne, "_3plots.pdf")), "\n")
cat("参数文件已保存：", file.path(out_dir, paste0(name_params, ".txt")), "\n")
cat("当前 R 环境中的 seurat 对象已更新。\n")