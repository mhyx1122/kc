suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(stringr)
})

# 1. 检查全局环境中是否存在 seurat 对象
if (!exists("seurat", envir = .GlobalEnv)) {
  stop("全局环境中没有 seurat 对象")
}

srt <- get("seurat", envir = .GlobalEnv)


# 2. 细胞周期分析参数模块

# 2.1 输出目录参数
out_dir_cycle <- "1.4细胞周期分析"

if (!dir.exists(out_dir_cycle)) {
  dir.create(out_dir_cycle, recursive = TRUE)
}

# 2.2 物种参数
# 可选值："human" 或 "mouse"
scRNA_species <- "human"

# 2.3 PCA 图保存参数
pt_size_cycle_pca <- 0.5
w_cycle_pca <- 7.5
h_cycle_pca <- 5.5
name_cycle_pca <- "1.PCA图（按细胞周期分组）"

palette_cycle_text <- "#E64B35FF,#4DBBD5FF,#00A087FF,#3C5488FF,#F39B7FFF,#8491B4FF,#91D1C2FF,#7E6148FF,#0073C2FF,#EFC000FF,#868686FF,#CD534CFF,#7AA6DCFF,#003C67FF,#8F7700FF,#3B3B3BFF,#A73030FF,#374E55FF,#DF8F44FF,#00A1D5FF,#B24745FF,#79AF97FF,#6A6599FF,#80796BFF,#BC3C29FF,#0072B5FF,#E18727FF,#20854EFF,#7876B1FF,#6F99ADFF,#FFDC91FF,#EE4C97FF"

palette_cycle <- unlist(strsplit(palette_cycle_text, ","))
palette_cycle <- trimws(palette_cycle)
palette_cycle <- palette_cycle[palette_cycle != ""]

# 2.4 UMAP 图保存参数
pt_size_cycle_umap <- 0.5
w_cycle_umap <- 7.5
h_cycle_umap <- 5.5
name_cycle_umap <- "2.UMAP图（按细胞周期分组）"

palette_cycle_umap_text <- "#E64B35FF,#4DBBD5FF,#00A087FF,#3C5488FF,#F39B7FFF,#8491B4FF,#91D1C2FF,#7E6148FF,#0073C2FF,#EFC000FF,#868686FF,#CD534CFF,#7AA6DCFF,#003C67FF,#8F7700FF,#3B3B3BFF,#A73030FF,#374E55FF,#DF8F44FF,#00A1D5FF,#B24745FF,#79AF97FF,#6A6599FF,#80796BFF,#BC3C29FF,#0072B5FF,#E18727FF,#20854EFF,#7876B1FF,#6F99ADFF,#FFDC91FF,#EE4C97FF"

palette_cycle_umap <- unlist(strsplit(palette_cycle_umap_text, ","))
palette_cycle_umap <- trimws(palette_cycle_umap)
palette_cycle_umap <- palette_cycle_umap[palette_cycle_umap != ""]


# 3. 加载细胞周期基因集模块

if (scRNA_species == "human") {
  
  s.genes <- Seurat::cc.genes.updated.2019$s.genes
  g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes
  
} else if (scRNA_species == "mouse") {
  
  s.genes <- Seurat::cc.genes.updated.2019$s.genes
  g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes
  
  s.genes <- stringr::str_to_title(s.genes)
  g2m.genes <- stringr::str_to_title(g2m.genes)
  
} else {
  
  stop("无法识别的物种类型，请设置 scRNA_species 为 human 或 mouse")
}


# 4. 细胞周期评分模块

# 4.1 进行细胞周期评分
srt <- CellCycleScoring(
  srt,
  s.features = s.genes,
  g2m.features = g2m.genes,
  set.ident = TRUE
)

# 4.2 写回当前 R 环境中的 seurat 对象
seurat <- srt


# 5. 保存细胞周期 PCA 图模块

# 5.1 检查 PCA 是否存在
if (!"pca" %in% Reductions(srt)) {
  stop("seurat 中没有 reduction pca，请先 RunPCA")
}

# 5.2 生成 PCA 图
p_cycle_pca <- suppressWarnings(
  DimPlot(
    srt,
    reduction = "pca",
    group.by = "Phase",
    pt.size = pt_size_cycle_pca,
    raster = FALSE,
    cols = palette_cycle
  ) +
    theme_classic() +
    theme(
      panel.border = element_rect(fill = NA, color = "black", size = 0.5),
      legend.position = "right"
    )
)

# 5.3 保存 PCA 图
ggsave(
  filename = file.path(out_dir_cycle, paste0(name_cycle_pca, ".pdf")),
  plot = p_cycle_pca,
  width = w_cycle_pca,
  height = h_cycle_pca,
  device = "pdf"
)


# 6. 保存细胞周期 UMAP 图模块

# 6.1 检查 UMAP 是否存在
if (!"umap" %in% Reductions(srt)) {
  stop("seurat 中没有 reduction umap，请先 RunUMAP")
}

# 6.2 生成 UMAP 图
p_cycle_umap <- suppressWarnings(
  DimPlot(
    srt,
    reduction = "umap",
    group.by = "Phase",
    pt.size = pt_size_cycle_umap,
    raster = FALSE,
    cols = palette_cycle_umap
  ) +
    theme_classic() +
    theme(
      panel.border = element_rect(fill = NA, color = "black", size = 0.5),
      legend.position = "right"
    )
)

# 6.3 保存 UMAP 图
ggsave(
  filename = file.path(out_dir_cycle, paste0(name_cycle_umap, ".pdf")),
  plot = p_cycle_umap,
  width = w_cycle_umap,
  height = h_cycle_umap,
  device = "pdf"
)


# 7. 细胞周期效应去除参数模块

# 7.1 输出目录参数
out_dir_regress <- "1.5细胞周期校正后降维作图"

if (!dir.exists(out_dir_regress)) {
  dir.create(out_dir_regress, recursive = TRUE)
}

# 7.2 是否执行去除细胞周期效应
run_regress_step <- TRUE

# 7.3 去除后重新筛选高变基因参数
hvg_nfeatures <- 2000
hvg_top_n <- 10

# 7.4 去除后 UMAP 参数
umap_dims_max <- 20

# 7.5 重新筛选高变基因图保存参数
w_hvg_regress <- 7
h_hvg_regress <- 7
name_hvg_regress <- "1.重新筛选高变基因"

# 7.6 去除后 UMAP 图保存参数
pt_size_regress_umap <- 0.5
w_regress_umap <- 7.5
h_regress_umap <- 5.5
name_regress_umap <- "降维umap图（按细胞周期分组）"

palette_regress_umap_text <- "#E64B35FF,#4DBBD5FF,#00A087FF,#3C5488FF,#F39B7FFF,#8491B4FF,#91D1C2FF,#7E6148FF,#0073C2FF,#EFC000FF,#868686FF,#CD534CFF,#7AA6DCFF,#003C67FF,#8F7700FF,#3B3B3BFF,#A73030FF,#374E55FF,#DF8F44FF,#00A1D5FF,#B24745FF,#79AF97FF,#6A6599FF,#80796BFF,#BC3C29FF,#0072B5FF,#E18727FF,#20854EFF,#7876B1FF,#6F99ADFF,#FFDC91FF,#EE4C97FF"

palette_regress_umap <- unlist(strsplit(palette_regress_umap_text, ","))
palette_regress_umap <- trimws(palette_regress_umap)
palette_regress_umap <- palette_regress_umap[palette_regress_umap != ""]

# 7.7 去除后 PCA 图保存参数
pt_size_regress_pca <- 0.5
w_regress_pca <- 7.5
h_regress_pca <- 5.5
name_regress_pca <- "降维pca图（按细胞周期分组）"

palette_regress_pca_text <- "#E64B35FF,#4DBBD5FF,#00A087FF,#3C5488FF,#F39B7FFF,#8491B4FF,#91D1C2FF,#7E6148FF,#0073C2FF,#EFC000FF,#868686FF,#CD534CFF,#7AA6DCFF,#003C67FF,#8F7700FF,#3B3B3BFF,#A73030FF,#374E55FF,#DF8F44FF,#00A1D5FF,#B24745FF,#79AF97FF,#6A6599FF,#80796BFF,#BC3C29FF,#0072B5FF,#E18727FF,#20854EFF,#7876B1FF,#6F99ADFF,#FFDC91FF,#EE4C97FF"

palette_regress_pca <- unlist(strsplit(palette_regress_pca_text, ","))
palette_regress_pca <- trimws(palette_regress_pca)
palette_regress_pca <- palette_regress_pca[palette_regress_pca != ""]


# 8. 去除细胞周期效应模块

if (run_regress_step) {
  
  # 8.1 检查细胞周期评分列是否存在
  if (!all(c("S.Score", "G2M.Score", "Phase") %in% colnames(srt@meta.data))) {
    stop("seurat@meta.data 中缺少 S.Score、G2M.Score 或 Phase，请先完成 CellCycleScoring")
  }
  
  # 8.2 ScaleData 回归细胞周期效应
  all_features <- rownames(srt)
  
  srt <- ScaleData(
    srt,
    features = all_features,
    vars.to.regress = c("S.Score", "G2M.Score")
  )
  
  # 8.3 重新筛选高变基因
  srt <- FindVariableFeatures(
    srt,
    selection.method = "vst",
    nfeatures = hvg_nfeatures
  )
  
  # 8.4 生成重新筛选高变基因图
  top_hvg <- head(VariableFeatures(srt), hvg_top_n)
  
  p_hvg_regress <- VariableFeaturePlot(srt)
  
  p_hvg_regress <- LabelPoints(
    plot = p_hvg_regress,
    points = top_hvg,
    repel = TRUE,
    xnudge = 0,
    ynudge = 0
  )
  
  # 8.5 保存重新筛选高变基因图
  ggsave(
    filename = file.path(out_dir_regress, paste0(name_hvg_regress, ".pdf")),
    plot = p_hvg_regress,
    width = w_hvg_regress,
    height = h_hvg_regress,
    device = "pdf"
  )
  
  # 8.6 重新运行 PCA
  srt <- RunPCA(
    srt,
    features = VariableFeatures(object = srt)
  )
  
  # 8.7 重新运行 UMAP
  srt <- RunUMAP(
    srt,
    dims = 1:umap_dims_max
  )
  
  # 8.8 生成去除后 UMAP 图
  p_regress_umap <- suppressWarnings(
    DimPlot(
      srt,
      reduction = "umap",
      group.by = "Phase",
      pt.size = pt_size_regress_umap,
      raster = FALSE,
      cols = palette_regress_umap
    ) +
      theme_classic() +
      theme(
        panel.border = element_rect(fill = NA, color = "black", size = 0.5),
        legend.position = "right"
      )
  )
  
  # 8.9 保存去除后 UMAP 图
  ggsave(
    filename = file.path(out_dir_regress, paste0(name_regress_umap, ".pdf")),
    plot = p_regress_umap,
    width = w_regress_umap,
    height = h_regress_umap,
    device = "pdf"
  )
  
  # 8.10 生成去除后 PCA 图
  p_regress_pca <- suppressWarnings(
    DimPlot(
      srt,
      reduction = "pca",
      group.by = "Phase",
      pt.size = pt_size_regress_pca,
      raster = FALSE,
      cols = palette_regress_pca
    ) +
      theme_classic() +
      theme(
        panel.border = element_rect(fill = NA, color = "black", size = 0.5),
        legend.position = "right"
      )
  )
  
  # 8.11 保存去除后 PCA 图
  ggsave(
    filename = file.path(out_dir_regress, paste0(name_regress_pca, ".pdf")),
    plot = p_regress_pca,
    width = w_regress_pca,
    height = h_regress_pca,
    device = "pdf"
  )
  
  # 8.12 写回当前 R 环境中的 seurat 对象
  seurat <- srt
}


# 9. 参数记录模块

# 9.1 参数文件保存参数
name_params <- "cell_cycle_parameters"

# 9.2 生成参数记录文本
param_text <- paste0(
  "本次分析参数总结：\n",
  "- 物种：", scRNA_species, "\n",
  "- 第一步：进行细胞周期评分，并生成 PCA / UMAP 按细胞周期分组图\n",
  "- 细胞周期评分函数：CellCycleScoring\n",
  "- set.ident：TRUE\n",
  "- 是否执行去除细胞周期效应：", run_regress_step, "\n",
  "- 去除细胞周期效应时回归变量：S.Score, G2M.Score\n",
  "- 去除后重新筛选高变基因数：", hvg_nfeatures, "\n",
  "- 去除后标注前Top基因数：", hvg_top_n, "\n",
  "- 去除后 UMAP dims：1:", umap_dims_max, "\n"
)

# 9.3 保存参数记录
writeLines(
  param_text,
  con = file.path(out_dir_regress, paste0(name_params, ".txt"))
)