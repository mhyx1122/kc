# 1. 加载必要的 R 包

suppressPackageStartupMessages({
  library(Seurat)
  library(SCENIC)
  library(AUCell)
  library(ggplot2)
  library(reshape2)
  library(pheatmap)
})


# 2. 检查分析所需对象

if (!exists("seurat", envir = .GlobalEnv)) {
  stop("全局环境中没有找到 seurat 对象。")
}

if (!exists("cellInfo", envir = .GlobalEnv)) {
  stop("全局环境中没有找到 cellInfo 对象，请先运行 SCENIC 主分析代码。")
}

if (!exists("scenicOptions", envir = .GlobalEnv)) {
  stop("全局环境中没有找到 scenicOptions 对象，请先运行 SCENIC 主分析代码。")
}


# 3. 创建结果保存文件夹

output_dir <- "SCENIC"

if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}


# 4. 读取并保存 Regulon AUC 数据

regulon_auc_file <- "./int/3.4_regulonAUC.Rds"

if (!file.exists(regulon_auc_file)) {
  stop(
    paste0(
      "没有找到 Regulon AUC 文件：",
      regulon_auc_file
    )
  )
}

regulonAUC <- readRDS(
  regulon_auc_file
)

regulonAUC <- as.matrix(
  regulonAUC@assays@data@listData[["AUC"]]
)

write.csv(
  regulonAUC,
  file = file.path(
    output_dir,
    "regulonAUC.csv"
  )
)


# 5. 创建 Regulon AUC Seurat 对象

regulonAUC_seurat <- CreateSeuratObject(
  counts = regulonAUC,
  project = "regulonAUC"
)

regulonAUC_seurat$celltype <- cellInfo$CellType

regulonAUC_seurat <- NormalizeData(
  object = regulonAUC_seurat,
  verbose = FALSE
)

regulonAUC_seurat <- FindVariableFeatures(
  object = regulonAUC_seurat,
  selection.method = "vst",
  nfeatures = 20
)

regulonAUC_seurat <- ScaleData(
  regulonAUC_seurat
)

regulonAUC_seurat <- RunPCA(
  regulonAUC_seurat,
  features = VariableFeatures(
    object = regulonAUC_seurat
  )
)

regulonAUC_seurat <- RunUMAP(
  regulonAUC_seurat,
  dims = 1:3
)


# 6. 使用原 Seurat 对象中的 UMAP 坐标

if (!"umap" %in% names(seurat@reductions)) {
  stop("seurat 对象中不存在 umap 降维结果。")
}

umap_coordinates <- seurat@reductions[["umap"]]@cell.embeddings

sampled_cells <- Cells(
  regulonAUC_seurat
)

umap_coordinates_sampled <- umap_coordinates[
  sampled_cells,
  ,
  drop = FALSE
]

regulonAUC_seurat@reductions[["umap"]]@cell.embeddings <-
  umap_coordinates_sampled


# 7. 绘制全部 Regulon 的 FeaturePlot

plot_width <- 8
plot_height <- 6

featureplot_colors <- c(
  "#440154",
  "#3B528B",
  "#21908C",
  "#5DC863",
  "#FDE725"
)

pdf(
  file = file.path(
    output_dir,
    "1.FeaturePlot_regulonAUC_umap.pdf"
  ),
  width = plot_width,
  height = plot_height
)

print(
  FeaturePlot(
    regulonAUC_seurat,
    features = rownames(regulonAUC_seurat),
    cols = featureplot_colors,
    reduction = "umap"
  )
)

dev.off()


# 8. 分别保存每个 Regulon 的 FeaturePlot

for (feature in rownames(regulonAUC_seurat)) {
  
  pdf_filename <- file.path(
    output_dir,
    paste0(
      "FeaturePlot_regulonAUC_umap_",
      feature,
      ".pdf"
    )
  )
  
  pdf(
    file = pdf_filename,
    width = plot_width,
    height = plot_height
  )
  
  print(
    FeaturePlot(
      object = regulonAUC_seurat,
      features = feature,
      cols = featureplot_colors,
      reduction = "umap"
    )
  )
  
  dev.off()
}


# 9. 读取 Regulon AUC 并去除重复 Regulon

ifuni4 <- TRUE

regulonAUC <- loadInt(
  scenicOptions,
  "aucell_regulonAUC"
)

if (ifuni4) {
  
  regulonAUC <- regulonAUC[
    onlyNonDuplicatedExtended(
      rownames(regulonAUC)
    ),
    ,
    drop = FALSE
  ]
  
} else {
  
  message("跳过去除重复 Regulon 的操作。")
}


# 10. 计算每种细胞类型的平均 Regulon 活性

regulonActivity_byCellType <- sapply(
  split(
    rownames(cellInfo),
    cellInfo$CellType
  ),
  function(cells) {
    rowMeans(
      getAUC(regulonAUC)[
        ,
        cells,
        drop = FALSE
      ]
    )
  }
)

write.csv(
  regulonActivity_byCellType,
  file = file.path(
    output_dir,
    "regulonActivity_byCellType.csv"
  )
)


# 11. 绘制细胞类型平均 Regulon 活性热图

heatmap_colors <- c(
  "blue",
  "white",
  "red"
)

regulonActivity_byCellType_Scaled <- t(
  scale(
    t(regulonActivity_byCellType),
    center = TRUE,
    scale = TRUE
  )
)

pdf(
  file = file.path(
    output_dir,
    "2.regulonActivity_byCellType_Scaled_heatmap.pdf"
  ),
  width = plot_width,
  height = plot_height
)

pheatmap::pheatmap(
  regulonActivity_byCellType_Scaled,
  fontsize_row = 8,
  color = colorRampPalette(
    heatmap_colors
  )(100),
  breaks = seq(
    -3,
    3,
    length.out = 100
  ),
  treeheight_row = 10,
  treeheight_col = 10,
  border_color = NA
)

dev.off()


# 12. 绘制全部细胞的 Regulon 活性热图

regulonAUC_matrix <- as.matrix(
  regulonAUC@assays@data@listData$AUC
)

pdf(
  file = file.path(
    output_dir,
    "2.regulonActivity_heatmap.pdf"
  ),
  width = plot_width,
  height = plot_height
)

pheatmap::pheatmap(
  regulonAUC_matrix,
  fontsize_row = 8,
  color = colorRampPalette(
    heatmap_colors
  )(100),
  breaks = seq(
    -3,
    3,
    length.out = 100
  ),
  treeheight_row = 10,
  treeheight_col = 10,
  border_color = NA,
  annotation_col = cellInfo,
  scale = "row",
  show_rownames = TRUE,
  show_colnames = FALSE,
  cluster_cols = TRUE
)

dev.off()


# 13. 提取各细胞类型中的高活性 Regulon

topRegulators <- reshape2::melt(
  apply(
    regulonActivity_byCellType_Scaled,
    2,
    function(x) {
      cbind(
        sort(
          x[x > 0],
          decreasing = TRUE
        )
      )
    }
  )
)[
  c(
    "L1",
    "Var1",
    "value"
  )
]

colnames(topRegulators) <- c(
  "CellType",
  "Regulon",
  "RelativeActivity"
)

write.csv(
  topRegulators,
  file = file.path(
    output_dir,
    "celltype_topRegulators.csv"
  )
)


# 14. 读取并绘制二值化 Regulon 活性热图

binary_colors <- c(
  "white",
  "red"
)

binary_regulon_file <- "./int/4.1_binaryRegulonActivity.Rds"

if (!file.exists(binary_regulon_file)) {
  stop(
    paste0(
      "没有找到二值化 Regulon 活性文件：",
      binary_regulon_file
    )
  )
}

binary <- readRDS(
  binary_regulon_file
)

pdf(
  file = file.path(
    output_dir,
    "3.binary_regulon_activity_heatmap.pdf"
  ),
  width = plot_width,
  height = plot_height
)

pheatmap::pheatmap(
  binary,
  color = colorRampPalette(
    colors = binary_colors
  )(100),
  breaks = unique(
    seq(
      0,
      1,
      length = 100
    )
  ),
  annotation_col = cellInfo,
  scale = "row",
  show_rownames = TRUE,
  show_colnames = FALSE
)

dev.off()


# 15. 按细胞类型计算二值化 Regulon 活性比例

minPerc4 <- 0.5

binarized_colors <- c(
  "white",
  "pink",
  "red"
)

binaryRegulonActivity <- loadInt(
  scenicOptions,
  "aucell_binary_nonDupl"
)

cellInfo_binarizedCells <- cellInfo[
  which(
    rownames(cellInfo) %in%
      colnames(binaryRegulonActivity)
  ),
  ,
  drop = FALSE
]

regulonActivity_byCellType_Binarized <- sapply(
  split(
    rownames(cellInfo_binarizedCells),
    cellInfo_binarizedCells$CellType
  ),
  function(cells) {
    rowMeans(
      binaryRegulonActivity[
        ,
        cells,
        drop = FALSE
      ]
    )
  }
)

binaryActPerc_subset <- regulonActivity_byCellType_Binarized[
  which(
    rowSums(
      regulonActivity_byCellType_Binarized >
        minPerc4
    ) > 0
  ),
  ,
  drop = FALSE
]


# 16. 绘制细胞类型二值化 Regulon 活性热图

pdf(
  file = file.path(
    output_dir,
    "4.binary_regulon_activity_binarized_heatmap.pdf"
  ),
  width = plot_width,
  height = plot_height
)

pheatmap::pheatmap(
  binaryActPerc_subset,
  fontsize_row = 8,
  color = binarized_colors
)

dev.off()


# 17. 导出高比例二值化 Regulon

topRegulators_binary <- reshape2::melt(
  regulonActivity_byCellType_Binarized
)

colnames(topRegulators_binary) <- c(
  "Regulon",
  "CellType",
  "RelativeActivity"
)

topRegulators_binary <- topRegulators_binary[
  which(
    topRegulators_binary$RelativeActivity >
      minPerc4
  ),
  ,
  drop = FALSE
]

write.csv(
  topRegulators_binary,
  file = file.path(
    output_dir,
    "celltype_topRegulators_binary.csv"
  )
)


# 18. 计算 Regulon 特异性得分

rss_colors <- c(
  "grey90",
  "darkolivegreen3",
  "darkgreen"
)

rss <- calcRSS(
  AUC = getAUC(regulonAUC),
  cellAnnotation = cellInfo[
    colnames(regulonAUC),
    "CellType"
  ]
)

rssPlot <- plotRSS(
  rss,
  zThreshold = NULL,
  col.low = rss_colors[1],
  col.mid = rss_colors[2],
  col.high = rss_colors[3]
)


# 19. 保存全部细胞类型 RSS 图

ggsave(
  filename = file.path(
    output_dir,
    "5.rss_plot.pdf"
  ),
  plot = rssPlot[[1]],
  width = plot_width,
  height = plot_height
)


# 20. 绘制并保存指定细胞类型的 Regulon

specific_cell_type <- unique(
  cellInfo$CellType
)[1]

plotRSS_oneSet(
  rss,
  setName = specific_cell_type
)

ggsave(
  filename = file.path(
    output_dir,
    "5.specific_regulators_plot.pdf"
  ),
  width = plot_width,
  height = plot_height
)