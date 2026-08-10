library(Seurat)
library(ggplot2)
library(plyr)
library(SCpubr)

# 1. 检查全局环境中是否存在 seurat 对象
if (!exists("seurat", envir = .GlobalEnv)) {
  stop("全局环境中没有 seurat 对象")
}

srt <- get("seurat", envir = .GlobalEnv)

if (!"cluster" %in% colnames(srt@meta.data)) {
  stop("seurat@meta.data 中没有 cluster 列")
}

if (!exists("do_DotPlot")) {
  stop("当前 R 环境中没有 do_DotPlot 函数，请确认 SCpubr 包是否正常加载")
}


# 2. 输出目录模块

# 2.1 输出目录参数
out_dir <- "3.2手动注释结果"

# 2.2 创建输出目录
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}


# 3. marker CSV 读取模块

# 3.1 marker CSV 文件路径
# CSV 格式要求：
# 每一列是一个细胞类型或分组名称
# 每一列下面填写对应 marker 基因
marker_csv_file <- "marker.csv"

# 3.2 读取 marker CSV
marker_df <- read.csv(
  marker_csv_file,
  header = TRUE,
  check.names = FALSE
)

# 3.3 整理 marker 基因列表
features_list <- list()

for (nm in names(marker_df)) {
  genes <- as.character(marker_df[[nm]])
  genes <- genes[genes != "" & !is.na(genes)]
  features_list[[nm]] <- genes
}


# 4. 手动注释模块

# 4.1 获取 cluster 顺序
cluster_levels <- sort(unique(as.character(srt@meta.data$cluster)))

if (all(!is.na(suppressWarnings(as.numeric(cluster_levels))))) {
  cluster_levels <- as.character(sort(as.numeric(cluster_levels)))
}

# 4.2 设置每个 cluster 对应的细胞类型名称
# cellnames 的顺序必须和 cluster_levels 的顺序一一对应
# 例如 cluster_levels 是 0、1、2、3，则 cellnames 第1个对应 cluster 0，第2个对应 cluster 1

cellnames <- c(
  "T cells",
  "B cells",
  "Macrophages",
  "Epithelial cells",
  "Fibroblasts",
  "Endothelial cells",
  "NK cells",
  "Monocytes",
  "Plasma cells",
  "Mast cells",
  "Dendritic cells",
  "Neutrophils",
  "Smooth muscle cells"
)

if (length(cellnames) != length(cluster_levels)) {
  stop("cellnames 的数量必须和 cluster_levels 的数量一致")
}

cluster_annotations <- setNames(cellnames, cluster_levels)

# 4.3 去掉原来的 cellType 列
srt2 <- srt

if ("cellType" %in% colnames(srt2@meta.data)) {
  srt2@meta.data$cellType <- NULL
}

# 4.4 按 cluster 设置 Idents
Idents(srt2) <- "cluster"

cur_lv <- levels(Idents(srt2))

if (all(!is.na(suppressWarnings(as.numeric(cur_lv))))) {
  Idents(srt2) <- factor(
    Idents(srt2),
    levels = as.character(sort(as.numeric(cur_lv)))
  )
}

# 4.5 将 cluster 映射为手动注释的 cellType
mapped <- plyr::mapvalues(
  x = as.character(Idents(srt2)),
  from = names(cluster_annotations),
  to = as.character(cluster_annotations),
  warn_missing = FALSE
)

Idents(srt2) <- mapped
srt2$cellType <- Idents(srt2)

# 4.6 写回当前 R 环境中的 seurat 对象
seurat <- srt2


# 5. 注释后 DotPlot 模块

# 5.1 DotPlot 绘图参数
dotplot_dot_scale <- 12
dotplot_legend_framewidth <- 2
dotplot_font_size <- 20

# 5.2 DotPlot 保存参数
dotplot_width <- 30
dotplot_height <- 15
dotplot_filename <- "第一次注释气泡结果图"

# 5.3 生成 DotPlot
p_dot <- do_DotPlot(
  sample = srt2,
  features = features_list,
  dot.scale = dotplot_dot_scale,
  legend.framewidth = dotplot_legend_framewidth,
  font.size = dotplot_font_size
)

# 5.4 保存 DotPlot
pdf(
  file = file.path(out_dir, paste0(dotplot_filename, ".pdf")),
  width = dotplot_width,
  height = dotplot_height
)

print(p_dot)

dev.off()


# 6. UMAP 分类图模块

# 6.1 UMAP 绘图参数
umap_width <- 8
umap_height <- 6
umap_label_size <- 3.5
umap_pt_size <- 0.5
raster_umap <- FALSE

umap_cols_text <- "#E64B35FF,#4DBBD5FF,#00A087FF,#3C5488FF,#F39B7FFF,#8491B4FF,#91D1C2FF,#7E6148FF,#0073C2FF,#EFC000FF,#868686FF,#CD534CFF,#7AA6DCFF,#003C67FF,#8F7700FF,#3B3B3BFF,#A73030FF,#374E55FF,#DF8F44FF,#00A1D5FF,#B24745FF,#79AF97FF,#6A6599FF,#80796BFF,#BC3C29FF,#0072B5FF,#E18727FF,#20854EFF,#7876B1FF,#6F99ADFF,#FFDC91FF,#EE4C97FF"

umap_cols <- trimws(unlist(strsplit(umap_cols_text, ",")))
umap_cols <- umap_cols[umap_cols != ""]

# 6.2 检查 UMAP 是否存在
if (!"umap" %in% Reductions(srt2)) {
  stop("seurat 中没有 reduction umap")
}

# 6.3 按细胞类型给色的 UMAP
p_umap_cell <- DimPlot(
  srt2,
  reduction = "umap",
  label = TRUE,
  cols = umap_cols,
  label.size = umap_label_size,
  pt.size = umap_pt_size,
  raster = raster_umap
) +
  theme_classic() +
  theme(
    panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"),
    legend.position = "right"
  )

pdf(
  file = file.path(out_dir, "分类UMAP（按细胞给色）.pdf"),
  width = umap_width,
  height = umap_height
)

print(p_umap_cell)

dev.off()

# 6.4 按 group 给色的 UMAP
if ("group" %in% colnames(srt2@meta.data)) {
  
  p_umap_group <- DimPlot(
    srt2,
    reduction = "umap",
    label = TRUE,
    group.by = "group",
    cols = umap_cols,
    label.size = umap_label_size,
    pt.size = umap_pt_size,
    raster = raster_umap
  ) +
    theme_classic() +
    theme(
      panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"),
      legend.position = "right"
    )
  
  pdf(
    file = file.path(out_dir, "分类UMAP（按分组给色）.pdf"),
    width = umap_width,
    height = umap_height
  )
  
  print(p_umap_group)
  
  dev.off()
}

# 6.5 按样本给色的 UMAP
if ("orig.ident" %in% colnames(srt2@meta.data)) {
  
  p_umap_sample <- DimPlot(
    srt2,
    reduction = "umap",
    label = TRUE,
    group.by = "orig.ident",
    cols = umap_cols,
    label.size = umap_label_size,
    pt.size = umap_pt_size,
    raster = raster_umap
  ) +
    theme_classic() +
    theme(
      panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"),
      legend.position = "right"
    )
  
  pdf(
    file = file.path(out_dir, "分类UMAP（按样本给色）.pdf"),
    width = umap_width,
    height = umap_height
  )
  
  print(p_umap_sample)
  
  dev.off()
}


# 7. tSNE 分类图模块

# 7.1 tSNE 绘图参数
tsne_width <- 8
tsne_height <- 6
tsne_label_size <- 3.5
tsne_pt_size <- 0.5
raster_tsne <- FALSE

tsne_cols_text <- "#E64B35FF,#4DBBD5FF,#00A087FF,#3C5488FF,#F39B7FFF,#8491B4FF,#91D1C2FF,#7E6148FF,#0073C2FF,#EFC000FF,#868686FF,#CD534CFF,#7AA6DCFF,#003C67FF,#8F7700FF,#3B3B3BFF,#A73030FF,#374E55FF,#DF8F44FF,#00A1D5FF,#B24745FF,#79AF97FF,#6A6599FF,#80796BFF,#BC3C29FF,#0072B5FF,#E18727FF,#20854EFF,#7876B1FF,#6F99ADFF,#FFDC91FF,#EE4C97FF"

tsne_cols <- trimws(unlist(strsplit(tsne_cols_text, ",")))
tsne_cols <- tsne_cols[tsne_cols != ""]

# 7.2 检查 tSNE 是否存在
if (!"tsne" %in% Reductions(srt2)) {
  stop("seurat 中没有 reduction tsne")
}

# 7.3 按细胞类型给色的 tSNE
p_tsne_cell <- DimPlot(
  srt2,
  reduction = "tsne",
  label = TRUE,
  cols = tsne_cols,
  label.size = tsne_label_size,
  pt.size = tsne_pt_size,
  raster = raster_tsne
) +
  theme_classic() +
  theme(
    panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"),
    legend.position = "right"
  )

pdf(
  file = file.path(out_dir, "分类TSNE（按细胞给色）.pdf"),
  width = tsne_width,
  height = tsne_height
)

print(p_tsne_cell)

dev.off()

# 7.4 按 group 给色的 tSNE
if ("group" %in% colnames(srt2@meta.data)) {
  
  p_tsne_group <- DimPlot(
    srt2,
    reduction = "tsne",
    label = TRUE,
    group.by = "group",
    cols = tsne_cols,
    label.size = tsne_label_size,
    pt.size = tsne_pt_size,
    raster = raster_tsne
  ) +
    theme_classic() +
    theme(
      panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"),
      legend.position = "right"
    )
  
  pdf(
    file = file.path(out_dir, "分类TSNE（按分组给色）.pdf"),
    width = tsne_width,
    height = tsne_height
  )
  
  print(p_tsne_group)
  
  dev.off()
}

# 7.5 按样本给色的 tSNE
if ("orig.ident" %in% colnames(srt2@meta.data)) {
  
  p_tsne_sample <- DimPlot(
    srt2,
    reduction = "tsne",
    label = TRUE,
    group.by = "orig.ident",
    cols = tsne_cols,
    label.size = tsne_label_size,
    pt.size = tsne_pt_size,
    raster = raster_tsne
  ) +
    theme_classic() +
    theme(
      panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"),
      legend.position = "right"
    )
  
  pdf(
    file = file.path(out_dir, "分类TSNE（按样本给色）.pdf"),
    width = tsne_width,
    height = tsne_height
  )
  
  print(p_tsne_sample)
  
  dev.off()
}