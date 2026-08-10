suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

# 1. 检查全局环境中是否存在 seurat 对象
if (!exists("seurat", envir = .GlobalEnv)) {
  stop("全局环境中没有 seurat 对象")
}

srt <- get("seurat", envir = .GlobalEnv)

if (!exists("do_DotPlot", envir = .GlobalEnv)) {
  stop("全局环境中没有 do_DotPlot 函数")
}

do_fun <- get("do_DotPlot", envir = .GlobalEnv)


# 2. 输出目录模块

# 2.1 输出目录参数
out_dir <- "7.基因表达量可视化"

# 2.2 创建输出目录
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}


# 3. 基因输入模块

# 3.1 输入基因参数
genes_text <- "ICOS, GZMK, CD6, BCL11B, TRAT1"

genes <- unlist(strsplit(genes_text, ","))
genes <- trimws(genes)
genes <- genes[genes != ""]

if (length(genes) == 0) {
  stop("请输入至少 1 个基因")
}


# 4. 批量保存参数模块

# 4.1 批量保存宽高
img_width <- 8
img_height <- 6


# 5. DotPlot 模块

# 5.1 DotPlot 参数
dot_scale <- 8
dot_font_size <- 12
dot_legend_length <- 10
dot_legend_width <- 1
dot_seq_palette <- "Blues"
dot_seq_palette_group <- "Blues"

# 5.2 DotPlot 当前保存参数
w_dot <- 10
h_dot <- 8
name_dot <- "DotPlot_selected_genes"

# 5.3 生成 cellType DotPlot
p_dot_celltype <- do_fun(
  sample = srt,
  features = genes,
  group.by = "cellType",
  dot.scale = dot_scale,
  legend.framewidth = 1,
  font.size = dot_font_size,
  legend.length = dot_legend_length,
  legend.width = dot_legend_width,
  sequential.palette = dot_seq_palette,
  sequential.direction = 1
)

# 5.4 生成 group DotPlot
p_dot_group <- do_fun(
  sample = srt,
  features = genes,
  group.by = "group",
  dot.scale = dot_scale,
  legend.framewidth = 1,
  font.size = dot_font_size,
  legend.length = dot_legend_length,
  legend.width = dot_legend_width,
  sequential.palette = dot_seq_palette_group,
  sequential.direction = 1
)

# 5.5 保存 DotPlot
ggsave(
  filename = file.path(out_dir, paste0(name_dot, "_by_cellType.pdf")),
  plot = p_dot_celltype,
  width = w_dot,
  height = h_dot,
  device = "pdf"
)

ggsave(
  filename = file.path(out_dir, paste0(name_dot, "_by_group.pdf")),
  plot = p_dot_group,
  width = w_dot,
  height = h_dot,
  device = "pdf"
)


# 6. 单基因 FeaturePlot 模块

# 6.1 FeaturePlot 参数
selected_gene_feature <- genes[1]
feature_plot_pt_size <- 0.5

feature_plot_cols_celltype_text <- "lightgrey,#c83732"
feature_plot_cols_celltype <- unlist(strsplit(feature_plot_cols_celltype_text, ","))
feature_plot_cols_celltype <- trimws(feature_plot_cols_celltype)
feature_plot_cols_celltype <- feature_plot_cols_celltype[feature_plot_cols_celltype != ""]

feature_plot_cols_group_text <- "lightgrey,#c83732"
feature_plot_cols_group <- unlist(strsplit(feature_plot_cols_group_text, ","))
feature_plot_cols_group <- trimws(feature_plot_cols_group)
feature_plot_cols_group <- feature_plot_cols_group[feature_plot_cols_group != ""]

if (length(feature_plot_cols_celltype) < 2) {
  stop("FeaturePlot by_cellType 颜色至少需要 2 个")
}

if (length(feature_plot_cols_group) < 2) {
  stop("FeaturePlot by_group 颜色至少需要 2 个")
}

# 6.2 FeaturePlot 当前保存参数
w_feature <- 8
h_feature <- 6
name_feature <- "FeaturePlot"

# 6.3 生成 by_cellType FeaturePlot
p_feature_celltype <- FeaturePlot(
  srt,
  features = selected_gene_feature,
  pt.size = feature_plot_pt_size,
  raster = FALSE,
  cols = feature_plot_cols_celltype
) +
  ggtitle(paste("FeaturePlot of", selected_gene_feature, "by cellType"))

# 6.4 生成 by_group FeaturePlot
p_feature_group <- FeaturePlot(
  srt,
  features = selected_gene_feature,
  pt.size = feature_plot_pt_size,
  raster = FALSE,
  cols = feature_plot_cols_group,
  split.by = "group"
) +
  ggtitle(paste("FeaturePlot of", selected_gene_feature, "by group"))

# 6.5 保存单基因 FeaturePlot
ggsave(
  filename = file.path(out_dir, paste0(name_feature, "_", selected_gene_feature, "_by_cellType.pdf")),
  plot = p_feature_celltype,
  width = w_feature,
  height = h_feature,
  device = "pdf"
)

ggsave(
  filename = file.path(out_dir, paste0(name_feature, "_", selected_gene_feature, "_by_group.pdf")),
  plot = p_feature_group,
  width = w_feature,
  height = h_feature,
  device = "pdf"
)


# 7. 单基因 VlnPlot 模块

# 7.1 VlnPlot 参数
selected_gene_vln <- genes[1]
vln_plot_pt_size <- 0.1

vln_fill_celltype_text <- "#E64B35FF,#4DBBD5FF,#00A087FF,#3C5488FF,#F39B7FFF,#8491B4FF,#91D1C2FF,#7E6148FF,#0073C2FF,#EFC000FF,#868686FF,#CD534CFF,#7AA6DCFF,#003C67FF,#8F7700FF,#3B3B3BFF,#A73030FF,#374E55FF,#DF8F44FF,#00A1D5FF,#B24745FF,#79AF97FF,#6A6599FF,#80796BFF,#BC3C29FF,#0072B5FF,#E18727FF,#20854EFF,#7876B1FF,#6F99ADFF,#FFDC91FF,#EE4C97FF"

vln_fill_celltype <- unlist(strsplit(vln_fill_celltype_text, ","))
vln_fill_celltype <- trimws(vln_fill_celltype)
vln_fill_celltype <- vln_fill_celltype[vln_fill_celltype != ""]

vln_fill_group_text <- "#E64B35FF,#4DBBD5FF,#00A087FF,#3C5488FF,#F39B7FFF,#8491B4FF,#91D1C2FF,#7E6148FF,#0073C2FF,#EFC000FF,#868686FF,#CD534CFF,#7AA6DCFF,#003C67FF,#8F7700FF,#3B3B3BFF,#A73030FF,#374E55FF,#DF8F44FF,#00A1D5FF,#B24745FF,#79AF97FF,#6A6599FF,#80796BFF,#BC3C29FF,#0072B5FF,#E18727FF,#20854EFF,#7876B1FF,#6F99ADFF,#FFDC91FF,#EE4C97FF"

vln_fill_group <- unlist(strsplit(vln_fill_group_text, ","))
vln_fill_group <- trimws(vln_fill_group)
vln_fill_group <- vln_fill_group[vln_fill_group != ""]

# 7.2 VlnPlot 当前保存参数
w_vln <- 8
h_vln <- 6
name_vln <- "VlnPlot"

# 7.3 生成 by_cellType VlnPlot
p_vln_celltype <- VlnPlot(
  srt,
  features = selected_gene_vln,
  pt.size = vln_plot_pt_size,
  group.by = "cellType"
) +
  ggtitle(paste("VlnPlot of", selected_gene_vln, "by cellType"))

if (length(vln_fill_celltype) > 0) {
  p_vln_celltype <- p_vln_celltype +
    scale_fill_manual(values = vln_fill_celltype)
}

# 7.4 生成 by_group VlnPlot
p_vln_group <- VlnPlot(
  srt,
  features = selected_gene_vln,
  pt.size = vln_plot_pt_size,
  group.by = "group"
) +
  ggtitle(paste("VlnPlot of", selected_gene_vln, "by group"))

if (length(vln_fill_group) > 0) {
  p_vln_group <- p_vln_group +
    scale_fill_manual(values = vln_fill_group)
}

# 7.5 保存单基因 VlnPlot
ggsave(
  filename = file.path(out_dir, paste0(name_vln, "_", selected_gene_vln, "_by_cellType.pdf")),
  plot = p_vln_celltype,
  width = w_vln,
  height = h_vln,
  device = "pdf"
)

ggsave(
  filename = file.path(out_dir, paste0(name_vln, "_", selected_gene_vln, "_by_group.pdf")),
  plot = p_vln_group,
  width = w_vln,
  height = h_vln,
  device = "pdf"
)


# 8. 批量保存 FeaturePlot 和 VlnPlot 模块

# 8.1 批量保存 DotPlot
ggsave(
  filename = file.path(out_dir, "DotPlot_selected_genes_by_cellType.pdf"),
  plot = p_dot_celltype,
  width = img_width,
  height = img_height,
  device = "pdf"
)

ggsave(
  filename = file.path(out_dir, "DotPlot_selected_genes_by_group.pdf"),
  plot = p_dot_group,
  width = img_width,
  height = img_height,
  device = "pdf"
)

# 8.2 批量保存每个基因的 FeaturePlot 和 VlnPlot
for (gene in genes) {
  
  p_feature_celltype_loop <- FeaturePlot(
    srt,
    features = gene,
    pt.size = feature_plot_pt_size,
    raster = FALSE,
    cols = feature_plot_cols_celltype
  ) +
    ggtitle(paste("FeaturePlot of", gene, "by cellType"))
  
  ggsave(
    filename = file.path(out_dir, paste0("FeaturePlot_", gene, "_by_cellType.pdf")),
    plot = p_feature_celltype_loop,
    width = img_width,
    height = img_height,
    device = "pdf"
  )
  
  p_feature_group_loop <- FeaturePlot(
    srt,
    features = gene,
    pt.size = feature_plot_pt_size,
    raster = FALSE,
    cols = feature_plot_cols_group,
    split.by = "group"
  ) +
    ggtitle(paste("FeaturePlot of", gene, "by group"))
  
  ggsave(
    filename = file.path(out_dir, paste0("FeaturePlot_", gene, "_by_group.pdf")),
    plot = p_feature_group_loop,
    width = img_width,
    height = img_height,
    device = "pdf"
  )
  
  p_vln_celltype_loop <- VlnPlot(
    srt,
    features = gene,
    pt.size = vln_plot_pt_size,
    group.by = "cellType"
  ) +
    ggtitle(paste("VlnPlot of", gene, "by cellType"))
  
  if (length(vln_fill_celltype) > 0) {
    p_vln_celltype_loop <- p_vln_celltype_loop +
      scale_fill_manual(values = vln_fill_celltype)
  }
  
  ggsave(
    filename = file.path(out_dir, paste0("VlnPlot_", gene, "_by_cellType.pdf")),
    plot = p_vln_celltype_loop,
    width = img_width,
    height = img_height,
    device = "pdf"
  )
  
  p_vln_group_loop <- VlnPlot(
    srt,
    features = gene,
    pt.size = vln_plot_pt_size,
    group.by = "group"
  ) +
    ggtitle(paste("VlnPlot of", gene, "by group"))
  
  if (length(vln_fill_group) > 0) {
    p_vln_group_loop <- p_vln_group_loop +
      scale_fill_manual(values = vln_fill_group)
  }
  
  ggsave(
    filename = file.path(out_dir, paste0("VlnPlot_", gene, "_by_group.pdf")),
    plot = p_vln_group_loop,
    width = img_width,
    height = img_height,
    device = "pdf"
  )
}


# 9. 参数记录模块

# 9.1 参数文件保存参数
name_params <- "gene_plots_parameters"

# 9.2 生成参数记录文本
param_text <- paste0(
  "本次分析参数总结：\n",
  "- genes：", paste(genes, collapse = ", "), "\n",
  "- img_width：", img_width, "\n",
  "- img_height：", img_height, "\n",
  "- DotPlot dot.scale：", dot_scale, "\n",
  "- DotPlot font.size：", dot_font_size, "\n",
  "- DotPlot legend.length：", dot_legend_length, "\n",
  "- DotPlot legend.width：", dot_legend_width, "\n",
  "- DotPlot cellType sequential.palette：", dot_seq_palette, "\n",
  "- DotPlot group sequential.palette：", dot_seq_palette_group, "\n",
  "- FeaturePlot pt.size：", feature_plot_pt_size, "\n",
  "- FeaturePlot cellType cols：", feature_plot_cols_celltype_text, "\n",
  "- FeaturePlot group cols：", feature_plot_cols_group_text, "\n",
  "- VlnPlot pt.size：", vln_plot_pt_size, "\n",
  "- VlnPlot cellType fill：", vln_fill_celltype_text, "\n",
  "- VlnPlot group fill：", vln_fill_group_text, "\n"
)

# 9.3 保存参数记录
writeLines(
  param_text,
  con = file.path(out_dir, paste0(name_params, ".txt"))
)