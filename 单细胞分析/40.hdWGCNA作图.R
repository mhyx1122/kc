# 1. 加载必要的 R 包

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(hdWGCNA)
})

options(future.globals.maxSize = 40 * 1024^3)


# 2. 检查 hdWGCNA 对象并创建输出文件夹

out_dir <- "hdWGCNA共表达网络分析"

if (!exists("seurat_hdWGCNA", envir = .GlobalEnv)) {
  stop("全局环境中没有 seurat_hdWGCNA 对象，请先完成软阈值筛选分析。")
}

srt <- get(
  "seurat_hdWGCNA",
  envir = .GlobalEnv
)

if (!inherits(srt, "Seurat")) {
  stop("全局环境中的 seurat_hdWGCNA 不是 Seurat 对象。")
}

if (!dir.exists(out_dir)) {
  dir.create(
    out_dir,
    recursive = TRUE
  )
}


# 3. 构建 hdWGCNA 共表达网络

tom_name <- "TCells"

# 设置为 NULL 时，使用前面软阈值筛选步骤确定的软阈值
# 也可以手动设置为数值，例如：
# soft_power <- 6
soft_power <- NULL

overwrite_tom <- TRUE

message("正在构建 hdWGCNA 共表达网络")

srt <- ConstructNetwork(
  srt,
  soft_power = soft_power,
  tom_name = tom_name,
  overwrite_tom = overwrite_tom
)


# 4. 提取并保存初始模块信息

modules <- GetModules(srt)

write.csv(
  modules,
  file = file.path(
    out_dir,
    "模块信息_完整表.csv"
  ),
  row.names = FALSE
)

modules_nogrey <- modules %>%
  subset(module != "grey")

write.csv(
  modules_nogrey,
  file = file.path(
    out_dir,
    "模块和基因对照表.csv"
  ),
  row.names = FALSE
)


# 5. 计算模块特征基因

# NULL 表示不按照批次变量进行协调
# 如需按照样本来源进行协调，可设置：
# group_by_vars <- "orig.ident"
group_by_vars <- NULL

message("正在计算模块特征基因")

srt <- ModuleEigengenes(
  srt,
  group.by.vars = group_by_vars
)


# 6. 提取并保存模块特征基因矩阵

hMEs <- GetMEs(srt)

MEs_raw <- GetMEs(
  srt,
  harmonized = FALSE
)

write.csv(
  hMEs,
  file = file.path(
    out_dir,
    "harmonized_module_eigengenes.csv"
  ),
  row.names = TRUE
)

write.csv(
  MEs_raw,
  file = file.path(
    out_dir,
    "raw_module_eigengenes.csv"
  ),
  row.names = TRUE
)


# 7. 计算模块连接度

connect_group_by <- "cellType"

connect_group_name <- "TCells"

message("正在计算模块连接度")

srt <- ModuleConnectivity(
  srt,
  group.by = connect_group_by,
  group_name = connect_group_name
)


# 8. 重命名模块

reset_module_name <- "TCells-M"

srt <- ResetModuleNames(
  srt,
  new_name = reset_module_name
)


# 9. 保存模块重命名后的模块和基因对照表

modules2 <- GetModules(srt)

modules2_nogrey <- modules2 %>%
  subset(module != "grey")

write.csv(
  modules2_nogrey,
  file = file.path(
    out_dir,
    "模块和基因对照表_重命名后.csv"
  ),
  row.names = FALSE
)


# 10. 提取并保存每个模块的 Hub 基因

n_hubs <- 10

hub_df <- GetHubGenes(
  srt,
  n_hubs = n_hubs
)

write.csv(
  hub_df,
  file = file.path(
    out_dir,
    "每个模块的核心基因.csv"
  ),
  row.names = FALSE
)


# 11. 计算模块表达评分

expr_score_genes <- 25

expr_score_method <- "UCell"

message("正在计算模块表达评分")

srt <- ModuleExprScore(
  srt,
  n_genes = expr_score_genes,
  method = expr_score_method
)


# 12. 根据 cellType 构建 cluster 字段

if ("cellType" %in% colnames(srt@meta.data)) {
  srt$cluster <- do.call(
    rbind,
    strsplit(
      as.character(srt$cellType),
      " "
    )
  )[, 1]
}


# 13. 将协调后的模块特征基因加入元数据

MEs2 <- GetMEs(
  srt,
  harmonized = TRUE
)

srt@meta.data <- cbind(
  srt@meta.data,
  MEs2
)


# 14. 运行模块 UMAP

run_module_umap_n_hubs <- 10

run_module_umap_n_neighbors <- 15

run_module_umap_min_dist <- 0.1

message("正在运行模块 UMAP")

srt <- RunModuleUMAP(
  srt,
  n_hubs = run_module_umap_n_hubs,
  n_neighbors = run_module_umap_n_neighbors,
  min_dist = run_module_umap_min_dist
)


# 15. 提取并保存模块 UMAP 坐标

umap_df <- GetModuleUMAP(srt)

write.csv(
  umap_df,
  file = file.path(
    out_dir,
    "模块UMAP坐标表.csv"
  ),
  row.names = FALSE
)


# 16. 输出各模块的网络图文件

message("正在输出各模块网络图")

ModuleNetworkPlot(
  srt,
  outdir = file.path(
    out_dir,
    "ModuleNetworks"
  )
)


# 17. 保存更新后的 hdWGCNA 对象

seurat_hdWGCNA <- srt


# 18. 绘制并保存模块聚类树状图

dendrogram_main <- "hdWGCNA Dendrogram"

dendrogram_width <- 10

dendrogram_height <- 8

dendrogram_file <- "1.Module_Dendrogram.pdf"

pdf(
  file = file.path(
    out_dir,
    dendrogram_file
  ),
  width = dendrogram_width,
  height = dendrogram_height
)

PlotDendrogram(
  srt,
  main = dendrogram_main
)

dev.off()


# 19. 绘制并保存 kME 排名图

kme_ncol <- 5

kme_width <- 14

kme_height <- 10

kme_file <- "2.Module_kME_Plot.pdf"

p_kme <- PlotKMEs(
  srt,
  ncol = kme_ncol
)

pdf(
  file = file.path(
    out_dir,
    kme_file
  ),
  width = kme_width,
  height = kme_height
)

print(p_kme)

dev.off()


# 20. 绘制并保存模块特征图

hmes_ncol <- 6

hmes_width <- 18

hmes_height <- 10

hmes_file <- "3.ModuleFeaturePlot_hMEs.pdf"

hmes_plot_list <- ModuleFeaturePlot(
  srt,
  features = "hMEs",
  order = TRUE
)

p_hmes <- wrap_plots(
  hmes_plot_list,
  ncol = hmes_ncol
)

pdf(
  file = file.path(
    out_dir,
    hmes_file
  ),
  width = hmes_width,
  height = hmes_height
)

print(p_hmes)

dev.off()


# 21. 绘制并保存模块评分图

# 可选：
# scores_order <- "shuffle"
# scores_order <- TRUE
scores_order <- "shuffle"

scores_ncol <- 6

scores_width <- 18

scores_height <- 10

scores_file <- "4.ModuleFeaturePlot_scores.pdf"

scores_plot_list <- ModuleFeaturePlot(
  srt,
  features = "scores",
  order = scores_order,
  ucell = TRUE
)

p_scores <- wrap_plots(
  scores_plot_list,
  ncol = scores_ncol
)

pdf(
  file = file.path(
    out_dir,
    scores_file
  ),
  width = scores_width,
  height = scores_height
)

print(p_scores)

dev.off()


# 22. 筛选用于模块雷达图的细胞

radar_group_by <- "cluster"

radar_celltype_filter <- "TCells"

radar_axis_label_size <- 4

radar_grid_label_size <- 4

radar_width <- 10

radar_height <- 8

radar_file <- "5.ModuleRadarPlot.pdf"

if (!"cellType" %in% colnames(srt@meta.data)) {
  stop("seurat_hdWGCNA@meta.data 中不存在 cellType 列。")
}

if (!radar_group_by %in% colnames(srt@meta.data)) {
  stop(
    paste0(
      "seurat_hdWGCNA@meta.data 中不存在雷达图分组列：",
      radar_group_by
    )
  )
}

barcode_use <- srt@meta.data %>%
  subset(cellType == radar_celltype_filter) %>%
  rownames()

if (length(barcode_use) == 0) {
  stop(
    paste0(
      "没有筛选到 cellType 为 ",
      radar_celltype_filter,
      " 的细胞。"
    )
  )
}


# 23. 绘制并保存模块雷达图

p_radar <- ModuleRadarPlot(
  srt,
  group.by = radar_group_by,
  barcodes = barcode_use,
  axis.label.size = radar_axis_label_size,
  grid.label.size = radar_grid_label_size
)

pdf(
  file = file.path(
    out_dir,
    radar_file
  ),
  width = radar_width,
  height = radar_height
)

print(p_radar)

dev.off()


# 24. 绘制并保存模块相关图

correlogram_width <- 10

correlogram_height <- 8

correlogram_file <- "6.ModuleCorrelogram.pdf"

pdf(
  file = file.path(
    out_dir,
    correlogram_file
  ),
  width = correlogram_width,
  height = correlogram_height
)

ModuleCorrelogram(srt)

dev.off()


# 25. 绘制并保存全部模块的 Hub 基因网络图

hub_all_n_hubs <- 5

hub_all_n_other <- 3

hub_all_edge_prop <- 0.5

hub_all_edge_alpha <- 1

hub_all_width <- 16

hub_all_height <- 12

hub_all_file <- "7.HubGeneNetworkPlot_all_modules.pdf"

pdf(
  file = file.path(
    out_dir,
    hub_all_file
  ),
  width = hub_all_width,
  height = hub_all_height
)

par(
  mar = c(3, 3, 2, 1),
  family = "sans"
)

HubGeneNetworkPlot(
  srt,
  n_hubs = hub_all_n_hubs,
  n_other = hub_all_n_other,
  edge_prop = hub_all_edge_prop,
  edge.alpha = hub_all_edge_alpha,
  mods = "all"
)

dev.off()


# 26. 提取用于前几个模块网络图的模块名称

hub_top5_n_hubs <- 10

hub_top5_n_other <- 3

hub_top5_edge_prop <- 0.75

hub_top5_n_mods <- 5

hub_top5_width <- 16

hub_top5_height <- 12

hub_top5_file <- "8.HubGeneNetworkPlot_top5_modules.pdf"

modules_plot <- GetModules(srt)

mods <- levels(
  modules_plot$module
)

mods <- mods[
  mods != "grey"
]

if (length(mods) == 0) {
  stop("没有可用于绘制 Hub 基因网络图的非 grey 模块。")
}

mods_use <- mods[
  1:min(
    hub_top5_n_mods,
    length(mods)
  )
]


# 27. 绘制并保存前几个模块的 Hub 基因网络图

pdf(
  file = file.path(
    out_dir,
    hub_top5_file
  ),
  width = hub_top5_width,
  height = hub_top5_height
)

par(
  mar = c(3, 3, 2, 1),
  family = "sans"
)

HubGeneNetworkPlot(
  srt,
  n_hubs = hub_top5_n_hubs,
  n_other = hub_top5_n_other,
  edge_prop = hub_top5_edge_prop,
  mods = mods_use
)

dev.off()


# 28. 设置模块 UMAP 散点图颜色和点大小

palette_module_umap_text <- paste0(
  "#E64B35FF,#4DBBD5FF,#00A087FF,#3C5488FF,",
  "#F39B7FFF,#8491B4FF,#91D1C2FF,#7E6148FF,",
  "#0073C2FF,#EFC000FF,#868686FF,#CD534CFF,",
  "#7AA6DCFF,#003C67FF,#8F7700FF,#3B3B3BFF,",
  "#A73030FF,#374E55FF,#DF8F44FF,#00A1D5FF,",
  "#B24745FF,#79AF97FF,#6A6599FF,#80796BFF,",
  "#BC3C29FF,#0072B5FF,#E18727FF,#20854EFF,",
  "#7876B1FF,#6F99ADFF,#FFDC91FF,#EE4C97FF"
)

module_umap_pt_scale <- 2

module_umap_scatter_width <- 12

module_umap_scatter_height <- 10

module_umap_scatter_file <- "9.ModuleUMAP_scatter.pdf"

palette_module_umap <- unlist(
  strsplit(
    palette_module_umap_text,
    ","
  )
)

palette_module_umap <- trimws(
  palette_module_umap
)

palette_module_umap <- palette_module_umap[
  palette_module_umap != ""
]


# 29. 整理模块 UMAP 散点图数据

module_umap_plot_df <- umap_df

if (nrow(module_umap_plot_df) == 0) {
  stop("模块 UMAP 坐标表中没有数据。")
}

if ("module" %in% colnames(module_umap_plot_df)) {
  
  module_names <- unique(
    as.character(
      module_umap_plot_df$module
    )
  )
  
  module_color_map <- setNames(
    rep(
      palette_module_umap,
      length.out = length(module_names)
    ),
    module_names
  )
  
  module_umap_plot_df$plot_color <- module_color_map[
    as.character(
      module_umap_plot_df$module
    )
  ]
  
} else if ("color" %in% colnames(module_umap_plot_df)) {
  
  module_umap_plot_df$plot_color <-
    module_umap_plot_df$color
  
} else {
  
  module_umap_plot_df$plot_color <-
    palette_module_umap[1]
}


# 30. 绘制并保存模块 UMAP 散点图

p_module_umap_scatter <- ggplot(
  module_umap_plot_df,
  aes(
    x = UMAP1,
    y = UMAP2
  )
) +
  geom_point(
    color = module_umap_plot_df$plot_color,
    size = module_umap_plot_df$kME *
      module_umap_pt_scale
  ) +
  theme_classic() +
  labs(
    x = "UMAP1",
    y = "UMAP2"
  )

ggsave(
  filename = file.path(
    out_dir,
    module_umap_scatter_file
  ),
  plot = p_module_umap_scatter,
  width = module_umap_scatter_width,
  height = module_umap_scatter_height,
  device = "pdf"
)


# 31. 绘制并保存模块 UMAP 网络图

module_umap_edge_alpha <- 0.25

module_umap_sample_edges <- TRUE

module_umap_edge_prop <- 0.1

module_umap_label_hubs <- 2

module_umap_keep_grey_edges <- FALSE

module_umap_network_width <- 16

module_umap_network_height <- 12

module_umap_network_file <- "10.ModuleUMAPPlot.pdf"

pdf(
  file = file.path(
    out_dir,
    module_umap_network_file
  ),
  width = module_umap_network_width,
  height = module_umap_network_height
)

ModuleUMAPPlot(
  srt,
  edge.alpha = module_umap_edge_alpha,
  sample_edges = module_umap_sample_edges,
  edge_prop = module_umap_edge_prop,
  label_hubs = module_umap_label_hubs,
  keep_grey_edges = module_umap_keep_grey_edges
)

dev.off()


# 32. 保存本次 hdWGCNA 分析参数

params_file <- "hdWGCNA_parameters.txt"

param_text <- paste0(
  "运行流程说明：\n",
  "1. 使用 hdWGCNA::ConstructNetwork 构建基因共表达网络，并生成 TOM 相关结果。\n",
  "2. 使用 hdWGCNA::GetModules 提取模块信息，并保存完整模块信息表与非 grey 模块表。\n",
  "3. 使用 hdWGCNA::ModuleEigengenes 计算模块特征基因；",
  "若 group.by.vars 设置为 orig.ident，则按照 orig.ident 进行协调；",
  "设置为 NULL 时不进行批次协调。\n",
  "4. 使用 hdWGCNA::GetMEs 获取协调后的 MEs 和原始 MEs，并保存为 CSV 文件。\n",
  "5. 使用 hdWGCNA::ModuleConnectivity 计算模块连接度。\n",
  "6. 使用 hdWGCNA::ResetModuleNames 重命名模块。\n",
  "7. 使用 hdWGCNA::GetHubGenes 提取 Hub 基因，并输出结果表。\n",
  "8. 使用 hdWGCNA::ModuleExprScore 计算模块评分。\n",
  "9. 若 meta.data 中存在 cellType，则基于 cellType 构建 cluster 字段。\n",
  "10. 使用 hdWGCNA::RunModuleUMAP 进行模块 UMAP 降维，",
  "并使用 hdWGCNA::GetModuleUMAP 提取坐标。\n",
  "11. 使用 hdWGCNA::ModuleNetworkPlot 输出各模块的网络图文件。\n",
  "12. 可视化结果包括树状图、kME排名图、模块特征图、模块评分图、",
  "模块雷达图、模块相关图、Hub基因网络图、模块UMAP散点图和模块UMAP网络图。\n",
  "\n",
  "本次运行参数：\n",
  "- tom_name：", tom_name, "\n",
  "- soft_power：",
  ifelse(
    is.null(soft_power),
    "自动设置（NULL）",
    as.character(soft_power)
  ),
  "\n",
  "- overwrite_tom：", overwrite_tom, "\n",
  "- group.by.vars：",
  ifelse(
    is.null(group_by_vars),
    "NULL（不去批次）",
    group_by_vars
  ),
  "\n",
  "- ModuleConnectivity group.by：",
  connect_group_by,
  "\n",
  "- ModuleConnectivity group_name：",
  connect_group_name,
  "\n",
  "- ResetModuleNames new_name：",
  reset_module_name,
  "\n",
  "- GetHubGenes n_hubs：",
  n_hubs,
  "\n",
  "- ModuleExprScore n_genes：",
  expr_score_genes,
  "\n",
  "- ModuleExprScore method：",
  expr_score_method,
  "\n",
  "- RunModuleUMAP n_hubs：",
  run_module_umap_n_hubs,
  "\n",
  "- RunModuleUMAP n_neighbors：",
  run_module_umap_n_neighbors,
  "\n",
  "- RunModuleUMAP min_dist：",
  run_module_umap_min_dist,
  "\n",
  "\n",
  "结果输出说明：\n",
  "- 模块信息_完整表.csv：完整模块注释结果。\n",
  "- 模块和基因对照表.csv：去除 grey 模块后的模块-基因对照表。\n",
  "- 模块和基因对照表_重命名后.csv：模块重命名后的模块-基因对照表。\n",
  "- 每个模块的核心基因.csv：Hub 基因结果表。\n",
  "- harmonized_module_eigengenes.csv：协调后的模块特征基因。\n",
  "- raw_module_eigengenes.csv：原始模块特征基因。\n",
  "- 模块UMAP坐标表.csv：模块 UMAP 坐标结果。\n",
  "- ModuleNetworks/：各模块网络图输出目录。\n"
)

writeLines(
  text = param_text,
  con = file.path(
    out_dir,
    params_file
  )
)

message(
  "hdWGCNA 共表达网络分析完成，结果已保存至：",
  out_dir
)