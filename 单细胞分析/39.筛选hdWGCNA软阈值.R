# 1. 加载必要的 R 包

suppressPackageStartupMessages({
  library(Seurat)
  library(hdWGCNA)
  library(ggplot2)
  library(patchwork)
})


# 2. 检查 Seurat 对象并创建输出文件夹

out_dir <- "hdWGCNA共表达网络分析"

if (!exists("seurat", envir = .GlobalEnv)) {
  stop("全局环境中没有找到 seurat 对象。")
}

srt <- get("seurat", envir = .GlobalEnv)

if (!inherits(srt, "Seurat")) {
  stop("全局环境中的 seurat 不是 Seurat 对象。")
}

if (!dir.exists(out_dir)) {
  dir.create(
    out_dir,
    recursive = TRUE
  )
}


# 3. 初始化 hdWGCNA 分析

gene_select <- "fraction"

fraction <- 0.05

wgcna_name <- "tutorial"

message("正在初始化 hdWGCNA 分析")

srt <- SetupForWGCNA(
  srt,
  gene_select = gene_select,
  fraction = fraction,
  wgcna_name = wgcna_name
)


# 4. 按分组信息构建 metacells

group_by <- "cellType,orig.ident"

k <- 25

max_shared <- 10

min_cells <- 100

ident_group <- "cellType"

reduction <- "umap"

group_by_vec <- unlist(
  strsplit(
    group_by,
    split = ","
  )
)

group_by_vec <- trimws(
  group_by_vec
)

message("正在构建 metacells，请耐心等待运行完成")

srt <- MetacellsByGroups(
  seurat_obj = srt,
  group.by = group_by_vec,
  k = k,
  max_shared = max_shared,
  min_cells = min_cells,
  ident.group = ident_group,
  reduction = reduction
)


# 5. 对 metacells 表达矩阵进行标准化

message("正在标准化 metacells 表达矩阵")

srt <- NormalizeMetacells(
  srt
)


# 6. 设置用于 WGCNA 分析的表达矩阵

group_name <- "TCells"

group_by_datexpr <- "cellType"

assay <- "RNA"

layer <- "data"

message(
  "正在提取 ",
  group_name,
  " 的表达矩阵"
)

srt <- SetDatExpr(
  srt,
  group_name = group_name,
  group.by = group_by_datexpr,
  assay = assay,
  layer = layer
)


# 7. 测试不同的软阈值

networkType <- "signed"

message("正在测试不同的软阈值")

srt <- TestSoftPowers(
  srt,
  networkType = networkType
)


# 8. 获取并保存软阈值评估图

softpower_width <- 12

softpower_height <- 6

softpower_file <- "1.软阈值评估图.pdf"

plot_list <- PlotSoftPowers(
  srt
)

p_softpower <- wrap_plots(
  plot_list,
  ncol = 2
)

ggsave(
  filename = file.path(
    out_dir,
    softpower_file
  ),
  plot = p_softpower,
  width = softpower_width,
  height = softpower_height,
  device = "pdf"
)

message(
  "软阈值评估图已保存至：",
  file.path(out_dir, softpower_file)
)


# 9. 获取并保存软阈值评估表

power_table_file <- "软阈值评估表.csv"

power_table <- GetPowerTable(
  srt
)

write.csv(
  power_table,
  file = file.path(
    out_dir,
    power_table_file
  ),
  row.names = FALSE
)

message(
  "软阈值评估表已保存至：",
  file.path(out_dir, power_table_file)
)


# 10. 保存 hdWGCNA 分析结果对象

seurat_hdWGCNA <- srt


# 11. 保存本次分析参数

params_file <- "hdWGCNA_parameters.txt"

param_text <- paste0(
  "hdWGCNA 软阈值筛选分析流程说明：\n\n",
  
  "1. SetupForWGCNA：使用 hdWGCNA::SetupForWGCNA() 函数初始化 WGCNA 分析，",
  "gene_select 参数设置为 ", gene_select, "，",
  "fraction 参数设置为 ", fraction, "，",
  "要求基因至少在 ", fraction * 100, "% 的细胞中表达。\n\n",
  
  "2. MetacellsByGroups：使用 hdWGCNA::MetacellsByGroups() 函数构建 metacells，",
  "按照 ", paste(group_by_vec, collapse = "、"), " 分组，",
  "每个 metacell 包含 ", k, " 个细胞，",
  "max_shared 参数为 ", max_shared, "，",
  "min_cells 参数为 ", min_cells, "，",
  "ident.group 设置为 ", ident_group, "，",
  "并基于 ", reduction, " 降维结果进行聚合。\n\n",
  
  "3. NormalizeMetacells：使用 hdWGCNA::NormalizeMetacells() 函数",
  "对 metacells 表达矩阵进行标准化处理。\n\n",
  
  "4. SetDatExpr：使用 hdWGCNA::SetDatExpr() 函数设置用于 WGCNA 的表达矩阵，",
  "group_name 设置为 ", group_name, "，",
  "group.by 设置为 ", group_by_datexpr, "，",
  "assay 设置为 ", assay, "，",
  "layer 设置为 ", layer, "。\n\n",
  
  "5. TestSoftPowers：使用 hdWGCNA::TestSoftPowers() 函数测试不同软阈值，",
  "networkType 设置为 ", networkType, "，",
  "用于评估网络拓扑特性。\n\n",
  
  "6. PlotSoftPowers：使用 hdWGCNA::PlotSoftPowers() 函数",
  "可视化软阈值评估结果，包括 Scale Independence 和 Mean Connectivity 两个指标。\n\n",
  
  "本次分析参数总结：\n",
  "- gene_select：", gene_select, "\n",
  "- fraction：", fraction, "\n",
  "- wgcna_name：", wgcna_name, "\n",
  "- group.by：", paste(group_by_vec, collapse = ", "), "\n",
  "- k：", k, "\n",
  "- max_shared：", max_shared, "\n",
  "- min_cells：", min_cells, "\n",
  "- ident.group：", ident_group, "\n",
  "- reduction：", reduction, "\n",
  "- group_name：", group_name, "\n",
  "- SetDatExpr group.by：", group_by_datexpr, "\n",
  "- assay：", assay, "\n",
  "- layer：", layer, "\n",
  "- networkType：", networkType, "\n",
  "- 输出文件夹：", out_dir, "\n",
  "- 结果对象：seurat_hdWGCNA\n\n",
  
  "写作提示：\n",
  "1. 通过 hdWGCNA 方法对单细胞数据进行加权基因共表达网络分析。\n",
  "2. 首先构建 metacells 以降低单细胞表达矩阵的稀疏性，",
  "随后进行软阈值筛选以确定网络构建参数。\n",
  "3. 软阈值选择需要综合考虑 Scale-free topology fit 和 Mean connectivity，",
  "通常选择拟合指数较高且平均连接度适中的最小阈值。\n",
  "4. 后续可基于选定的软阈值构建共表达网络并识别功能模块。\n"
)

writeLines(
  text = param_text,
  con = file.path(
    out_dir,
    params_file
  )
)

message(
  "hdWGCNA 软阈值筛选分析完成，结果已保存至：",
  out_dir
)