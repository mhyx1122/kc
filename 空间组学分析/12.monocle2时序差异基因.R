suppressPackageStartupMessages({
  library(Seurat)
  library(monocle)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

# 1. 检查并读取拟时序分析对象

# 1.1 参数设置

out_dir <- "Monocle2结果"

# 1.2 创建输出文件夹

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# 1.3 检查全局对象

if (!exists("mycds", envir = .GlobalEnv)) {
  stop("全局环境中没有 mycds 对象，请先完成上一阶段Monocle2拟时序分析。")
}

if (!exists("disp.genes", envir = .GlobalEnv)) {
  stop("全局环境中没有 disp.genes 对象，请先完成上一阶段Monocle2拟时序分析。")
}

if (!exists("Spatial_Data", envir = .GlobalEnv)) {
  stop("全局环境中没有 Spatial_Data 对象，请先准备空间转录组对象。")
}

# 1.4 读取对象

mycds_obj <- get("mycds", envir = .GlobalEnv)
disp_genes <- get("disp.genes", envir = .GlobalEnv)
seurat_obj <- get("Spatial_Data", envir = .GlobalEnv)

if (!inherits(mycds_obj, "CellDataSet")) {
  stop("mycds 不是有效的Monocle2 CellDataSet对象。")
}

if (!inherits(seurat_obj, "Seurat")) {
  stop("Spatial_Data 不是有效的Seurat对象。")
}

# 1.5 检查轨迹构建基因

disp_genes <- as.character(disp_genes)
disp_genes <- disp_genes[!is.na(disp_genes)]
disp_genes <- disp_genes[disp_genes != ""]
disp_genes <- unique(disp_genes)

missing_disp_genes <- setdiff(disp_genes, rownames(mycds_obj))

if (length(missing_disp_genes) > 0) {
  warning(
    paste0(
      "以下disp.genes不在mycds中，将不参与差异分析：",
      paste(missing_disp_genes, collapse = ", ")
    )
  )
}

disp_genes <- intersect(disp_genes, rownames(mycds_obj))

if (length(disp_genes) == 0) {
  stop("disp.genes中没有可用于拟时序差异分析的基因。")
}

# 2. 进行拟时序差异基因分析

# 2.1 参数设置

cores_used <- 1
num_top_genes <- 8
name_diff_csv <- "Pseudotime_Differential_Genes"

# 2.2 运行拟时序差异基因分析

Time_diff <- differentialGeneTest(
  mycds_obj[disp_genes, ],
  cores = cores_used,
  fullModelFormulaStr = "~sm.ns(Pseudotime)"
)

# 2.3 检查差异分析结果

required_columns <- c("gene_short_name", "qval")
missing_columns <- setdiff(required_columns, colnames(Time_diff))

if (length(missing_columns) > 0) {
  stop(
    paste0(
      "differentialGeneTest结果中缺少以下列：",
      paste(missing_columns, collapse = ", ")
    )
  )
}

if (nrow(Time_diff) == 0) {
  stop("differentialGeneTest没有返回基因结果。")
}

# 2.4 提取全部检测基因

Time_genes <- Time_diff %>%
  pull(gene_short_name) %>%
  as.character()

Time_genes <- Time_genes[!is.na(Time_genes)]
Time_genes <- Time_genes[Time_genes != ""]
Time_genes <- unique(Time_genes)
Time_genes <- intersect(Time_genes, rownames(mycds_obj))

if (length(Time_genes) == 0) {
  stop("没有获得可用于拟时序热图的基因。")
}

# 2.5 根据q值选取Top基因

top_genes <- Time_diff %>%
  slice_min(
    order_by = qval,
    n = num_top_genes,
    with_ties = FALSE
  ) %>%
  pull(gene_short_name) %>%
  as.character()

top_genes <- top_genes[!is.na(top_genes)]
top_genes <- top_genes[top_genes != ""]
top_genes <- unique(top_genes)
top_genes <- intersect(top_genes, rownames(mycds_obj))

if (length(top_genes) == 0) {
  stop("没有获得可用于作图的Top拟时序差异基因。")
}

# 2.6 保存拟时序差异基因结果

diff_csv_file <- file.path(
  out_dir,
  paste0(name_diff_csv, ".csv")
)

write.csv(
  Time_diff,
  file = diff_csv_file,
  row.names = TRUE
)

# 3. 绘制所有检测基因的拟时序热图

# 3.1 参数设置

save_heatmap <- TRUE
num_clusters <- 4

w_all_heatmap <- 8
h_all_heatmap <- 8
name_all_heatmap <- "2.所有基因的热图"

# 3.2 检查热图聚类数

if (save_heatmap && num_clusters > length(Time_genes)) {
  stop(
    paste0(
      "num_clusters不能大于所有检测基因数量。",
      "当前num_clusters为 ",
      num_clusters,
      "，基因数量为 ",
      length(Time_genes),
      "。"
    )
  )
}

# 3.3 绘制所有检测基因热图

if (save_heatmap) {
  plot_pseudotime_heatmap(
    mycds_obj[Time_genes, ],
    num_clusters = num_clusters,
    show_rownames = FALSE,
    return_heatmap = FALSE
  )
  
  # 3.4 保存所有检测基因热图
  
  all_heatmap_file <- file.path(
    out_dir,
    paste0(name_all_heatmap, ".pdf")
  )
  
  pdf(
    file = all_heatmap_file,
    width = w_all_heatmap,
    height = h_all_heatmap
  )
  
  plot_pseudotime_heatmap(
    mycds_obj[Time_genes, ],
    num_clusters = num_clusters,
    show_rownames = FALSE,
    return_heatmap = FALSE
  )
  
  dev.off()
} else {
  all_heatmap_file <- NA_character_
}

# 4. 绘制Top基因拟时序热图

# 4.1 参数设置

w_top_heatmap <- 8
h_top_heatmap <- 8
name_top_heatmap <- "3.top基因热图"

# 4.2 检查热图聚类数

if (num_clusters > length(top_genes)) {
  stop(
    paste0(
      "num_clusters不能大于Top基因数量。",
      "当前num_clusters为 ",
      num_clusters,
      "，Top基因数量为 ",
      length(top_genes),
      "。"
    )
  )
}

# 4.3 绘制Top基因热图

plot_pseudotime_heatmap(
  mycds_obj[top_genes, ],
  num_clusters = num_clusters,
  show_rownames = TRUE,
  return_heatmap = FALSE
)

# 4.4 保存Top基因热图

top_heatmap_file <- file.path(
  out_dir,
  paste0(name_top_heatmap, ".pdf")
)

pdf(
  file = top_heatmap_file,
  width = w_top_heatmap,
  height = h_top_heatmap
)

plot_pseudotime_heatmap(
  mycds_obj[top_genes, ],
  num_clusters = num_clusters,
  show_rownames = TRUE,
  return_heatmap = FALSE
)

dev.off()

# 5. 绘制Top基因拟时序表达图

# 5.1 参数设置

w_top_expr <- 8
h_top_expr <- 8
name_top_expr <- "4.top基因表达图"

# 5.2 提取Top基因对应的Monocle2对象

top_cds_subset <- mycds_obj[top_genes, ]

# 5.3 绘制按Pseudotime着色的表达趋势图

top_pseudotime_plot <- monocle::plot_genes_in_pseudotime(
  top_cds_subset,
  color_by = "Pseudotime"
) +
  scale_color_viridis_c(option = "viridis")

# 5.4 绘制按State着色的表达趋势图

top_state_plot <- monocle::plot_genes_in_pseudotime(
  top_cds_subset,
  color_by = "State"
)

# 5.5 合并Top基因表达趋势图

top_expr_plot <- top_pseudotime_plot + top_state_plot

print(top_expr_plot)

# 5.6 保存Top基因表达趋势图

top_expr_file <- file.path(
  out_dir,
  paste0(name_top_expr, ".pdf")
)

ggsave(
  filename = top_expr_file,
  plot = top_expr_plot,
  width = w_top_expr,
  height = h_top_expr,
  device = "pdf"
)

# 6. 绘制指定基因的Seurat表达分布图

# 6.1 参数设置

custom_genes_text <- "TFPI2, CCL19, S100A7, TMEM64, ATP9A, PARD6B"

w_feature <- 8
h_feature <- 8
name_feature <- "6.特定基因表达图"

# 6.2 解析指定基因

custom_genes <- trimws(
  unlist(
    strsplit(custom_genes_text, ",")
  )
)

custom_genes <- custom_genes[custom_genes != ""]
custom_genes <- unique(custom_genes)

if (length(custom_genes) == 0) {
  stop("请至少设置一个指定基因。")
}

# 6.3 检查指定基因是否存在于Seurat对象

missing_seurat_genes <- setdiff(
  custom_genes,
  rownames(seurat_obj)
)

if (length(missing_seurat_genes) > 0) {
  stop(
    paste0(
      "以下指定基因不在Spatial_Data中：",
      paste(missing_seurat_genes, collapse = ", ")
    )
  )
}

# 6.4 绘制指定基因表达分布图

feature_plot <- FeaturePlot(
  seurat_obj,
  features = custom_genes
)

print(feature_plot)

# 6.5 保存指定基因表达分布图

feature_file <- file.path(
  out_dir,
  paste0(name_feature, ".pdf")
)

ggsave(
  filename = feature_file,
  plot = feature_plot,
  width = w_feature,
  height = h_feature,
  device = "pdf"
)

# 7. 绘制指定基因拟时序表达图

# 7.1 参数设置

w_custom_expr <- 8
h_custom_expr <- 8
name_custom_expr <- "7.特定基因时序表达"

# 7.2 检查指定基因是否存在于Monocle2对象

missing_monocle_genes <- setdiff(
  custom_genes,
  rownames(mycds_obj)
)

if (length(missing_monocle_genes) > 0) {
  stop(
    paste0(
      "以下指定基因不在mycds中：",
      paste(missing_monocle_genes, collapse = ", ")
    )
  )
}

# 7.3 提取指定基因对应的Monocle2对象

custom_cds_subset <- mycds_obj[custom_genes, ]

# 7.4 绘制按Pseudotime着色的指定基因表达趋势图

custom_pseudotime_plot <- monocle::plot_genes_in_pseudotime(
  custom_cds_subset,
  color_by = "Pseudotime"
) +
  scale_color_viridis_c(option = "viridis")

# 7.5 绘制按State着色的指定基因表达趋势图

custom_state_plot <- monocle::plot_genes_in_pseudotime(
  custom_cds_subset,
  color_by = "State"
)

# 7.6 合并指定基因表达趋势图

custom_expr_plot <- custom_pseudotime_plot + custom_state_plot

print(custom_expr_plot)

# 7.7 保存指定基因表达趋势图

custom_expr_file <- file.path(
  out_dir,
  paste0(name_custom_expr, ".pdf")
)

ggsave(
  filename = custom_expr_file,
  plot = custom_expr_plot,
  width = w_custom_expr,
  height = h_custom_expr,
  device = "pdf"
)

# 8. 绘制指定基因拟时序热图

# 8.1 参数设置

w_custom_heatmap <- 8
h_custom_heatmap <- 8
name_custom_heatmap <- "8.特定基因时序热图"

# 8.2 检查指定基因数量与聚类数

if (num_clusters > length(custom_genes)) {
  stop(
    paste0(
      "num_clusters不能大于指定基因数量。",
      "当前num_clusters为 ",
      num_clusters,
      "，指定基因数量为 ",
      length(custom_genes),
      "。"
    )
  )
}

# 8.3 绘制指定基因拟时序热图

plot_pseudotime_heatmap(
  mycds_obj[custom_genes, ],
  num_clusters = num_clusters,
  show_rownames = TRUE,
  return_heatmap = FALSE
)

# 8.4 保存指定基因拟时序热图

custom_heatmap_file <- file.path(
  out_dir,
  paste0(name_custom_heatmap, ".pdf")
)

pdf(
  file = custom_heatmap_file,
  width = w_custom_heatmap,
  height = h_custom_heatmap
)

plot_pseudotime_heatmap(
  mycds_obj[custom_genes, ],
  num_clusters = num_clusters,
  show_rownames = TRUE,
  return_heatmap = FALSE
)

dev.off()

# 9. 保存参数记录

# 9.1 参数设置

name_params <- "pseudotime_plot_parameters"

# 9.2 生成参数记录

param_text <- paste0(
  "本次拟时序绘图参数总结：\n",
  "- output_folder：", out_dir, "\n",
  "- num_clusters：", num_clusters, "\n",
  "- num_top_genes：", num_top_genes, "\n",
  "- cores_used：", cores_used, "\n",
  "- save_heatmap：", save_heatmap, "\n",
  "- disp.genes数量：", length(disp_genes), "\n",
  "- Time_diff结果行数：", nrow(Time_diff), "\n",
  "- 所有检测基因数量：", length(Time_genes), "\n",
  "- Top基因：", paste(top_genes, collapse = ", "), "\n",
  "- custom_genes：", paste(custom_genes, collapse = ", "), "\n",
  "- 差异基因CSV：", diff_csv_file, "\n",
  "- 本轮作图依赖全局对象：mycds、disp.genes、Spatial_Data\n",
  "说明：\n",
  "1. 使用Monocle2 differentialGeneTest()识别与Pseudotime相关的基因。\n",
  "2. 使用全部检测基因和Top显著基因分别绘制拟时序热图。\n",
  "3. 对Top基因和指定基因绘制Pseudotime及State表达趋势图。\n",
  "4. 使用Seurat FeaturePlot展示指定基因的降维表达分布。"
)

# 9.3 保存参数记录

param_file <- file.path(
  out_dir,
  paste0(name_params, ".txt")
)

writeLines(
  param_text,
  con = param_file
)