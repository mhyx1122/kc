# 1. 加载必要的 R 包

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(Seurat)
  library(monocle)
})


# 2. 检查分析所需对象并创建输出文件夹

out_dir <- "Monocle2结果"

if (!exists("mycds", envir = .GlobalEnv)) {
  stop("全局环境中没有找到 mycds 对象，请先运行上一套 Monocle2 拟时序分析代码。")
}

if (!exists("disp.genes", envir = .GlobalEnv)) {
  stop("全局环境中没有找到 disp.genes 对象，请先运行上一套 Monocle2 拟时序分析代码。")
}

if (!exists("seurat", envir = .GlobalEnv)) {
  stop("全局环境中没有找到 seurat 对象，请先加载 Seurat 对象。")
}

mycds_obj <- get("mycds", envir = .GlobalEnv)
disp_genes <- get("disp.genes", envir = .GlobalEnv)
seurat_obj <- get("seurat", envir = .GlobalEnv)

if (!inherits(seurat_obj, "Seurat")) {
  stop("全局环境中的 seurat 不是 Seurat 对象。")
}

if (length(disp_genes) == 0) {
  stop("disp.genes 中没有轨迹构建基因。")
}

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}


# 3. 进行拟时序差异基因分析

cores_used <- 1

diff_csv_file <- "Pseudotime_Differential_Genes.csv"

message("正在进行拟时序差异基因分析，请耐心等待运行完毕")

Time_diff <- differentialGeneTest(
  mycds_obj[disp_genes, ],
  cores = cores_used,
  fullModelFormulaStr = "~sm.ns(Pseudotime)"
)

if (!"gene_short_name" %in% colnames(Time_diff)) {
  stop("differentialGeneTest() 结果中不存在 gene_short_name 列。")
}

if (!"qval" %in% colnames(Time_diff)) {
  stop("differentialGeneTest() 结果中不存在 qval 列。")
}

write.csv(
  Time_diff,
  file = file.path(out_dir, diff_csv_file),
  row.names = TRUE
)

message(
  "拟时序差异基因结果已保存至：",
  file.path(out_dir, diff_csv_file)
)


# 4. 提取拟时序差异分析中的全部基因

Time_genes <- Time_diff %>%
  dplyr::pull(gene_short_name) %>%
  as.character()

Time_genes <- Time_genes[
  !is.na(Time_genes) &
    Time_genes != ""
]

if (length(Time_genes) == 0) {
  stop("差异分析结果中没有可用于绘图的基因。")
}


# 5. 绘制并保存所有差异基因拟时序热图

save_heatmap <- TRUE

num_clusters <- 4

all_heatmap_width <- 8
all_heatmap_height <- 8
all_heatmap_file <- "2.所有基因的热图.pdf"

if (save_heatmap) {
  message("正在绘制所有差异基因拟时序热图")
  
  pdf(
    file = file.path(out_dir, all_heatmap_file),
    width = all_heatmap_width,
    height = all_heatmap_height
  )
  
  plot_pseudotime_heatmap(
    mycds_obj[Time_genes, ],
    num_clusters = num_clusters,
    show_rownames = FALSE,
    return_heatmap = FALSE
  )
  
  dev.off()
  
  message(
    "所有差异基因拟时序热图已保存至：",
    file.path(out_dir, all_heatmap_file)
  )
} else {
  message("save_heatmap 为 FALSE，已跳过所有差异基因拟时序热图。")
}


# 6. 筛选拟时序差异最显著的 Top 基因

num_top_genes <- 8

top_genes <- Time_diff %>%
  dplyr::slice_min(
    order_by = qval,
    n = num_top_genes
  ) %>%
  dplyr::pull(gene_short_name) %>%
  as.character()

top_genes <- top_genes[
  !is.na(top_genes) &
    top_genes != ""
]

if (length(top_genes) == 0) {
  stop("没有筛选到可用于绘图的 Top 基因。")
}


# 7. 绘制并保存 Top 基因拟时序热图

top_heatmap_width <- 8
top_heatmap_height <- 8
top_heatmap_file <- "3.top基因热图.pdf"

message("正在绘制 Top 基因拟时序热图")

pdf(
  file = file.path(out_dir, top_heatmap_file),
  width = top_heatmap_width,
  height = top_heatmap_height
)

plot_pseudotime_heatmap(
  mycds_obj[top_genes, ],
  num_clusters = num_clusters,
  show_rownames = TRUE,
  return_heatmap = FALSE
)

dev.off()

message(
  "Top 基因拟时序热图已保存至：",
  file.path(out_dir, top_heatmap_file)
)


# 8. 绘制并保存 Top 基因时序表达图

top_expr_width <- 8
top_expr_height <- 8
top_expr_file <- "4.top基因表达图.pdf"

top_genes_cds <- mycds_obj[top_genes, ]

p_top_pseudotime <- monocle::plot_genes_in_pseudotime(
  top_genes_cds,
  color_by = "Pseudotime"
) +
  scale_color_viridis_c(option = "viridis")

p_top_state <- monocle::plot_genes_in_pseudotime(
  top_genes_cds,
  color_by = "State"
)

p_top_expr <- p_top_pseudotime + p_top_state

pdf(
  file = file.path(out_dir, top_expr_file),
  width = top_expr_width,
  height = top_expr_height
)

print(p_top_expr)

dev.off()

message(
  "Top 基因时序表达图已保存至：",
  file.path(out_dir, top_expr_file)
)


# 9. 设置并检查指定基因

custom_genes_text <- "TFPI2, CCL19, S100A7, TMEM64, ATP9A, PARD6B"

custom_genes <- trimws(
  unlist(
    strsplit(
      custom_genes_text,
      split = ","
    )
  )
)

custom_genes <- custom_genes[custom_genes != ""]

if (length(custom_genes) == 0) {
  stop("没有设置指定基因。")
}

missing_genes_seurat <- setdiff(
  custom_genes,
  rownames(seurat_obj)
)

if (length(missing_genes_seurat) > 0) {
  stop(
    paste0(
      "以下指定基因不存在于 seurat 对象中：",
      paste(missing_genes_seurat, collapse = ", ")
    )
  )
}

missing_genes_monocle <- setdiff(
  custom_genes,
  rownames(mycds_obj)
)

if (length(missing_genes_monocle) > 0) {
  stop(
    paste0(
      "以下指定基因不存在于 mycds 对象中：",
      paste(missing_genes_monocle, collapse = ", ")
    )
  )
}


# 10. 绘制并保存指定基因 FeaturePlot

feature_width <- 8
feature_height <- 8
feature_file <- "6.特定基因表达图.pdf"

p_feature <- FeaturePlot(
  seurat_obj,
  features = custom_genes
)

pdf(
  file = file.path(out_dir, feature_file),
  width = feature_width,
  height = feature_height
)

print(p_feature)

dev.off()

message(
  "指定基因表达图已保存至：",
  file.path(out_dir, feature_file)
)


# 11. 绘制并保存指定基因时序表达图

custom_expr_width <- 8
custom_expr_height <- 8
custom_expr_file <- "7.特定基因时序表达.pdf"

custom_genes_cds <- mycds_obj[custom_genes, ]

p_custom_pseudotime <- monocle::plot_genes_in_pseudotime(
  custom_genes_cds,
  color_by = "Pseudotime"
) +
  scale_color_viridis_c(option = "viridis")

p_custom_state <- monocle::plot_genes_in_pseudotime(
  custom_genes_cds,
  color_by = "State"
)

p_custom_expr <- p_custom_pseudotime + p_custom_state

pdf(
  file = file.path(out_dir, custom_expr_file),
  width = custom_expr_width,
  height = custom_expr_height
)

print(p_custom_expr)

dev.off()

message(
  "指定基因时序表达图已保存至：",
  file.path(out_dir, custom_expr_file)
)


# 12. 绘制并保存指定基因拟时序热图

custom_heatmap_width <- 8
custom_heatmap_height <- 8
custom_heatmap_file <- "8.特定基因时序热图.pdf"

message("正在绘制指定基因拟时序热图")

pdf(
  file = file.path(out_dir, custom_heatmap_file),
  width = custom_heatmap_width,
  height = custom_heatmap_height
)

plot_pseudotime_heatmap(
  mycds_obj[custom_genes, ],
  num_clusters = num_clusters,
  show_rownames = TRUE,
  return_heatmap = FALSE
)

dev.off()

message(
  "指定基因拟时序热图已保存至：",
  file.path(out_dir, custom_heatmap_file)
)


# 13. 保存本次拟时序绘图参数

params_file <- "pseudotime_plot_parameters.txt"

param_text <- paste0(
  "本次拟时序绘图参数总结：\n",
  "- output_folder：", out_dir, "\n",
  "- num_clusters：", num_clusters, "\n",
  "- num_top_genes：", num_top_genes, "\n",
  "- cores_used：", cores_used, "\n",
  "- save_heatmap：", save_heatmap, "\n",
  "- custom_genes：", paste(custom_genes, collapse = ", "), "\n",
  "- 差异基因CSV：", file.path(out_dir, diff_csv_file), "\n",
  "- 本轮作图依赖全局对象：mycds / disp.genes / seurat\n",
  "说明：\n",
  "1. Monocle2基因的可视化。\n"
)

writeLines(
  text = param_text,
  con = file.path(out_dir, params_file)
)

message(
  "拟时序绘图参数已保存至：",
  file.path(out_dir, params_file)
)