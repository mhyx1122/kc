# 1. 加载必要的 R 包

suppressPackageStartupMessages({
  library(Seurat)
  library(scTenifoldKnk)
  library(dplyr)
  library(ggplot2)
})


# 2. 检查 Seurat 对象并创建输出文件夹

out_dir <- "scTenifoldKnk_模拟敲除"

if (!exists("seurat", envir = .GlobalEnv)) {
  stop("全局环境中没有名为 seurat 的对象。")
}

seurat_obj <- get("seurat", envir = .GlobalEnv)

if (!inherits(seurat_obj, "Seurat")) {
  stop("全局环境中的 seurat 不是 Seurat 对象。")
}

if (!"RNA" %in% names(seurat_obj@assays)) {
  stop("seurat 对象中不存在 RNA assay。")
}

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}


# 3. 提取 RNA counts 矩阵

raw_counts <- tryCatch(
  {
    LayerData(
      object = seurat_obj,
      assay = "RNA",
      layer = "counts"
    )
  },
  error = function(e) {
    NULL
  }
)

if (is.null(raw_counts)) {
  raw_counts <- tryCatch(
    {
      GetAssayData(
        object = seurat_obj,
        assay = "RNA",
        slot = "counts"
      )
    },
    error = function(e) {
      NULL
    }
  )
}

if (is.null(raw_counts)) {
  stop("无法从 Seurat 对象中提取 RNA counts 矩阵。")
}


# 4. 设置敲除基因和基因过滤参数

target_gene <- "FLT1"

# 可选：
# filter_method <- "pct"
# filter_method <- "hvg"
filter_method <- "hvg"

# 当 filter_method = "pct" 时使用
min_pct <- 0.05

# 当 filter_method = "hvg" 时使用
top_hvg_n <- 2000

if (!(target_gene %in% rownames(raw_counts))) {
  stop(
    paste0(
      "目标基因 ",
      target_gene,
      " 在原始 counts 矩阵中不存在。"
    )
  )
}


# 5. 根据设定的方式过滤基因

if (filter_method == "pct") {
  
  gene_keep_idx <- rowMeans(raw_counts != 0) >= min_pct
  
  countMatrix <- raw_counts[
    gene_keep_idx,
    ,
    drop = FALSE
  ]
  
} else if (filter_method == "hvg") {
  
  seurat_tmp <- seurat_obj
  
  DefaultAssay(seurat_tmp) <- "RNA"
  
  seurat_tmp <- FindVariableFeatures(
    object = seurat_tmp,
    selection.method = "vst",
    nfeatures = top_hvg_n,
    verbose = FALSE
  )
  
  hvg_genes <- VariableFeatures(seurat_tmp)
  
  hvg_genes <- intersect(
    hvg_genes,
    rownames(raw_counts)
  )
  
  countMatrix <- raw_counts[
    hvg_genes,
    ,
    drop = FALSE
  ]
  
} else {
  stop("filter_method 只能设置为 \"pct\" 或 \"hvg\"。")
}


# 6. 检查并补回目标基因

if (!(target_gene %in% rownames(countMatrix))) {
  countMatrix <- raw_counts[
    unique(
      c(
        rownames(countMatrix),
        target_gene
      )
    ),
    ,
    drop = FALSE
  ]
}

message(
  "进入 scTenifoldKnk 分析的基因数量：",
  nrow(countMatrix)
)

message(
  "进入 scTenifoldKnk 分析的细胞数量：",
  ncol(countMatrix)
)


# 7. 运行 scTenifoldKnk 模拟敲除分析

nc_nNet <- 10
nc_nCells <- 100
nc_nComp <- 3
td_K <- 3
nCores <- 8
random_seed <- 1234

set.seed(random_seed)

message(
  "正在运行 ",
  target_gene,
  " 的 scTenifoldKnk 模拟敲除分析"
)

knk_result <- scTenifoldKnk(
  countMatrix = countMatrix,
  qc = FALSE,
  gKO = target_gene,
  nc_nNet = nc_nNet,
  nc_nCells = nc_nCells,
  nc_nComp = nc_nComp,
  td_K = td_K,
  nCores = nCores
)

if (is.null(knk_result$diffRegulation)) {
  stop("scTenifoldKnk 未返回 diffRegulation 结果。")
}


# 8. 整理模拟敲除结果

result_tbl <- as.data.frame(
  knk_result$diffRegulation,
  stringsAsFactors = FALSE
)

needed_cols <- c(
  "gene",
  "distance",
  "Z",
  "FC",
  "p.value",
  "p.adj"
)

miss_cols <- setdiff(
  needed_cols,
  colnames(result_tbl)
)

if (length(miss_cols) > 0) {
  stop(
    paste0(
      "diffRegulation 缺少以下列：",
      paste(miss_cols, collapse = ", ")
    )
  )
}

result_tbl <- result_tbl %>%
  dplyr::arrange(desc(FC))


# 9. 提取 Top 基因结果

numGene <- 20

top_tbl <- head(
  result_tbl,
  numGene
)


# 10. 保存全部结果和 Top 基因结果

name_result_csv <- "total_genes"
name_top_csv <- "top_genes"

total_file <- file.path(
  out_dir,
  paste0(name_result_csv, ".csv")
)

top_file <- file.path(
  out_dir,
  paste0(name_top_csv, ".csv")
)

write.csv(
  result_tbl,
  file = total_file,
  row.names = FALSE
)

write.csv(
  top_tbl,
  file = top_file,
  row.names = FALSE
)

message(
  "全部模拟敲除结果已保存至：",
  total_file
)

message(
  "Top 基因结果已保存至：",
  top_file
)


# 11. 设置网络扰动图参数

x_var <- "Z"
y_var <- "p.adj"

use_neglog10_y <- TRUE
include_target_gene <- TRUE

top_label_n <- 15

custom_label_genes_text <- ""

knk_width <- 8
knk_height <- 7
knk_file <- "01_网络扰动图.pdf"


# 12. 准备网络扰动图数据

plot_df <- result_tbl

if (!include_target_gene) {
  plot_df <- subset(
    plot_df,
    gene != target_gene
  )
}

sig_col <- if (y_var == "p.value") {
  "p.value"
} else {
  "p.adj"
}

plot_df <- plot_df %>%
  dplyr::arrange(desc(FC)) %>%
  dplyr::mutate(
    sig_group = ifelse(
      .data[[sig_col]] < 0.05,
      "Significant",
      "Not Significant"
    )
  )

top_n_genes <- head(
  plot_df$gene,
  top_label_n
)

custom_label_genes <- trimws(
  unlist(
    strsplit(
      custom_label_genes_text,
      split = ","
    )
  )
)

custom_label_genes <- custom_label_genes[
  custom_label_genes != ""
]

label_genes <- unique(
  c(
    top_n_genes,
    custom_label_genes
  )
)

if (include_target_gene) {
  label_genes <- unique(
    c(
      target_gene,
      label_genes
    )
  )
}

plot_df$gene_label <- ifelse(
  plot_df$gene %in% label_genes,
  plot_df$gene,
  NA
)

plot_df$x_plot <- plot_df[[x_var]]

if (use_neglog10_y) {
  plot_df$y_plot <- -log10(
    plot_df[[y_var]] + 1e-300
  )
} else {
  plot_df$y_plot <- plot_df[[y_var]]
}

x_lab <- x_var

if (use_neglog10_y) {
  y_lab <- paste0(
    "-log10(",
    y_var,
    ")"
  )
} else {
  y_lab <- y_var
}


# 13. 绘制网络扰动图

p_knk <- ggplot(
  plot_df,
  aes(
    x = x_plot,
    y = y_plot
  )
) +
  geom_point(
    aes(fill = sig_group),
    shape = 21,
    color = "white",
    size = 2.8,
    alpha = 0.9,
    stroke = 0.3
  ) +
  geom_text(
    data = subset(
      plot_df,
      !is.na(gene_label)
    ),
    aes(label = gene_label),
    size = 3.8,
    color = "#5C5C5C",
    vjust = -0.6,
    check_overlap = TRUE
  ) +
  scale_fill_manual(
    values = c(
      "Significant" = "#F4A7B9",
      "Not Significant" = "#BFD7EA"
    )
  ) +
  labs(
    title = paste0(
      target_gene,
      " Network Perturbation Analysis"
    ),
    subtitle = "scTenifoldKnk simulated knockout",
    x = x_lab,
    y = y_lab,
    fill = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      color = "#4F4F4F",
      size = 16
    ),
    plot.subtitle = element_text(
      hjust = 0.5,
      color = "#7A7A7A",
      size = 12
    ),
    axis.title = element_text(
      color = "#5A5A5A",
      face = "bold"
    ),
    axis.text = element_text(
      color = "#6B6B6B"
    ),
    legend.position = "top",
    panel.grid = element_blank(),
    axis.line = element_line(
      color = "#7A7A7A",
      linewidth = 0.6
    ),
    axis.ticks = element_line(
      color = "#7A7A7A",
      linewidth = 0.5
    ),
    plot.background = element_rect(
      fill = "#FFFDFB",
      color = NA
    ),
    panel.background = element_rect(
      fill = "#FFFDFB",
      color = NA
    )
  )

if (
  use_neglog10_y &&
  y_var %in% c("p.adj", "p.value")
) {
  p_knk <- p_knk +
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed",
      color = "#B7B7B7",
      linewidth = 0.6
    )
}


# 14. 保存网络扰动图

ggsave(
  filename = file.path(
    out_dir,
    knk_file
  ),
  plot = p_knk,
  width = knk_width,
  height = knk_height,
  device = "pdf"
)

message(
  "网络扰动图已保存至：",
  file.path(out_dir, knk_file)
)


# 15. 保存本次模拟敲除分析参数

params_file <- "scTenifoldKnk_参数记录.txt"

if (filter_method == "pct") {
  filter_desc <- paste0(
    "按表达细胞比例过滤（至少在 ",
    min_pct * 100,
    "% 细胞中表达）"
  )
} else {
  filter_desc <- paste0(
    "按高变基因过滤（前 ",
    top_hvg_n,
    " 个高变基因）"
  )
}

sig_rule <- if (y_var == "p.value") {
  "p.value < 0.05"
} else {
  "p.adj < 0.05"
}

param_text <- paste0(
  "本次 scTenifoldKnk 模拟敲除分析参数总结：\n",
  "- 敲除基因：", target_gene, "\n",
  "- 基因过滤方式：", filter_desc, "\n",
  "- 进入分析的基因数量：", nrow(countMatrix), "\n",
  "- 进入分析的细胞数量：", ncol(countMatrix), "\n",
  "- qc：FALSE（固定）\n",
  "- nc_nNet：", nc_nNet, "\n",
  "- nc_nCells：", nc_nCells, "\n",
  "- nc_nComp：", nc_nComp, "\n",
  "- td_K：", td_K, "\n",
  "- nCores：", nCores, "\n",
  "- 随机种子：", random_seed, "\n",
  "- X轴变量：", x_var, "\n",
  "- Y轴变量：", y_var, "\n",
  "- Y轴是否进行 -log10 转换：",
  ifelse(use_neglog10_y, "是", "否"),
  "\n",
  "- Significant 判定规则：", sig_rule, "\n",
  "- 作图是否纳入目标基因：",
  ifelse(include_target_gene, "是", "否"),
  "\n",
  "- 默认标注前N基因：", top_label_n, "\n",
  "- 额外标注基因：",
  ifelse(
    custom_label_genes_text == "",
    "无",
    custom_label_genes_text
  ),
  "\n",
  "- top_genes 导出数量：", numGene, "\n",
  "- 输出文件夹：", out_dir, "\n",
  "- 全部结果文件：", total_file, "\n",
  "- Top结果文件：", top_file, "\n",
  "- 网络扰动图：", file.path(out_dir, knk_file), "\n",
  "\n",
  "运行说明：\n",
  "1. 从 Seurat 对象的 RNA assay 中提取 counts 矩阵。\n",
  "2. 根据设定的方式进行表达比例过滤或高变基因过滤。\n",
  "3. 如果目标基因不在过滤后的矩阵中，则从原始 counts 矩阵中补回。\n",
  "4. 使用 scTenifoldKnk() 对目标基因进行模拟敲除分析。\n",
  "5. 提取 diffRegulation 结果，并按照 FC 从高到低排序。\n",
  "6. 保存全部基因结果、Top 基因结果和网络扰动图。\n"
)

writeLines(
  text = param_text,
  con = file.path(
    out_dir,
    params_file
  )
)

message(
  "分析参数已保存至：",
  file.path(out_dir, params_file)
)