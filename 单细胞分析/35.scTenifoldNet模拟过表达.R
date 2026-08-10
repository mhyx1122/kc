# 1. 加载必要的 R 包

suppressPackageStartupMessages({
  library(Seurat)
  library(scTenifoldNet)
  library(dplyr)
  library(ggplot2)
  library(MASS)
})


# 2. 检查 Seurat 对象并创建输出文件夹

out_dir <- "scTenifoldNet_模拟过表达"

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


# 4. 设置目标基因和基因过滤参数

target_gene <- "FLT1"

# 可选：
# filter_method <- "pct"
# filter_method <- "hvg"
filter_method <- "hvg"

# filter_method = "pct" 时使用
min_pct <- 0.05

# filter_method = "hvg" 时使用
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


# 5. 根据设置过滤基因

if (filter_method == "pct") {
  
  gene_keep_idx <- rowMeans(raw_counts != 0) >= min_pct
  
  expr_mat <- raw_counts[
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
  
  expr_mat <- raw_counts[
    hvg_genes,
    ,
    drop = FALSE
  ]
  
} else {
  stop("filter_method 只能设置为 \"pct\" 或 \"hvg\"。")
}


# 6. 检查并补回目标基因

if (!(target_gene %in% rownames(expr_mat))) {
  expr_mat <- raw_counts[
    unique(
      c(
        rownames(expr_mat),
        target_gene
      )
    ),
    ,
    drop = FALSE
  ]
}

message(
  "进入网络分析的基因数量：",
  nrow(expr_mat)
)

message(
  "进入网络分析的细胞数量：",
  ncol(expr_mat)
)


# 7. 构建野生型基因调控网络

nNet <- 10
nCells <- 500
nComp <- 3
nCores <- 8
random_seed <- 1234

nCells_use <- min(
  nCells,
  ncol(expr_mat)
)

set.seed(random_seed)

message("正在构建野生型基因调控网络")

wt_net_list <- makeNetworks(
  X = expr_mat,
  nNet = nNet,
  nCells = nCells_use,
  nComp = nComp,
  nCores = nCores
)


# 8. 对野生型网络进行张量分解

K <- 3

message("正在对野生型网络进行张量分解")

wt_tensor_fit <- tensorDecomposition(
  xList = wt_net_list,
  K = K
)

wt_grn <- wt_tensor_fit[[1]]


# 9. 构建目标基因模拟过表达网络

amplify_factor <- 10

if (!(target_gene %in% rownames(wt_grn))) {
  stop(
    paste0(
      "目标基因 ",
      target_gene,
      " 不存在于张量分解后的基因调控网络中。"
    )
  )
}

perturb_grn <- wt_grn

perturb_grn[target_gene, ] <- perturb_grn[target_gene, ] * amplify_factor


# 10. 对野生型网络和模拟过表达网络进行流形对齐

message("正在进行野生型网络与模拟过表达网络的流形对齐")

ma_embed <- manifoldAlignment(
  t(as.matrix(wt_grn)),
  t(as.matrix(perturb_grn))
)


# 11. 从流形对齐结果中提取基因名称

all_node_names <- rownames(ma_embed)

ref_genes <- all_node_names[
  grepl("^X_", all_node_names)
]

ref_genes <- gsub(
  "^X_",
  "",
  ref_genes
)

n_gene <- length(ref_genes)

if (n_gene == 0) {
  stop("未能从 manifoldAlignment 结果中解析到基因名称。")
}


# 12. 计算每个基因在两个网络之间的位移距离

gene_shift <- numeric(n_gene)

for (i in seq_len(n_gene)) {
  gene_shift[i] <- as.numeric(
    dist(
      rbind(
        ma_embed[i, ],
        ma_embed[i + n_gene, ]
      )
    )
  )
}


# 13. 对基因位移距离进行 Box-Cox 转换和标准化

positive_shift <- gene_shift[
  gene_shift > 0
]

bc_lambda <- try(
  {
    bc_obj <- boxcox(
      positive_shift ~ 1,
      plotit = FALSE
    )
    
    bc_obj$x[
      which.max(bc_obj$y)
    ]
  },
  silent = TRUE
)

if (
  inherits(bc_lambda, "try-error") ||
  length(positive_shift) == 0
) {
  
  shift_transformed <- gene_shift
  
} else {
  
  if (abs(bc_lambda) < 1e-6) {
    
    shift_transformed <- log(
      gene_shift + 1e-8
    )
    
  } else {
    
    shift_transformed <- (
      (gene_shift + 1e-8)^bc_lambda - 1
    ) / bc_lambda
  }
}

z_metric <- as.numeric(
  scale(shift_transformed)
)


# 14. 计算扰动统计量和显著性

background_idx <- which(
  ref_genes != target_gene
)

bg_mean_sq <- mean(
  gene_shift[background_idx]^2
)

score_ratio <- (
  gene_shift^2
) / bg_mean_sq

p_raw <- pchisq(
  q = score_ratio,
  df = 1,
  lower.tail = FALSE
)

p_fdr <- p.adjust(
  p_raw,
  method = "fdr"
)


# 15. 整理模拟过表达结果

result_tbl <- data.frame(
  gene = ref_genes,
  shift_distance = as.numeric(gene_shift),
  z_score = z_metric,
  FC = score_ratio,
  p_value = p_raw,
  p_adj = p_fdr,
  stringsAsFactors = FALSE
) %>%
  dplyr::arrange(desc(z_score))

top_20_genes <- head(
  result_tbl,
  20
)


# 16. 保存网络扰动结果

result_csv_file <- "02_网络扰动结果.csv"

write.csv(
  result_tbl,
  file = file.path(
    out_dir,
    result_csv_file
  ),
  row.names = FALSE
)

message(
  "网络扰动结果已保存至：",
  file.path(out_dir, result_csv_file)
)


# 17. 设置网络扰动图参数

x_var <- "z_score"
y_var <- "p_adj"

use_neglog10_y <- TRUE
include_target_gene_in_plot <- TRUE

top_label_n <- 15

custom_label_genes_text <- ""

perturb_width <- 8
perturb_height <- 7
perturb_file <- "01_网络扰动图.pdf"


# 18. 准备网络扰动图数据

plot_df <- result_tbl

if (!include_target_gene_in_plot) {
  plot_df <- subset(
    plot_df,
    gene != target_gene
  )
}

sig_col <- if (y_var == "p_value") {
  "p_value"
} else {
  "p_adj"
}

plot_df <- plot_df %>%
  dplyr::arrange(desc(z_score)) %>%
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


# 19. 绘制网络扰动图

p_perturb <- ggplot(
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
    subtitle = paste0(
      "Outgoing edge weights amplified by ",
      amplify_factor,
      "×"
    ),
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
  y_var %in% c("p_adj", "p_value")
) {
  p_perturb <- p_perturb +
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed",
      color = "#B7B7B7",
      linewidth = 0.6
    )
}


# 20. 保存网络扰动图

ggsave(
  filename = file.path(
    out_dir,
    perturb_file
  ),
  plot = p_perturb,
  width = perturb_width,
  height = perturb_height,
  device = "pdf"
)

message(
  "网络扰动图已保存至：",
  file.path(out_dir, perturb_file)
)


# 21. 保存本次模拟过表达分析参数

params_file <- "scTenifoldNet_OE_参数记录.txt"

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

sig_rule <- if (y_var == "p_value") {
  "p_value < 0.05"
} else {
  "p_adj < 0.05"
}

param_text <- paste0(
  "本次 scTenifoldNet 模拟过表达分析参数总结：\n",
  "- 目标基因：", target_gene, "\n",
  "- 基因过滤方式：", filter_desc, "\n",
  "- 进入分析的基因数量：", nrow(expr_mat), "\n",
  "- 进入分析的细胞数量：", ncol(expr_mat), "\n",
  "- 出边放大倍数 amplify_factor：", amplify_factor, "\n",
  "- nNet：", nNet, "\n",
  "- nCells（实际使用）：", nCells_use, "\n",
  "- nComp：", nComp, "\n",
  "- K：", K, "\n",
  "- nCores：", nCores, "\n",
  "- 随机种子：", random_seed, "\n",
  "- X轴变量：", x_var, "\n",
  "- Y轴变量：", y_var, "\n",
  "- Y轴是否进行 -log10 转换：",
  ifelse(use_neglog10_y, "是", "否"),
  "\n",
  "- Significant 判定规则：", sig_rule, "\n",
  "- 作图是否纳入目标基因本身：",
  ifelse(include_target_gene_in_plot, "是", "否"),
  "\n",
  "- 默认标记前N基因：", top_label_n, "\n",
  "- 额外标注基因：",
  ifelse(
    custom_label_genes_text == "",
    "无",
    custom_label_genes_text
  ),
  "\n",
  "- 输出文件夹：", out_dir, "\n",
  "- 网络扰动结果：",
  file.path(out_dir, result_csv_file),
  "\n",
  "- 网络扰动图：",
  file.path(out_dir, perturb_file),
  "\n",
  "\n",
  "运行说明：\n",
  "1. 从 Seurat 对象的 RNA assay 中提取 counts 矩阵。\n",
  "2. 根据设置进行表达细胞比例过滤或高变基因过滤。\n",
  "3. 如果目标基因不在过滤后的矩阵中，则从原始 counts 矩阵中补回。\n",
  "4. 使用 makeNetworks() 构建野生型基因调控网络集合。\n",
  "5. 使用 tensorDecomposition() 进行张量分解并获得去噪网络。\n",
  "6. 将目标基因的出边权重放大指定倍数，构建模拟过表达网络。\n",
  "7. 使用 manifoldAlignment() 对野生型网络和模拟过表达网络进行流形对齐。\n",
  "8. 计算基因位移距离、标准化扰动分数、统计量和显著性。\n",
  "9. 保存全基因扰动结果和网络扰动图。\n"
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