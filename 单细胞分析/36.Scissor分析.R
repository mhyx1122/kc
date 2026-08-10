# 1. 加载必要的 R 包

suppressPackageStartupMessages({
  library(Seurat)
  library(Scissor)
  library(limma)
  library(dplyr)
  library(ggplot2)
  library(scCustomize)
})


# 2. 设置输入文件并检查分析对象

bulk_expr_file <- "bulk_expr.csv"
bulk_pheno_file <- "bulk_pheno.csv"

out_dir <- "Scissor分析"

if (!file.exists(bulk_expr_file)) {
  stop(
    paste0(
      "没有找到 bulk 表达矩阵文件：",
      bulk_expr_file
    )
  )
}

if (!file.exists(bulk_pheno_file)) {
  stop(
    paste0(
      "没有找到 bulk 表型文件：",
      bulk_pheno_file
    )
  )
}

if (!exists("seurat", envir = .GlobalEnv)) {
  stop("全局环境中没有 seurat 对象。")
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


# 3. 读取 bulk 表达矩阵

bulk_dataset <- read.csv(
  bulk_expr_file,
  row.names = 1,
  check.names = FALSE
)

rt2 <- bulk_dataset


# 4. 根据 TCGA 样本编号筛选肿瘤样本

# 当所有样本名称长度均大于等于15位时，
# 提取样本编号第14至15位，并保留样本类型编码小于10的样本

if (all(nchar(colnames(rt2)) >= 15)) {
  
  sample_type_code <- as.numeric(
    stringr::str_sub(
      colnames(rt2),
      14,
      15
    )
  )
  
  tumor_sample_index <- which(
    sample_type_code < 10
  )
  
  exp_data_T <- rt2 %>%
    dplyr::select(
      tumor_sample_index
    )
  
} else {
  
  exp_data_T <- rt2
}


# 5. 转置 bulk 表达矩阵并处理样本编号

tumorData <- as.matrix(
  exp_data_T
)

tumorData <- t(
  tumorData
)

if (all(nchar(rownames(tumorData)) > 15)) {
  rownames(tumorData) <- substr(
    rownames(tumorData),
    1,
    12
  )
}

# 对截取样本编号后出现的重复样本进行平均合并

data22 <- limma::avereps(
  tumorData
)


# 6. 读取并整理 bulk 表型数据

phenotype <- read.csv(
  bulk_pheno_file,
  row.names = 1,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

phenotype <- na.omit(
  phenotype
)


# 7. 对齐 bulk 表达矩阵和表型数据中的样本

sameSample <- intersect(
  rownames(data22),
  rownames(phenotype)
)

if (length(sameSample) == 0) {
  stop(
    paste0(
      "bulk 表达矩阵与表型数据没有共同样本，",
      "请检查两个文件中的样本ID。"
    )
  )
}

bulk_data <- data22[
  sameSample,
  ,
  drop = FALSE
]

phenotype_cli <- phenotype[
  sameSample,
  ,
  drop = FALSE
]

bulk_data <- t(
  bulk_data
)

bulk_data <- as.matrix(
  bulk_data
)

storage.mode(bulk_data) <- "double"


# 8. 设置单细胞元数据和降维参数

celltype_col <- "cellType"

reduction <- "umap"

if (!celltype_col %in% colnames(srt@meta.data)) {
  stop(
    paste0(
      "seurat@meta.data 中不存在列：",
      celltype_col
    )
  )
}

if (!reduction %in% names(srt@reductions)) {
  stop(
    paste0(
      "seurat 对象中不存在降维结果：",
      reduction
    )
  )
}


# 9. 将 Seurat RNA assay 转换为 Scissor 适用格式

seurat_V4 <- scCustomize::Convert_Assay(
  seurat_object = srt,
  convert_to = "V3",
  assay = "RNA"
)


# 10. 运行 Scissor 分析

alpha <- 0.05

# 可选：
# family <- "cox"
# family <- "binomial"
# family <- "gaussian"
family <- "cox"

Save_file <- "Scissor_LUAD_survival.RData"

message("正在运行 Scissor 分析，请耐心等待运行完成")

infos1 <- Scissor(
  bulk_dataset = bulk_data,
  sc_dataset = seurat_V4,
  phenotype = phenotype_cli,
  alpha = alpha,
  family = family,
  Save_file = Save_file
)


# 11. 将 Scissor 分类结果写入 Seurat 元数据

Scissor_select <- rep(
  "0",
  ncol(seurat_V4)
)

names(Scissor_select) <- colnames(
  seurat_V4
)

Scissor_select[
  infos1$Scissor_pos
] <- "Scissor+"

Scissor_select[
  infos1$Scissor_neg
] <- "Scissor-"

seurat_V4 <- AddMetaData(
  object = seurat_V4,
  metadata = Scissor_select,
  col.name = "scissor"
)

seurat <- seurat_V4


# 12. 保存每个细胞的 Scissor 分类结果

scissor_cell_result <- data.frame(
  Cell = colnames(seurat_V4),
  Scissor = seurat_V4$scissor,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

write.csv(
  scissor_cell_result,
  file = file.path(
    out_dir,
    "Scissor_细胞分类结果.csv"
  ),
  row.names = FALSE
)


# 13. 生成并保存 Scissor 与细胞类型列联表

df_tile <- as.data.frame(
  table(
    seurat_V4$scissor,
    seurat_V4[[
      celltype_col,
      drop = TRUE
    ]]
  )
)

colnames(df_tile) <- c(
  "Scissor",
  "CellType",
  "Freq"
)

write.csv(
  df_tile,
  file = file.path(
    out_dir,
    "Scissor_CellType_列联表.csv"
  ),
  row.names = FALSE
)


# 14. 生成并保存细胞类型组成比例原始表

df_bar <- as.data.frame(
  table(
    seurat_V4[[
      celltype_col,
      drop = TRUE
    ]],
    seurat_V4$scissor
  )
)

colnames(df_bar) <- c(
  "CellType",
  "Scissor",
  "Freq"
)

write.csv(
  df_bar,
  file = file.path(
    out_dir,
    "CellType_Scissor_组成比例原始表.csv"
  ),
  row.names = FALSE
)


# 15. 设置并整理 Scissor UMAP 参数

pt_size_umap <- 0.1

scissor_order_text <- "Scissor+,Scissor-,0"

palette_umap_text <- paste0(
  "#E64B35FF,#4DBBD5FF,#00A087FF,#3C5488FF,",
  "#F39B7FFF,#8491B4FF,#91D1C2FF,#7E6148FF,",
  "#0073C2FF,#EFC000FF,#868686FF,#CD534CFF,",
  "#7AA6DCFF,#003C67FF,#8F7700FF,#3B3B3BFF,",
  "#A73030FF,#374E55FF,#DF8F44FF,#00A1D5FF,",
  "#B24745FF,#79AF97FF,#6A6599FF,#80796BFF,",
  "#BC3C29FF,#0072B5FF,#E18727FF,#20854EFF,",
  "#7876B1FF,#6F99ADFF,#FFDC91FF,#EE4C97FF"
)

umap_width <- 8
umap_height <- 7
umap_file <- "1.Scissor_UMAP.pdf"

palette_umap <- trimws(
  unlist(
    strsplit(
      palette_umap_text,
      split = ","
    )
  )
)

palette_umap <- palette_umap[
  palette_umap != ""
]

scissor_order <- trimws(
  unlist(
    strsplit(
      scissor_order_text,
      split = ","
    )
  )
)

scissor_order <- scissor_order[
  scissor_order != ""
]

scissor_levels <- unique(
  as.character(
    seurat_V4$scissor
  )
)

order_use <- scissor_order[
  scissor_order %in% scissor_levels
]

remaining_levels <- setdiff(
  scissor_levels,
  order_use
)

final_scissor_levels <- c(
  order_use,
  remaining_levels
)

seurat_V4$scissor <- factor(
  as.character(seurat_V4$scissor),
  levels = final_scissor_levels
)

seurat$scissor <- seurat_V4$scissor

if (length(palette_umap) < length(final_scissor_levels)) {
  stop("Scissor UMAP 的颜色数量不足。")
}

umap_colors <- palette_umap[
  seq_along(final_scissor_levels)
]

names(umap_colors) <- final_scissor_levels


# 16. 绘制并保存 Scissor UMAP

p_umap <- suppressWarnings(
  DimPlot(
    seurat_V4,
    reduction = reduction,
    group.by = "scissor",
    cols = umap_colors,
    pt.size = pt_size_umap,
    raster = FALSE,
    order = final_scissor_levels
  ) +
    theme_classic() +
    theme(
      panel.border = element_rect(
        fill = NA,
        color = "black",
        linewidth = 0.5
      ),
      legend.position = "right"
    )
)

ggsave(
  filename = file.path(
    out_dir,
    umap_file
  ),
  plot = p_umap,
  width = umap_width,
  height = umap_height,
  device = "pdf"
)


# 17. 设置 Scissor 与细胞类型热图参数

tile_low <- "#f7fbff"
tile_high <- "#2171b5"

tile_text_size <- 4
tile_angle <- 45

tile_width <- 10
tile_height <- 6
tile_file <- "2.Scissor_CellType_热图.pdf"


# 18. 绘制并保存 Scissor 与细胞类型热图

p_tile <- ggplot(
  df_tile,
  aes(
    x = CellType,
    y = Scissor,
    fill = Freq
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.5
  ) +
  geom_text(
    aes(label = Freq),
    size = tile_text_size
  ) +
  scale_fill_gradient(
    low = tile_low,
    high = tile_high
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = tile_angle,
      hjust = 1
    ),
    panel.grid = element_blank(),
    axis.title = element_blank()
  )

ggsave(
  filename = file.path(
    out_dir,
    tile_file
  ),
  plot = p_tile,
  width = tile_width,
  height = tile_height,
  device = "pdf"
)


# 19. 设置细胞类型组成比例图参数

bar_angle <- 45

bar_order_text <- "0,Scissor-,Scissor+"

bar_y_title <- "Proportion"

palette_bar_text <- paste0(
  "#E64B35FF,#4DBBD5FF,#00A087FF,#3C5488FF,",
  "#F39B7FFF,#8491B4FF,#91D1C2FF,#7E6148FF,",
  "#0073C2FF,#EFC000FF,#868686FF,#CD534CFF,",
  "#7AA6DCFF,#003C67FF,#8F7700FF,#3B3B3BFF,",
  "#A73030FF,#374E55FF,#DF8F44FF,#00A1D5FF,",
  "#B24745FF,#79AF97FF,#6A6599FF,#80796BFF,",
  "#BC3C29FF,#0072B5FF,#E18727FF,#20854EFF,",
  "#7876B1FF,#6F99ADFF,#FFDC91FF,#EE4C97FF"
)

bar_width <- 10
bar_height <- 6
bar_file <- "3.CellType组成比例.pdf"

palette_bar <- trimws(
  unlist(
    strsplit(
      palette_bar_text,
      split = ","
    )
  )
)

palette_bar <- palette_bar[
  palette_bar != ""
]

bar_order <- trimws(
  unlist(
    strsplit(
      bar_order_text,
      split = ","
    )
  )
)

bar_order <- bar_order[
  bar_order != ""
]

bar_scissor_levels <- unique(
  as.character(
    df_bar$Scissor
  )
)

bar_order_use <- bar_order[
  bar_order %in% bar_scissor_levels
]

bar_remaining_levels <- setdiff(
  bar_scissor_levels,
  bar_order_use
)

bar_final_levels <- c(
  bar_order_use,
  bar_remaining_levels
)

df_bar$Scissor <- factor(
  as.character(df_bar$Scissor),
  levels = bar_final_levels
)

if (length(palette_bar) < length(bar_final_levels)) {
  stop("细胞类型组成比例图的颜色数量不足。")
}

bar_colors <- palette_bar[
  seq_along(bar_final_levels)
]

names(bar_colors) <- bar_final_levels


# 20. 绘制并保存细胞类型组成比例图

p_bar <- ggplot(
  df_bar,
  aes(
    x = CellType,
    y = Freq,
    fill = Scissor
  )
) +
  geom_bar(
    stat = "identity",
    position = "fill",
    width = 0.8
  ) +
  scale_fill_manual(
    values = bar_colors,
    drop = FALSE
  ) +
  scale_y_continuous(
    labels = scales::percent
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = bar_angle,
      hjust = 1
    ),
    panel.grid = element_blank()
  ) +
  labs(
    x = NULL,
    y = bar_y_title
  )

ggsave(
  filename = file.path(
    out_dir,
    bar_file
  ),
  plot = p_bar,
  width = bar_width,
  height = bar_height,
  device = "pdf"
)


# 21. 保存本次 Scissor 分析参数

params_file <- "Scissor_parameters.txt"

scissor_counts <- table(
  seurat_V4$scissor
)

scissor_counts_text <- paste(
  names(scissor_counts),
  as.integer(scissor_counts),
  sep = ": ",
  collapse = "\n"
)

param_text <- paste0(
  "运行流程说明：\n",
  "1. 读取 bulk 表达矩阵与表型数据。\n",
  "2. 当 bulk 样本名长度均大于等于15位时，",
  "按 TCGA 样本命名规则提取第14至15位样本类型编码，",
  "并保留编码小于10的样本。\n",
  "3. 将 bulk 表达矩阵转置；当样本名长度大于15位时，",
  "将样本名截断为前12位用于后续与表型数据匹配。\n",
  "4. 对截取样本编号后出现的重复样本进行平均合并，",
  "并与表型数据按照共同样本ID进行对齐。\n",
  "5. 将 Seurat 对象中的 RNA assay 转换为适用于 Scissor 的 V3 格式。\n",
  "6. 使用 Scissor 将 bulk 表型信号映射到单细胞数据中，",
  "识别 Scissor+、Scissor- 和未归类细胞。\n",
  "7. 将 Scissor 结果写入 Seurat 元数据中的 scissor 列，",
  "并输出 UMAP、列联热图和组成比例图。\n",
  "8. 保存细胞分类结果表、列联表、组成比例原始表和参数记录。\n",
  "\n",
  "本次分析参数总结：\n",
  "- bulk表达矩阵：", bulk_expr_file, "\n",
  "- bulk表型文件：", bulk_pheno_file, "\n",
  "- alpha：", alpha, "\n",
  "- family：", family, "\n",
  "- Save_file：", Save_file, "\n",
  "- reduction：", reduction, "\n",
  "- cellType列：", celltype_col, "\n",
  "- bulk与表型共同样本数：", length(sameSample), "\n",
  "- 单细胞总数：", ncol(seurat_V4), "\n",
  "- UMAP点大小：", pt_size_umap, "\n",
  "- UMAP显示顺序：",
  paste(final_scissor_levels, collapse = ", "),
  "\n",
  "- 比例图堆叠顺序：",
  paste(bar_final_levels, collapse = ", "),
  "\n",
  "- Scissor结果统计：\n",
  scissor_counts_text,
  "\n",
  "\n",
  "写作提示词：\n",
  "1. 将 bulk 临床或表型信息与单细胞转录组数据进行整合，",
  "利用 Scissor 识别与目标表型显著相关的细胞亚群。\n",
  "2. 通过降维图展示 Scissor+、Scissor- 和未归类细胞",
  "在单细胞图谱中的分布特征。\n",
  "3. 结合细胞类型列联热图和组成比例图，",
  "评估目标表型相关细胞在不同细胞群中的富集模式。\n"
)

writeLines(
  text = param_text,
  con = file.path(
    out_dir,
    params_file
  )
)

message(
  "Scissor 分析完成，结果已保存至：",
  out_dir
)