suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(scop)
  library(grid)
  library(RColorBrewer)
})

# 1. 检查全局环境中是否存在 seurat 对象
if (!exists("seurat", envir = .GlobalEnv)) {
  stop("全局环境中没有 seurat 对象")
}

srt <- get("seurat", envir = .GlobalEnv)


# 2. 输出目录模块

# 2.1 输出目录参数
out_dir <- "4.5注释图美化"

# 2.2 创建输出目录
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}


# 3. GroupHeatmap 参数模块

# 3.1 特征基因参数
features_text <- "ICOS, GZMK, CD6, BCL11B, TRAT1, CLEC10A, TLR8, CLEC4E, FCN1, LGALS2, ADGRL4, CLEC14A, MYCT1, SOX17, NOVA2, FHL5, OR51E1, NGF, LRRC10B, LYPD1, SAA1, KRT5, ELF5, ANKRD30A, KRT17, SFRP2, COL10A1, DPT, OMD, CILP, MKI67, GTSE1, DLGAP5, HJURP, KIF15, IGHG1, IGHGP, IGHG4, IGLC2, IGLC3"

features <- unlist(strsplit(features_text, ","))
features <- trimws(features)
features <- features[features != ""]

# 3.2 固定参数
feature_split_by <- "cellType"
group_by <- "cellType"
layer_use <- "data"

# 3.3 拆分依据参数
split_by <- "group"
# 可选：
# split_by <- NULL
# split_by <- "cellType"
# split_by <- "group"

# 3.4 表达计算方法参数
exp_method <- "zscore"
# 可选：
# exp_method <- "raw"
# exp_method <- "fc"
# exp_method <- "log2fc"
# exp_method <- "log1p"

# 3.5 细胞注释参数
cell_annotation_text <- "Phase, G2M.Score, CD8B"

cell_annotation <- unlist(strsplit(cell_annotation_text, ","))
cell_annotation <- trimws(cell_annotation)
cell_annotation <- cell_annotation[cell_annotation != ""]

# 3.6 细胞注释调色板参数
cell_annotation_palette_sel <- "Dark2"
cell_annotation_palette <- rep(
  cell_annotation_palette_sel,
  length.out = max(1, length(cell_annotation))
)

# 3.7 分组调色板参数
group_palette <- "Paired"

# 3.8 点和网格参数
add_dot <- TRUE
dotsize <- 10
add_reticle <- FALSE

# 3.9 行名和列名参数
show_row_names <- TRUE
show_column_names <- FALSE
row_names_side <- "right"

# 3.10 热图调色板参数
heatmap_palette <- "OrRd"

# 3.11 保存 PDF 参数
w_pdf <- 10
h_pdf <- 8
name_pdf <- "GroupHeatmap"


# 4. 参数检查模块

# 4.1 检查固定分组列是否存在
if (!(feature_split_by %in% colnames(srt@meta.data))) {
  stop(paste0("seurat@meta.data 中不存在 feature_split_by 列：", feature_split_by))
}

if (!(group_by %in% colnames(srt@meta.data))) {
  stop(paste0("seurat@meta.data 中不存在 group.by 列：", group_by))
}

# 4.2 检查 split.by 是否存在
if (!is.null(split_by) && !(split_by %in% colnames(srt@meta.data))) {
  stop(paste0("seurat@meta.data 中不存在 split.by 列：", split_by))
}

# 4.3 检查 cell_annotation 中的列是否存在
missing_annotation <- setdiff(cell_annotation, colnames(srt@meta.data))

if (length(missing_annotation) > 0) {
  stop(paste0(
    "seurat@meta.data 中不存在以下 cell_annotation 列：",
    paste(missing_annotation, collapse = ", ")
  ))
}


# 5. GroupHeatmap 绘图模块

# 5.1 生成 GroupHeatmap
p_heatmap <- suppressWarnings(
  GroupHeatmap(
    srt = srt,
    features = features,
    feature_split_by = feature_split_by,
    group.by = group_by,
    split.by = split_by,
    layer = layer_use,
    exp_method = exp_method,
    cell_annotation = cell_annotation,
    cell_annotation_palette = cell_annotation_palette,
    add_dot = add_dot,
    dot_size = grid::unit(dotsize, "mm"),
    add_reticle = add_reticle,
    heatmap_palette = heatmap_palette,
    group_palette = group_palette,
    row_names_side = row_names_side,
    show_row_names = show_row_names,
    show_column_names = show_column_names
  )
)

# 5.2 保存 GroupHeatmap
pdf(
  file = file.path(out_dir, paste0(name_pdf, ".pdf")),
  width = w_pdf,
  height = h_pdf
)

suppressWarnings(suppressMessages(print(p_heatmap)))

dev.off()


# 6. RColorBrewer 调色板预览模块

# 6.1 调色板预览保存参数
w_brewer <- 12
h_brewer <- 10
name_brewer <- "RColorBrewer全部调色板预览"

# 6.2 整理 RColorBrewer 调色板数据
brewer_info <- RColorBrewer::brewer.pal.info
brewer_info$pal <- rownames(brewer_info)
brewer_info$category <- factor(
  brewer_info$category,
  levels = c("qual", "seq", "div")
)

brewer_info <- brewer_info[
  order(brewer_info$category, brewer_info$pal),
  ,
  drop = FALSE
]

max_show <- 12

brewer_dat <- do.call(
  rbind,
  lapply(seq_len(nrow(brewer_info)), function(i) {
    pal <- brewer_info$pal[i]
    cat <- as.character(brewer_info$category[i])
    maxn <- brewer_info$maxcolors[i]
    n_use <- min(maxn, max_show)
    cols <- RColorBrewer::brewer.pal(n_use, pal)
    
    data.frame(
      category = cat,
      pal = pal,
      idx = seq_along(cols),
      col = cols,
      stringsAsFactors = FALSE
    )
  })
)

brewer_dat$category_label <- dplyr::recode(
  brewer_dat$category,
  "qual" = "Qualitative（定性）",
  "seq" = "Sequential（顺序）",
  "div" = "Diverging（发散）"
)

# 6.3 生成 RColorBrewer 调色板预览图
p_brewer_all <- ggplot(
  brewer_dat,
  aes(x = idx, y = pal, fill = col)
) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_identity() +
  facet_grid(
    category_label ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  scale_x_continuous(breaks = 1:max_show) +
  theme_minimal(base_size = 12) +
  theme(
    axis.title = element_blank(),
    panel.grid = element_blank(),
    strip.text.y = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 9),
    plot.margin = margin(10, 10, 10, 10)
  )

# 6.4 保存 RColorBrewer 调色板预览图
ggsave(
  filename = file.path(out_dir, paste0(name_brewer, ".pdf")),
  plot = p_brewer_all,
  width = w_brewer,
  height = h_brewer,
  device = "pdf"
)


# 7. 参数记录模块

# 7.1 参数文件保存参数
name_params <- "GroupHeatmap_parameters"

# 7.2 生成参数记录文本
param_text <- paste0(
  "本次分析参数总结：\n",
  "- features：", features_text, "\n",
  "- feature_split_by：", feature_split_by, "\n",
  "- group.by：", group_by, "\n",
  "- split.by：", ifelse(is.null(split_by), "NULL", split_by), "\n",
  "- layer：", layer_use, "\n",
  "- exp_method：", exp_method, "\n",
  "- cell_annotation：", paste(cell_annotation, collapse = ", "), "\n",
  "- cell_annotation_palette：", paste(cell_annotation_palette, collapse = ", "), "\n",
  "- add_dot：", add_dot, "\n",
  "- dotsize(mm)：", dotsize, "\n",
  "- add_reticle：", add_reticle, "\n",
  "- show_row_names：", show_row_names, "\n",
  "- show_column_names：", show_column_names, "\n",
  "- row_names_side：", row_names_side, "\n",
  "- heatmap_palette：", heatmap_palette, "\n",
  "- group_palette：", group_palette, "\n"
)

# 7.3 保存参数记录

writeLines(
  param_text,
  con = file.path(out_dir, paste0(name_params, ".txt"))
)