suppressPackageStartupMessages({
  library(Seurat)
  library(scCustomize)
  library(plot1cell)
})

# 1. 检查全局环境中是否存在 seurat 对象
if (!exists("seurat", envir = .GlobalEnv)) {
  stop("全局环境中没有 seurat 对象")
}

srt <- get("seurat", envir = .GlobalEnv)

if (!"group" %in% colnames(srt@meta.data)) {
  stop("seurat@meta.data 中没有 group 列")
}

if (!"orig.ident" %in% colnames(srt@meta.data)) {
  stop("seurat@meta.data 中没有 orig.ident 列")
}


# 2. 输出目录模块

# 2.1 输出目录参数
out_dir <- "4.4降维图带圆环"

# 2.2 创建输出目录
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}


# 3. Seurat Assay 转换模块

# 3.1 Convert_Assay 参数
convert_to <- "V3"
assay_use <- "RNA"

# 3.2 转换 Assay 到 V3
srt_v3 <- Convert_Assay(
  seurat_object = srt,
  convert_to = convert_to,
  assay = assay_use
)


# 4. 准备 circlize 数据模块

# 4.1 prepare_circlize_data 参数
scale_val <- 0.8

# 4.2 生成 circlize 数据
circ_data <- prepare_circlize_data(
  srt_v3,
  scale = scale_val
)


# 5. 颜色参数模块

# 5.1 cluster 颜色参数
cluster_cols_text <- "#E64B35FF,#4DBBD5FF,#00A087FF,#3C5488FF,#F39B7FFF,#8491B4FF,#91D1C2FF,#7E6148FF,#0073C2FF,#EFC000FF,#868686FF,#CD534CFF,#7AA6DCFF,#003C67FF,#8F7700FF,#3B3B3BFF,#A73030FF,#374E55FF,#DF8F44FF,#00A1D5FF,#B24745FF,#79AF97FF,#6A6599FF,#80796BFF,#BC3C29FF,#0072B5FF,#E18727FF,#20854EFF,#7876B1FF,#6F99ADFF,#FFDC91FF,#EE4C97FF"

cluster_cols <- trimws(unlist(strsplit(cluster_cols_text, ",")))
cluster_cols <- cluster_cols[cluster_cols != ""]

# 5.2 group 颜色参数
group_cols_text <- "#E64B35FF,#4DBBD5FF,#00A087FF,#3C5488FF,#F39B7FFF,#8491B4FF,#91D1C2FF,#7E6148FF,#0073C2FF,#EFC000FF,#868686FF,#CD534CFF,#7AA6DCFF,#003C67FF,#8F7700FF,#3B3B3BFF,#A73030FF,#374E55FF,#DF8F44FF,#00A1D5FF,#B24745FF,#79AF97FF,#6A6599FF,#80796BFF,#BC3C29FF,#0072B5FF,#E18727FF,#20854EFF,#7876B1FF,#6F99ADFF,#FFDC91FF,#EE4C97FF"

group_cols <- trimws(unlist(strsplit(group_cols_text, ",")))
group_cols <- group_cols[group_cols != ""]

# 5.3 orig.ident 颜色参数
orig_cols_text <- "#E64B35FF,#4DBBD5FF,#00A087FF,#3C5488FF,#F39B7FFF,#8491B4FF,#91D1C2FF,#7E6148FF,#0073C2FF,#EFC000FF,#868686FF,#CD534CFF,#7AA6DCFF,#003C67FF,#8F7700FF,#3B3B3BFF,#A73030FF,#374E55FF,#DF8F44FF,#00A1D5FF,#B24745FF,#79AF97FF,#6A6599FF,#80796BFF,#BC3C29FF,#0072B5FF,#E18727FF,#20854EFF,#7876B1FF,#6F99ADFF,#FFDC91FF,#EE4C97FF"

orig_cols <- trimws(unlist(strsplit(orig_cols_text, ",")))
orig_cols <- orig_cols[orig_cols != ""]

if (length(cluster_cols) == 0) {
  stop("cluster 颜色集合为空")
}

if (length(group_cols) == 0) {
  stop("group 颜色集合为空")
}

if (length(orig_cols) == 0) {
  stop("orig.ident 颜色集合为空")
}


# 6. 分组水平模块

# 6.1 获取 cluster 水平
cluster_lv <- levels(Idents(srt_v3))

if (is.null(cluster_lv) || length(cluster_lv) == 0) {
  cluster_lv <- sort(unique(as.character(Idents(srt_v3))))
}

# 6.2 获取 group 水平
group_lv <- names(table(srt_v3$group))

# 6.3 获取 orig.ident 水平
orig_lv <- names(table(srt_v3$orig.ident))

# 6.4 扩展颜色数量
cluster_colors <- rep(cluster_cols, length.out = length(cluster_lv))
group_colors <- rep(group_cols, length.out = length(group_lv))
orig_colors <- rep(orig_cols, length.out = length(orig_lv))


# 7. Circlize 图参数模块

# 7.1 plot_circlize 参数
do_label <- TRUE
pt_size <- 1
bg_color <- "white"
kde2d_n <- 200
repel <- TRUE
label_cex <- 0.6


# 8. 保存 PNG 模块

# 8.1 PNG 保存参数
w_png <- 6
h_png <- 6
res_png <- 300
name_png <- "circlize_plot"

# 8.2 保存 PNG
png(
  filename = file.path(out_dir, paste0(name_png, ".png")),
  width = w_png,
  height = h_png,
  units = "in",
  res = res_png
)

plot_circlize(
  circ_data,
  do.label = do_label,
  pt.size = pt_size,
  col.use = cluster_colors,
  bg.color = bg_color,
  kde2d.n = kde2d_n,
  repel = repel,
  label.cex = label_cex
)

add_track(
  circ_data,
  group = "group",
  colors = group_colors,
  track_num = 2
)

add_track(
  circ_data,
  group = "orig.ident",
  colors = orig_colors,
  track_num = 3
)

dev.off()


# 9. 保存 PDF 模块

# 9.1 PDF 保存参数
w_pdf <- 6
h_pdf <- 6
name_pdf <- "circlize_plot"

# 9.2 保存 PDF
pdf(
  file = file.path(out_dir, paste0(name_pdf, ".pdf")),
  width = w_pdf,
  height = h_pdf
)

plot_circlize(
  circ_data,
  do.label = do_label,
  pt.size = pt_size,
  col.use = cluster_colors,
  bg.color = bg_color,
  kde2d.n = kde2d_n,
  repel = repel,
  label.cex = label_cex
)

add_track(
  circ_data,
  group = "group",
  colors = group_colors,
  track_num = 2
)

add_track(
  circ_data,
  group = "orig.ident",
  colors = orig_colors,
  track_num = 3
)

dev.off()


# 10. 参数记录模块

# 10.1 参数文件保存参数
name_params <- "circlize_parameters"

# 10.2 生成参数记录文本
param_text <- paste0(
  "本次分析参数总结：\n",
  "- Convert_Assay：convert_to = ", convert_to, ", assay = ", assay_use, "\n",
  "- prepare_circlize_data scale：", scale_val, "\n",
  "- do.label：", do_label, "\n",
  "- pt.size：", pt_size, "\n",
  "- bg.color：", bg_color, "\n",
  "- kde2d.n：", kde2d_n, "\n",
  "- repel：", repel, "\n",
  "- label.cex：", label_cex, "\n",
  "- cluster 颜色集合：", cluster_cols_text, "\n",
  "- group 颜色集合：", group_cols_text, "\n",
  "- orig.ident 颜色集合：", orig_cols_text, "\n",
  "- cluster 数量：", length(cluster_lv), "\n",
  "- group 数量：", length(group_lv), "\n",
  "- orig.ident 数量：", length(orig_lv), "\n"
)

# 10.3 保存参数记录
writeLines(
  param_text,
  con = file.path(out_dir, paste0(name_params, ".txt"))
)