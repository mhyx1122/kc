suppressPackageStartupMessages({
  library(ggplot2)
  library(scRNAtoolVis)
})

# 1. 检查全局环境中是否存在 all_markers 对象
if (!exists("all_markers", envir = .GlobalEnv)) {
  stop("全局环境中没有 all_markers 对象")
}

deg_markers <- get("all_markers", envir = .GlobalEnv)


# 2. 输出目录模块

# 2.1 输出目录参数
out_dir <- "7.差异基因可视化图"

# 2.2 创建输出目录
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}


# 3. 公共参数模块

# 3.1 火山图公共参数
top_gene_n <- 5
log2fc_cutoff <- 0.5
volcano_col_type <- "updown"
back_col <- "white"

# 3.2 颜色参数
color_input <- "#E64B35FF,#4DBBD5FF,#00A087FF,#3C5488FF,#F39B7FFF,#8491B4FF,#91D1C2FF,#7E6148FF,#0073C2FF,#EFC000FF,#868686FF,#CD534CFF,#7AA6DCFF,#003C67FF,#8F7700FF,#3B3B3BFF,#A73030FF,#374E55FF,#DF8F44FF,#00A1D5FF,#B24745FF,#79AF97FF,#6A6599FF,#80796BFF,#BC3C29FF,#0072B5FF,#E18727FF,#20854EFF,#7876B1FF,#6F99ADFF,#FFDC91FF,#EE4C97FF"

pal <- unlist(strsplit(color_input, ","))
pal <- trimws(pal)
pal <- pal[pal != ""]

if (length(pal) == 0) {
  stop("颜色集合不能为空")
}


# 4. 差异基因火山图模块

# 4.1 图1参数
point_size1 <- 5
celltype_size1 <- 6
polar1 <- TRUE

# 4.2 图1保存参数
w_volcano1 <- 15
h_volcano1 <- 12
name_volcano1 <- "1.差异基因火山图"

# 4.3 生成图1
p_volcano1 <- jjVolcano(
  diffData = deg_markers,
  col.type = volcano_col_type,
  log2FC.cutoff = log2fc_cutoff,
  tile.col = pal,
  back.col = back_col,
  size = point_size1,
  celltypeSize = celltype_size1,
  topGeneN = top_gene_n,
  polar = polar1
)

# 4.4 保存图1
ggsave(
  filename = file.path(out_dir, paste0(name_volcano1, ".pdf")),
  plot = p_volcano1,
  width = w_volcano1,
  height = h_volcano1,
  device = "pdf"
)


# 5. 指定差异基因火山图模块

# 5.1 指定差异基因参数
selected_markers_text <- "FOLR2,TRAC,CCR7,GNLY"

mygenes <- unlist(strsplit(selected_markers_text, ","))
mygenes <- trimws(mygenes)
mygenes <- mygenes[mygenes != ""]

if (length(mygenes) == 0) {
  stop("指定差异基因不能为空")
}

# 5.2 图2参数
point_size2 <- 3
celltype_size2 <- 4
polar2 <- TRUE

# 5.3 图2保存参数
w_volcano2 <- 15
h_volcano2 <- 12
name_volcano2 <- "2.指定差异基因的火山图"

# 5.4 生成图2
p_volcano2 <- jjVolcano(
  diffData = deg_markers,
  myMarkers = mygenes,
  log2FC.cutoff = log2fc_cutoff,
  tile.col = pal,
  col.type = volcano_col_type,
  back.col = back_col,
  size = point_size2,
  celltypeSize = celltype_size2,
  polar = polar2
)

# 5.5 保存图2
ggsave(
  filename = file.path(out_dir, paste0(name_volcano2, ".pdf")),
  plot = p_volcano2,
  width = w_volcano2,
  height = h_volcano2,
  device = "pdf"
)


# 6. markerVolcano 模块

# 6.1 markerVolcano 保存参数
w_marker <- 12
h_marker <- 10
name_marker <- "3.markerVolcano"

# 6.2 生成 markerVolcano 图
p_marker_volcano <- markerVolcano(
  markers = deg_markers,
  topn = top_gene_n,
  labelCol = pal
)

# 6.3 保存 markerVolcano 图
ggsave(
  filename = file.path(out_dir, paste0(name_marker, ".pdf")),
  plot = p_marker_volcano,
  width = w_marker,
  height = h_marker,
  device = "pdf"
)


# 7. 参数记录模块

# 7.1 参数文件保存参数
name_params <- "volcano_parameters"

# 7.2 生成参数记录文本
param_text <- paste0(
  "本次分析参数总结：\n",
  "【公共参数】\n",
  "- topGeneN：", top_gene_n, "\n",
  "- log2FC.cutoff：", log2fc_cutoff, "\n",
  "- col.type：", volcano_col_type, "\n",
  "- back.col：", back_col, "\n",
  "- 颜色集合：", paste(pal, collapse = ", "), "\n",
  "【图1参数】\n",
  "- size：", point_size1, "\n",
  "- celltypeSize：", celltype_size1, "\n",
  "- polar：", polar1, "\n",
  "【图2参数】\n",
  "- 指定差异基因：", paste(mygenes, collapse = ", "), "\n",
  "- size：", point_size2, "\n",
  "- celltypeSize：", celltype_size2, "\n",
  "- polar：", polar2, "\n",
  "【markerVolcano参数】\n",
  "- topn：", top_gene_n, "\n",
  "- 输出目录：", out_dir, "\n"
)

# 7.3 保存参数记录
writeLines(
  param_text,
  con = file.path(out_dir, paste0(name_params, ".txt"))
)