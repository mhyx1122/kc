suppressPackageStartupMessages({
  library(ggplot2)
  library(scRNAtoolVis)
})

# 1. 读取差异基因结果并设置公共参数

# 1.1 参数设置

out_dir <- "差异基因可视化图"

top_gene_n <- 5
log2fc_cutoff <- 0.5

# 可选："updown"、"adjustP"
volcano_col_type <- "updown"

back_col <- "white"

volcano_colors <- c(
  "#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF",
  "#91D1C2FF", "#7E6148FF", "#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF",
  "#7AA6DCFF", "#003C67FF", "#8F7700FF", "#3B3B3BFF", "#A73030FF", "#374E55FF",
  "#DF8F44FF", "#00A1D5FF", "#B24745FF", "#79AF97FF", "#6A6599FF", "#80796BFF",
  "#BC3C29FF", "#0072B5FF", "#E18727FF", "#20854EFF", "#7876B1FF", "#6F99ADFF",
  "#FFDC91FF", "#EE4C97FF"
)

# 1.2 创建输出文件夹

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# 1.3 检查并读取all_markers对象

if (!exists("all_markers", envir = .GlobalEnv)) {
  stop("全局环境中没有 all_markers 对象，请先完成差异表达分析。")
}

deg_markers <- get("all_markers", envir = .GlobalEnv)

if (!is.data.frame(deg_markers)) {
  stop("all_markers 不是数据框，请检查差异分析结果。")
}

if (nrow(deg_markers) == 0) {
  stop("all_markers 中没有差异基因结果。")
}

if (length(volcano_colors) == 0) {
  stop("颜色集合不能为空。")
}

if (!volcano_col_type %in% c("updown", "adjustP")) {
  stop("volcano_col_type 只能设置为 updown 或 adjustP。")
}

# 2. 绘制差异基因火山图

# 2.1 参数设置

point_size1 <- 5
celltype_size1 <- 6
polar1 <- TRUE

w_volcano1 <- 15
h_volcano1 <- 12
name_volcano1 <- "1.差异基因火山图"

# 2.2 绘制火山图

volcano_plot1 <- jjVolcano(
  diffData = deg_markers,
  col.type = volcano_col_type,
  log2FC.cutoff = log2fc_cutoff,
  tile.col = volcano_colors,
  back.col = back_col,
  size = point_size1,
  celltypeSize = celltype_size1,
  topGeneN = top_gene_n,
  polar = polar1
)

print(volcano_plot1)

# 2.3 保存火山图

ggsave(
  filename = file.path(out_dir, paste0(name_volcano1, ".pdf")),
  plot = volcano_plot1,
  width = w_volcano1,
  height = h_volcano1,
  device = "pdf"
)

# 3. 绘制指定差异基因火山图

# 3.1 参数设置

selected_markers <- "FOLR2,TRAC,CCR7,GNLY"

point_size2 <- 3
celltype_size2 <- 4
polar2 <- TRUE

w_volcano2 <- 15
h_volcano2 <- 12
name_volcano2 <- "2.指定差异基因的火山图"

# 3.2 解析指定差异基因

mygenes <- unlist(strsplit(selected_markers, ","))
mygenes <- trimws(mygenes)
mygenes <- mygenes[mygenes != ""]

if (length(mygenes) == 0) {
  stop("指定差异基因不能为空。")
}

# 3.3 绘制指定差异基因火山图

volcano_plot2 <- jjVolcano(
  diffData = deg_markers,
  myMarkers = mygenes,
  log2FC.cutoff = log2fc_cutoff,
  tile.col = volcano_colors,
  col.type = volcano_col_type,
  back.col = back_col,
  size = point_size2,
  celltypeSize = celltype_size2,
  polar = polar2
)

print(volcano_plot2)

# 3.4 保存指定差异基因火山图

ggsave(
  filename = file.path(out_dir, paste0(name_volcano2, ".pdf")),
  plot = volcano_plot2,
  width = w_volcano2,
  height = h_volcano2,
  device = "pdf"
)

# 4. 绘制markerVolcano图

# 4.1 参数设置

w_marker <- 12
h_marker <- 10
name_marker <- "3.markerVolcano"

# 4.2 绘制markerVolcano图

marker_volcano_plot <- markerVolcano(
  markers = deg_markers,
  topn = top_gene_n,
  labelCol = volcano_colors
)

print(marker_volcano_plot)

# 4.3 保存markerVolcano图

ggsave(
  filename = file.path(out_dir, paste0(name_marker, ".pdf")),
  plot = marker_volcano_plot,
  width = w_marker,
  height = h_marker,
  device = "pdf"
)

# 5. 保存参数记录

# 5.1 参数设置

name_params <- "volcano_parameters"

# 5.2 生成参数记录

param_text <- paste0(
  "本次分析参数总结：\n",
  
  "【公共参数】\n",
  "- topGeneN：", top_gene_n, "\n",
  "- log2FC.cutoff：", log2fc_cutoff, "\n",
  "- col.type：", volcano_col_type, "\n",
  "- back.col：", back_col, "\n",
  "- 颜色集合：", paste(volcano_colors, collapse = ", "), "\n\n",
  
  "【图1：差异基因火山图】\n",
  "- size：", point_size1, "\n",
  "- celltypeSize：", celltype_size1, "\n",
  "- polar：", polar1, "\n",
  "- 保存宽高：", w_volcano1, " × ", h_volcano1, " 英寸\n",
  "- 文件名：", name_volcano1, ".pdf\n\n",
  
  "【图2：指定差异基因火山图】\n",
  "- 指定差异基因：", paste(mygenes, collapse = ", "), "\n",
  "- size：", point_size2, "\n",
  "- celltypeSize：", celltype_size2, "\n",
  "- polar：", polar2, "\n",
  "- 保存宽高：", w_volcano2, " × ", h_volcano2, " 英寸\n",
  "- 文件名：", name_volcano2, ".pdf\n\n",
  
  "【图3：markerVolcano】\n",
  "- topn：", top_gene_n, "\n",
  "- 保存宽高：", w_marker, " × ", h_marker, " 英寸\n",
  "- 文件名：", name_marker, ".pdf\n",
  
  "- 输出目录：", out_dir, "\n"
)

# 5.3 保存参数记录

param_file <- file.path(
  out_dir,
  paste0(name_params, ".txt")
)

writeLines(
  param_text,
  con = param_file
)