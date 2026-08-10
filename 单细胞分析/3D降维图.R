library(Seurat)
library(plotly)
library(htmlwidgets)

# 1. 检查全局环境中是否存在 seurat 对象
if (!exists("seurat", envir = .GlobalEnv)) {
  stop("全局环境中没有 seurat 对象")
}

srt <- get("seurat", envir = .GlobalEnv)


# 2. 输出目录模块

# 2.1 输出目录参数
out_dir <- "4.1立体UMAP_TSNE_3D"

# 2.2 创建输出目录
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}


# 3. 3D 降维参数模块

# 3.1 选择 RunTSNE / RunUMAP 使用的 reduction
# 常用可选值："pca"、"harmony"
reduction_method <- "pca"

# 3.2 使用的维度数量
dims_num <- 9

# 3.3 分组信息来源
# 可选值需要存在于 seurat@meta.data 中
# 常用："cellType"、"group"、"orig.ident"
group_type <- "cellType"

# 3.4 是否保留 plotly 3D 网格线和背景盒
keep_grid <- FALSE


# 4. 参数检查模块

# 4.1 检查 reduction 是否存在
if (!(reduction_method %in% Reductions(srt))) {
  stop(paste0("seurat 中不存在 reduction：", reduction_method))
}

# 4.2 检查分组列是否存在
if (!(group_type %in% colnames(srt@meta.data))) {
  stop(paste0("seurat@meta.data 中不存在列：", group_type))
}


# 5. 运行 3D tSNE 和 3D UMAP 模块

# 5.1 运行 3D tSNE
srt <- RunTSNE(
  srt,
  reduction = reduction_method,
  dims = 1:dims_num,
  dim.embed = 3,
  check_duplicates = FALSE,
  verbose = FALSE
)

# 5.2 运行 3D UMAP
srt <- RunUMAP(
  srt,
  reduction = reduction_method,
  dims = 1:dims_num,
  n.components = 3,
  verbose = FALSE
)

# 5.3 提取 3D tSNE 和 3D UMAP 坐标
tsne_embed <- Embeddings(srt, "tsne")[, 1:3, drop = FALSE]
umap_embed <- Embeddings(srt, "umap")[, 1:3, drop = FALSE]

# 5.4 提取分组信息
group_vec <- srt@meta.data[[group_type]]


# 6. tSNE 3D 图模块

# 6.1 tSNE 3D 绘图参数
palette_tsne3d_text <- "#E64B35FF,#4DBBD5FF,#00A087FF,#3C5488FF,#F39B7FFF,#8491B4FF,#91D1C2FF,#7E6148FF,#0073C2FF,#EFC000FF,#868686FF,#CD534CFF,#7AA6DCFF,#003C67FF,#8F7700FF,#3B3B3BFF,#A73030FF,#374E55FF,#DF8F44FF,#00A1D5FF,#B24745FF,#79AF97FF,#6A6599FF,#80796BFF,#BC3C29FF,#0072B5FF,#E18727FF,#20854EFF,#7876B1FF,#6F99ADFF,#FFDC91FF,#EE4C97FF"

palette_tsne3d <- unlist(strsplit(palette_tsne3d_text, ","))
palette_tsne3d <- trimws(palette_tsne3d)
palette_tsne3d <- palette_tsne3d[palette_tsne3d != ""]

marker_size_tsne3d <- 1
name_tsne_html <- "tsne_3d_plot"

# 6.2 设置 tSNE 分组颜色
group_fac_tsne <- as.factor(group_vec)
group_levels_tsne <- levels(group_fac_tsne)

cols_tsne3d <- rep(palette_tsne3d, length.out = length(group_levels_tsne))
names(cols_tsne3d) <- group_levels_tsne

# 6.3 设置 tSNE 3D 坐标轴样式
if (isTRUE(keep_grid)) {
  
  scene_tsne <- list(
    xaxis = list(title = "tSNE 1"),
    yaxis = list(title = "tSNE 2"),
    zaxis = list(title = "tSNE 3")
  )
  
} else {
  
  scene_tsne <- list(
    xaxis = list(title = "tSNE 1", showgrid = FALSE, zeroline = FALSE, showbackground = FALSE),
    yaxis = list(title = "tSNE 2", showgrid = FALSE, zeroline = FALSE, showbackground = FALSE),
    zaxis = list(title = "tSNE 3", showgrid = FALSE, zeroline = FALSE, showbackground = FALSE)
  )
}

# 6.4 生成 tSNE 3D plotly 图
p_tsne3d <- plot_ly(
  x = tsne_embed[, 1],
  y = tsne_embed[, 2],
  z = tsne_embed[, 3],
  type = "scatter3d",
  mode = "markers",
  color = group_fac_tsne,
  colors = cols_tsne3d,
  marker = list(size = marker_size_tsne3d)
)

p_tsne3d <- layout(
  p_tsne3d,
  scene = scene_tsne
)

# 6.5 保存 tSNE 3D HTML
tsne_html_file <- file.path(out_dir, paste0(name_tsne_html, ".html"))

htmlwidgets::saveWidget(
  widget = p_tsne3d,
  file = tsne_html_file,
  selfcontained = TRUE,
  libdir = NULL,
  background = "white",
  title = name_tsne_html
)

tsne_lib_folder <- sub(".html$", "_files", tsne_html_file, fixed = TRUE)

if (dir.exists(tsne_lib_folder)) {
  unlink(tsne_lib_folder, recursive = TRUE, force = TRUE)
}


# 7. UMAP 3D 图模块

# 7.1 UMAP 3D 绘图参数
palette_umap3d_text <- "#E64B35FF,#4DBBD5FF,#00A087FF,#3C5488FF,#F39B7FFF,#8491B4FF,#91D1C2FF,#7E6148FF,#0073C2FF,#EFC000FF,#868686FF,#CD534CFF,#7AA6DCFF,#003C67FF,#8F7700FF,#3B3B3BFF,#A73030FF,#374E55FF,#DF8F44FF,#00A1D5FF,#B24745FF,#79AF97FF,#6A6599FF,#80796BFF,#BC3C29FF,#0072B5FF,#E18727FF,#20854EFF,#7876B1FF,#6F99ADFF,#FFDC91FF,#EE4C97FF"

palette_umap3d <- unlist(strsplit(palette_umap3d_text, ","))
palette_umap3d <- trimws(palette_umap3d)
palette_umap3d <- palette_umap3d[palette_umap3d != ""]

marker_size_umap3d <- 2
name_umap_html <- "umap_3d_plot"

# 7.2 设置 UMAP 分组颜色
group_fac_umap <- as.factor(group_vec)
group_levels_umap <- levels(group_fac_umap)

cols_umap3d <- rep(palette_umap3d, length.out = length(group_levels_umap))
names(cols_umap3d) <- group_levels_umap

# 7.3 设置 UMAP 3D 坐标轴样式
if (isTRUE(keep_grid)) {
  
  scene_umap <- list(
    xaxis = list(title = "UMAP 1"),
    yaxis = list(title = "UMAP 2"),
    zaxis = list(title = "UMAP 3")
  )
  
} else {
  
  scene_umap <- list(
    xaxis = list(title = "UMAP 1", showgrid = FALSE, zeroline = FALSE, showbackground = FALSE),
    yaxis = list(title = "UMAP 2", showgrid = FALSE, zeroline = FALSE, showbackground = FALSE),
    zaxis = list(title = "UMAP 3", showgrid = FALSE, zeroline = FALSE, showbackground = FALSE)
  )
}

# 7.4 生成 UMAP 3D plotly 图
p_umap3d <- plot_ly(
  x = umap_embed[, 1],
  y = umap_embed[, 2],
  z = umap_embed[, 3],
  type = "scatter3d",
  mode = "markers",
  color = group_fac_umap,
  colors = cols_umap3d,
  marker = list(size = marker_size_umap3d)
)

p_umap3d <- layout(
  p_umap3d,
  scene = scene_umap
)

# 7.5 保存 UMAP 3D HTML
umap_html_file <- file.path(out_dir, paste0(name_umap_html, ".html"))

htmlwidgets::saveWidget(
  widget = p_umap3d,
  file = umap_html_file,
  selfcontained = TRUE,
  libdir = NULL,
  background = "white",
  title = name_umap_html
)

umap_lib_folder <- sub(".html$", "_files", umap_html_file, fixed = TRUE)

if (dir.exists(umap_lib_folder)) {
  unlink(umap_lib_folder, recursive = TRUE, force = TRUE)
}


# 8. 参数记录模块

# 8.1 参数文件保存参数
name_params <- "seurat_3d_parameters"

# 8.2 生成参数记录文本
param_text <- paste0(
  "本次分析参数总结：\n",
  "- reduction_method：", reduction_method, "\n",
  "- dims：1 ~ ", dims_num, "\n",
  "- group_type：", group_type, "\n",
  "- keep_grid：", keep_grid, "\n",
  "- tSNE 点大小：", marker_size_tsne3d, "\n",
  "- UMAP 点大小：", marker_size_umap3d, "\n",
  "- tSNE 颜色集合：", palette_tsne3d_text, "\n",
  "- UMAP 颜色集合：", palette_umap3d_text, "\n",
  "\n写作提示词（自行组装语言，或借助AI组装）：\n",
  "1.使用指定降维空间（pca/harmony）的前n个PC作为输入，构建3D tSNE与3D UMAP嵌入。\n",
  "2.按照指定meta.data分组列进行着色展示，交互式三维散点图用于观察群体分离与连续变化。\n",
  "3.可选择隐藏网格线与背景盒，以获得更适合论文排版的简洁风格。\n",
  "4.结果支持导出HTML便于浏览器中旋转/缩放查看。"
)

# 8.3 保存参数记录
writeLines(
  param_text,
  con = file.path(out_dir, paste0(name_params, ".txt"))
)
