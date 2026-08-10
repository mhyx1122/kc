library(Seurat)
library(SCpubr)

# 1. 检查全局环境中是否存在 seurat 对象
if (!exists("seurat", envir = .GlobalEnv)) {
  stop("全局环境中没有 seurat 对象")
}

srt <- get("seurat", envir = .GlobalEnv)

# 2. 检查 do_DotPlot 函数是否存在
if (!exists("do_DotPlot")) {
  stop("当前 R 环境中没有 do_DotPlot 函数，请确认 SCpubr 包是否正常加载")
}

do_fun <- get("do_DotPlot")


# 3. 输出目录模块

# 3.1 输出目录参数
out_dir <- "3.1细胞注释"

# 3.2 创建输出目录
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}


# 4. marker CSV 读取模块

# 4.1 marker CSV 文件路径
# CSV 格式要求：
# 每一列是一个细胞类型或分组名称
# 每一列下面填写对应 marker 基因
marker_csv_file <- "marker.csv"

# 4.2 读取 marker CSV
marker_df <- read.csv(
  marker_csv_file,
  header = TRUE,
  check.names = FALSE
)

# 4.3 整理 marker 基因列表
features_list <- list()

for (nm in names(marker_df)) {
  genes <- as.character(marker_df[[nm]])
  genes <- genes[genes != "" & !is.na(genes)]
  features_list[[nm]] <- genes
}

# 4.4 查看整理后的 marker 分组
names(features_list)

# 4.5 查看每组 marker 数量
sapply(features_list, length)


# 5. DotPlot 气泡图模块

# 5.1 DotPlot 绘图参数
dot_scale <- 20
legend_framewidth <- 2
font_size <- 20

# 5.2 DotPlot 保存参数
w_dot <- 30
h_dot <- 15
name_dot <- "1.第一次注释气泡图"

# 5.3 生成 DotPlot
p_dot <- do_fun(
  sample = srt,
  features = features_list,
  dot.scale = dot_scale,
  legend.framewidth = legend_framewidth,
  font.size = font_size
)

# 5.4 显示 DotPlot
print(p_dot)

# 5.5 保存 DotPlot
pdf(
  file = file.path(out_dir, paste0(name_dot, ".pdf")),
  width = w_dot,
  height = h_dot
)

print(p_dot)

dev.off()


# 6. FeaturePlot UMAP 批量绘图模块

# 6.1 UMAP 绘图参数
cols_umap_text <- "grey,red"

cols_umap <- trimws(unlist(strsplit(cols_umap_text, ",")))
cols_umap <- cols_umap[cols_umap != ""]

if (length(cols_umap) < 2) {
  cols_umap <- c("grey", "red")
} else {
  cols_umap <- cols_umap[1:2]
}

pt_size_umap <- 0.1
ncol_umap <- 2
raster_umap <- FALSE

# 6.2 UMAP 保存参数
w_umap <- 12
h_umap <- 12
prefix_umap <- "2.基因分布UMAP"

# 6.3 批量保存每组 marker 的 UMAP FeaturePlot
saved_umap <- 0

for (grp in names(features_list)) {
  
  p_umap <- FeaturePlot(
    srt,
    features = features_list[[grp]],
    reduction = "umap",
    cols = cols_umap,
    min.cutoff = NA,
    max.cutoff = NA,
    ncol = ncol_umap,
    pt.size = pt_size_umap,
    raster = raster_umap
  )
  
  pdf(
    file = file.path(out_dir, paste0(prefix_umap, "_", grp, ".pdf")),
    width = w_umap,
    height = h_umap
  )
  
  print(p_umap)
  
  dev.off()
  
  saved_umap <- saved_umap + 1
}

cat("UMAP FeaturePlot 批量保存完成：", saved_umap, " 个 PDF\n")


# 7. FeaturePlot tSNE 批量绘图模块

# 7.1 tSNE 绘图参数
cols_tsne_text <- "grey,red"

cols_tsne <- trimws(unlist(strsplit(cols_tsne_text, ",")))
cols_tsne <- cols_tsne[cols_tsne != ""]

if (length(cols_tsne) < 2) {
  cols_tsne <- c("grey", "red")
} else {
  cols_tsne <- cols_tsne[1:2]
}

pt_size_tsne <- 0.1
ncol_tsne <- 2
raster_tsne <- FALSE

# 7.2 tSNE 保存参数
w_tsne <- 12
h_tsne <- 12
prefix_tsne <- "2.基因分布TSNE"

# 7.3 批量保存每组 marker 的 tSNE FeaturePlot
saved_tsne <- 0

for (grp in names(features_list)) {
  
  p_tsne <- FeaturePlot(
    srt,
    features = features_list[[grp]],
    reduction = "tsne",
    cols = cols_tsne,
    min.cutoff = NA,
    max.cutoff = NA,
    ncol = ncol_tsne,
    pt.size = pt_size_tsne,
    raster = raster_tsne
  )
  
  pdf(
    file = file.path(out_dir, paste0(prefix_tsne, "_", grp, ".pdf")),
    width = w_tsne,
    height = h_tsne
  )
  
  print(p_tsne)
  
  dev.off()
  
  saved_tsne <- saved_tsne + 1
}

cat("tSNE FeaturePlot 批量保存完成：", saved_tsne, " 个 PDF\n")


# 8. 完成提示
cat("\n3.1 细胞注释绘图完成。\n")
cat("结果保存目录：", out_dir, "\n")
cat("DotPlot 已保存：", file.path(out_dir, paste0(name_dot, ".pdf")), "\n")
cat("UMAP FeaturePlot 保存数量：", saved_umap, "\n")
cat("tSNE FeaturePlot 保存数量：", saved_tsne, "\n")