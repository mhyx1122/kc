suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(scop)
})

# 1. 检查全局环境中是否存在 seurat 对象
if (!exists("seurat", envir = .GlobalEnv)) {
  stop("全局环境中没有 seurat 对象")
}

srt <- get("seurat", envir = .GlobalEnv)

if (!exists("FeatureDimPlot")) {
  stop("当前 R 环境中没有 FeatureDimPlot 函数，请确认 scop 包是否正常加载")
}


# 2. 输出目录模块

# 2.1 输出目录参数
out_dir <- "7.基因表达量可视化"

# 2.2 创建输出目录
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}


# 3. FeatureDimPlot 参数模块

# 3.1 输入特征参数
features_text <- "CNN3,CD8A,MTF1,PRKCZ"

features <- unlist(strsplit(features_text, ","))
features <- trimws(features)
features <- features[features != ""]

if (length(features) == 0) {
  stop("请输入至少 1 个特征")
}

# 3.2 降维方式参数
reduction_by <- "UMAP"
# 可选：
# reduction_by <- "PCA"
# reduction_by <- "tSNE"

# 3.3 主题参数
theme_use <- "theme_blank"
# 可选：
# theme_use <- "theme_classic"

# 3.4 颜色参数
palcolor1 <- "grey"
palcolor2 <- "red"
palcolor <- c(palcolor1, palcolor2)

# 3.5 点参数
pt_size <- 0.5
pt_alpha <- 0.8

# 3.6 保存 PDF 参数
w_pdf <- 8
h_pdf <- 6
name_pdf <- "FeatureDimPlot"


# 4. 绘图模块

# 4.1 生成 FeatureDimPlot 图
p_featuredim <- FeatureDimPlot(
  srt = srt,
  features = features,
  reduction = reduction_by,
  theme_use = theme_use,
  palcolor = palcolor,
  pt.size = pt_size,
  pt.alpha = pt_alpha
)

# 4.2 保存 FeatureDimPlot 图
ggsave(
  filename = file.path(out_dir, paste0(name_pdf, ".pdf")),
  plot = p_featuredim,
  width = w_pdf,
  height = h_pdf,
  device = "pdf"
)


# 5. 参数记录模块

# 5.1 参数文件保存参数
name_params <- "FeatureDimPlot_parameters"

# 5.2 生成参数记录文本
param_text <- paste0(
  "本次分析参数总结：\n",
  "- features：", features_text, "\n",
  "- reduction：", reduction_by, "\n",
  "- theme_use：", theme_use, "\n",
  "- palcolor low：", palcolor1, "\n",
  "- palcolor high：", palcolor2, "\n",
  "- pt.size：", pt_size, "\n",
  "- pt.alpha：", pt_alpha, "\n"
)

# 5.3 保存参数记录
writeLines(
  param_text,
  con = file.path(out_dir, paste0(name_params, ".txt"))
)