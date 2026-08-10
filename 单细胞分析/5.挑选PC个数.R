library(Seurat)

# 1. 检查全局环境中是否存在 seurat 对象
if (!exists("seurat", envir = .GlobalEnv)) {
  stop("全局环境中没有 seurat 对象")
}

srt <- get("seurat", envir = .GlobalEnv)

if (!"pca" %in% Reductions(srt)) {
  stop("seurat 中没有 reduction pca，请先 RunPCA")
}


# 2. 主成分选择参数模块

# 2.1 JackStraw 重复次数
# 如果设置为 0，则不运行 JackStraw
num_replicate <- 0

# 2.2 PCA 维度范围
dims_end <- 30

# 2.3 PC 热图使用的细胞数
cells_to_use <- 1000


# 3. 输出目录模块

# 3.1 输出目录参数
out_dir <- "2.1确定主成分的个数"

# 3.2 创建输出目录
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}


# 4. JackStraw 分析模块

if (num_replicate > 0) {
  
  # 4.1 运行 JackStraw
  srt_js <- JackStraw(
    srt,
    num.replicate = num_replicate
  )
  
  # 4.2 计算 JackStraw 评分
  srt_js <- ScoreJackStraw(
    srt_js,
    dims = 1:dims_end
  )
  
} else {
  
  # 4.3 不运行 JackStraw
  srt_js <- NULL
}


# 5. PC 热图模块

# 5.1 PC 热图保存参数
w_heatmap <- 10
h_heatmap <- 7
name_heatmap <- "1.PC热图"

# 5.2 设置 PC 热图使用的细胞数
cells_use <- min(cells_to_use, ncol(srt))

# 5.3 显示 PC 热图
print(
  DimHeatmap(
    srt,
    reduction = "pca",
    dims = 1:dims_end,
    cells = cells_use,
    balanced = TRUE
  )
)

# 5.4 保存 PC 热图
pdf(
  file = file.path(out_dir, paste0(name_heatmap, ".pdf")),
  width = w_heatmap,
  height = h_heatmap
)

print(
  DimHeatmap(
    srt,
    reduction = "pca",
    dims = 1:dims_end,
    cells = cells_use,
    balanced = TRUE
  )
)

dev.off()


# 6. JackStrawPlot 模块

# 6.1 JackStrawPlot 保存参数
w_jack <- 10
h_jack <- 7
name_jack <- "2.JackStrawPlot"

# 6.2 如果 num_replicate > 0，则生成并保存 JackStrawPlot
if (num_replicate > 0) {
  
  p_jack <- JackStrawPlot(
    srt_js,
    dims = 1:dims_end
  )
  
  print(p_jack)
  
  pdf(
    file = file.path(out_dir, paste0(name_jack, ".pdf")),
    width = w_jack,
    height = h_jack
  )
  
  print(p_jack)
  
  dev.off()
  
} else {
  
  cat("num_replicate = 0，跳过 JackStrawPlot。\n")
}


# 7. ElbowPlot 模块

# 7.1 ElbowPlot 保存参数
w_elbow <- 10
h_elbow <- 7
name_elbow <- "3.ElbowPlot"

# 7.2 生成 ElbowPlot
p_elbow <- ElbowPlot(
  srt,
  ndims = dims_end,
  reduction = "pca"
)

# 7.3 显示 ElbowPlot
print(p_elbow)

# 7.4 保存 ElbowPlot
pdf(
  file = file.path(out_dir, paste0(name_elbow, ".pdf")),
  width = w_elbow,
  height = h_elbow
)

print(p_elbow)

dev.off()


# 8. 参数记录模块

# 8.1 参数文件保存参数
name_params <- "PC_choose_parameters"

# 8.2 生成参数记录文本
param_text <- paste0(
  "本次分析参数总结：\n",
  "- PCA 维度范围：1 ~ ", dims_end, "\n",
  "- PC 热图细胞数：", cells_to_use, "\n",
  "- JackStraw 重复次数：", num_replicate, "\n",
  "写作提示词（自行组装语言，或借助AI组装）：\n",
  "1.主成分数量的选择基于 PCA 的 ElbowPlot 及 PC 热图结果进行综合判断。\n",
  "2.当 JackStraw 重复次数大于 0 时，进一步结合 JackStraw 分析结果对主成分显著性进行评估；\n",
  "当 JackStraw 重复次数设为 0 时，则不进行 JackStraw 随机化检验。\n",
  "3.在保证主要生物学信号得以保留的前提下，最终选取前X个主成分用于后续的降维可视化和聚类分析。"
)

# 8.3 保存参数记录
writeLines(
  param_text,
  con = file.path(out_dir, paste0(name_params, ".txt"))
)


# 9. 完成提示
cat("\n2.1 主成分选择分析完成。\n")
cat("输出目录：", out_dir, "\n")
cat("PC 热图已保存：", file.path(out_dir, paste0(name_heatmap, ".pdf")), "\n")

if (num_replicate > 0) {
  cat("JackStrawPlot 已保存：", file.path(out_dir, paste0(name_jack, ".pdf")), "\n")
} else {
  cat("JackStrawPlot 未保存，因为 num_replicate = 0。\n")
}

cat("ElbowPlot 已保存：", file.path(out_dir, paste0(name_elbow, ".pdf")), "\n")
cat("参数文件已保存：", file.path(out_dir, paste0(name_params, ".txt")), "\n")