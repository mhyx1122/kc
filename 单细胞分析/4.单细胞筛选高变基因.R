library(Seurat)
library(ggplot2)

# 1. 检查全局环境中是否存在 seurat 对象
if (!exists("seurat", envir = .GlobalEnv)) {
  stop("全局环境中未找到对象 seurat。请先生成/加载 seurat。")
}

srt0 <- get("seurat", envir = .GlobalEnv)


# 2. Step2 参数设置模块

# 2.1 输出目录参数
scRNA_default_dir <- "1.2筛选高变基因"

# 2.2 高变基因和归一化参数
nfeatures <- 2000
scale.factor <- 10000

# 2.3 JoinLayers 顺序参数
# 可选值：
# "A"：先 JoinLayers，再 Normalize/HVG/Scale/PCA
# "B"：先 Normalize/HVG/Scale/PCA，最后 JoinLayers，去批次时用这个
scRNA_join_order <- "A"

# 2.4 ScaleData 标准化 features 参数
# 可选值：
# "all"：用全部基因 rownames(seurat)，更慢
# "hvg"：只用高变基因 VariableFeatures，节约内存，更快
scRNA_scale_features_mode <- "hvg"


# 3. 创建输出目录
if (!dir.exists(scRNA_default_dir)) {
  dir.create(scRNA_default_dir, recursive = TRUE)
}


# 4. 按 JoinLayers 顺序运行 Step2

if (scRNA_join_order == "A") {
  
  # 4.1 A 策略：先 JoinLayers，再 NormalizeData / FindVariableFeatures / ScaleData / RunPCA
  srt <- JoinLayers(srt0)
  
  srt <- NormalizeData(
    srt,
    normalization.method = "LogNormalize",
    scale.factor = scale.factor
  )
  
  srt <- FindVariableFeatures(
    srt,
    selection.method = "vst",
    nfeatures = nfeatures
  )
  
  # 4.2 提取高变基因
  hvg <- VariableFeatures(srt)
  
  # 4.3 生成高变基因图
  p_hvg <- suppressWarnings(
    VariableFeaturePlot(srt)
  )
  
  # 4.4 生成标记前10高变基因图
  top10 <- head(hvg, 10)
  
  p_hvg_top10 <- suppressWarnings(
    LabelPoints(
      plot = VariableFeaturePlot(srt),
      points = top10,
      repel = TRUE,
      xnudge = 0,
      ynudge = 0
    )
  )
  
  # 4.5 ScaleData
  if (scRNA_scale_features_mode == "all") {
    
    srt <- ScaleData(
      srt,
      features = rownames(srt)
    )
    
  } else {
    
    srt <- ScaleData(
      srt,
      features = VariableFeatures(srt)
    )
  }
  
  # 4.6 RunPCA
  srt <- RunPCA(
    srt,
    features = VariableFeatures(srt)
  )
  
} else {
  
  # 5.1 B 策略：先 NormalizeData / FindVariableFeatures / ScaleData / RunPCA，最后 JoinLayers
  srt <- srt0
  
  srt <- NormalizeData(
    srt,
    normalization.method = "LogNormalize",
    scale.factor = scale.factor
  )
  
  srt <- FindVariableFeatures(
    srt,
    selection.method = "vst",
    nfeatures = nfeatures
  )
  
  # 5.2 提取高变基因
  hvg <- VariableFeatures(srt)
  
  # 5.3 生成高变基因图
  p_hvg <- suppressWarnings(
    VariableFeaturePlot(srt)
  )
  
  # 5.4 生成标记前10高变基因图
  top10 <- head(hvg, 10)
  
  p_hvg_top10 <- suppressWarnings(
    LabelPoints(
      plot = VariableFeaturePlot(srt),
      points = top10,
      repel = TRUE,
      xnudge = 0,
      ynudge = 0
    )
  )
  
  # 5.5 ScaleData
  if (scRNA_scale_features_mode == "all") {
    
    srt <- ScaleData(
      srt,
      features = rownames(srt)
    )
    
  } else {
    
    srt <- ScaleData(
      srt,
      features = VariableFeatures(srt)
    )
  }
  
  # 5.6 RunPCA
  srt <- RunPCA(
    srt,
    features = VariableFeatures(srt)
  )
  
  # 5.7 最后 JoinLayers
  srt <- JoinLayers(srt)
}


# 6. 将处理后的对象写回 seurat
seurat <- srt


# 7. 查看 HVG 结果

# 7.1 输出 HVG 总数
cat("HVG 总数：", length(hvg), "\n\n")

# 7.2 输出前50个 HVG
cat("前50个 HVG：\n")
print(head(hvg, 50))


# 8. 保存高变基因图

# 8.1 高变基因图保存参数
scRNA_w_hvg <- 7
scRNA_h_hvg <- 7
scRNA_name_hvg <- "1_高变基因"

# 8.2 显示高变基因图
print(p_hvg)

# 8.3 保存高变基因图
ggsave(
  filename = file.path(scRNA_default_dir, paste0(scRNA_name_hvg, ".pdf")),
  plot = p_hvg,
  width = scRNA_w_hvg,
  height = scRNA_h_hvg,
  device = "pdf"
)


# 9. 保存标记前10高变基因图

# 9.1 标记前10高变基因图保存参数
scRNA_w_hvg_top10 <- 7
scRNA_h_hvg_top10 <- 7
scRNA_name_hvg_top10 <- "2_标记前10高变基因"

# 9.2 显示标记前10高变基因图
print(p_hvg_top10)

# 9.3 保存标记前10高变基因图
ggsave(
  filename = file.path(scRNA_default_dir, paste0(scRNA_name_hvg_top10, ".pdf")),
  plot = p_hvg_top10,
  width = scRNA_w_hvg_top10,
  height = scRNA_h_hvg_top10,
  device = "pdf"
)


# 10. 导出 HVG 列表

# 10.1 HVG 文件保存参数
hvg_csv <- file.path(scRNA_default_dir, "all_variable_features.csv")

# 10.2 写出 HVG 表格
write.csv(
  data.frame(highly_variable_gene = hvg),
  hvg_csv,
  row.names = FALSE
)


# 11. 完成提示
cat("\nStep2 完成：已更新当前 R 环境中的 seurat 对象。\n")
cat("Join顺序：", scRNA_join_order, "\n")
cat("ScaleData features：", scRNA_scale_features_mode, "\n")
cat("HVG数量：", length(hvg), "\n")
cat("结果保存目录：", scRNA_default_dir, "\n")