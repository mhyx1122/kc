library(Seurat)
library(ggplot2)

# 1. 检查 Seurat 对象是否存在
if (!exists("seurat", envir = .GlobalEnv)) {
  stop("全局环境中未找到对象 seurat。请先在 R 会话中加载 seurat。")
}

seurat_obj <- get("seurat", envir = .GlobalEnv)


# 2. Step1 数据质控过滤模块

# 2.1 质控过滤参数
scRNA_rb_pattern <- "^RPS|^RPL"
scRNA_mt_pattern <- "^MT-"

scRNA_nCount_RNA_min <- 1000
scRNA_nCount_RNA_max <- 50000

scRNA_nFeature_RNA_min <- 200
scRNA_nFeature_RNA_max <- 6000

scRNA_percent_mt_max <- 15

# 2.2 统计质控前每个样本的细胞数量
scRNA_counts_before <- table(seurat_obj$orig.ident)

cat("质控前每个样本细胞数量：\n")
print(scRNA_counts_before)

# 2.3 计算核糖体基因比例和线粒体基因比例
seurat_obj[["percent.rb"]] <- PercentageFeatureSet(
  seurat_obj,
  pattern = scRNA_rb_pattern
)

seurat_obj[["percent.mt"]] <- PercentageFeatureSet(
  seurat_obj,
  pattern = scRNA_mt_pattern
)

# 2.4 保存质控前对象
seurat_before <- seurat_obj

# 2.5 根据质控阈值筛选细胞
keep_count <- seurat_obj$nCount_RNA >= scRNA_nCount_RNA_min &
  seurat_obj$nCount_RNA <= scRNA_nCount_RNA_max

keep_feature <- seurat_obj$nFeature_RNA >= scRNA_nFeature_RNA_min &
  seurat_obj$nFeature_RNA <= scRNA_nFeature_RNA_max

keep_mt <- seurat_obj$percent.mt <= scRNA_percent_mt_max

keep_rb <- seurat_obj$percent.rb <= 100

keep_cells <- keep_count & keep_feature & keep_mt & keep_rb

# 2.6 生成质控过滤后的 Seurat 对象
seurat_filtered <- seurat_obj[, keep_cells]

# 2.7 统计质控后每个样本的细胞数量
scRNA_counts_after <- table(seurat_filtered$orig.ident)

cat("\n质控后每个样本细胞数量：\n")
print(scRNA_counts_after)

# 2.8 将过滤后的对象写回 seurat
seurat <- seurat_filtered


# 3. 输出目录模块

# 3.1 输出目录参数
scRNA_default_dir <- "1.1数据质控"

# 3.2 创建输出目录
if (!dir.exists(scRNA_default_dir)) {
  dir.create(scRNA_default_dir, recursive = TRUE)
}


# 4. 质控前小提琴图模块

# 4.1 质控前小提琴图参数
scRNA_palette <- "#E64B35FF,#4DBBD5FF,#00A087FF,#3C5488FF,#F39B7FFF,#8491B4FF,#91D1C2FF,#7E6148FF,#0073C2FF,#EFC000FF,#868686FF,#CD534CFF"

scRNA_vln_ptsize_before <- 0.1

scRNA_w_qc_before <- 8
scRNA_h_qc_before <- 6
scRNA_name_qc_before <- "1_数据质控前"

pal <- unlist(strsplit(scRNA_palette, ","))

# 4.2 生成质控前小提琴图
p_qc_before <- VlnPlot(
  seurat_before,
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb"),
  layer = "counts",
  ncol = 4,
  cols = pal,
  group.by = "orig.ident",
  pt.size = scRNA_vln_ptsize_before,
  raster = FALSE
)

print(p_qc_before)

# 4.3 保存质控前小提琴图
ggsave(
  filename = file.path(scRNA_default_dir, paste0(scRNA_name_qc_before, ".pdf")),
  plot = p_qc_before,
  width = scRNA_w_qc_before,
  height = scRNA_h_qc_before,
  device = "pdf"
)


# 5. 质控后小提琴图模块

# 5.1 质控后小提琴图参数
scRNA_vln_ptsize_after <- 0.1

scRNA_w_qc_after <- 8
scRNA_h_qc_after <- 6
scRNA_name_qc_after <- "2_数据质控后"

# 5.2 生成质控后小提琴图
p_qc_after <- VlnPlot(
  seurat,
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb"),
  layer = "counts",
  ncol = 4,
  cols = pal,
  group.by = "orig.ident",
  pt.size = scRNA_vln_ptsize_after,
  raster = FALSE
)

print(p_qc_after)

# 5.3 保存质控后小提琴图
ggsave(
  filename = file.path(scRNA_default_dir, paste0(scRNA_name_qc_after, ".pdf")),
  plot = p_qc_after,
  width = scRNA_w_qc_after,
  height = scRNA_h_qc_after,
  device = "pdf"
)


# 6. nCount_RNA vs nFeature_RNA 散点图模块

# 6.1 散点图参数
scRNA_w_scatter1 <- 7
scRNA_h_scatter1 <- 7
scRNA_name_scatter1 <- "3_nCount_vs_nFeature"

# 6.2 生成散点图
p_scatter1 <- FeatureScatter(
  seurat,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA",
  raster = FALSE
) + RotatedAxis()

print(p_scatter1)

# 6.3 保存散点图
ggsave(
  filename = file.path(scRNA_default_dir, paste0(scRNA_name_scatter1, ".pdf")),
  plot = p_scatter1,
  width = scRNA_w_scatter1,
  height = scRNA_h_scatter1,
  device = "pdf"
)


# 7. nCount_RNA vs percent.rb 散点图模块

# 7.1 散点图参数
scRNA_w_scatter2 <- 7
scRNA_h_scatter2 <- 7
scRNA_name_scatter2 <- "4_nCount_vs_percent_rb"

# 7.2 生成散点图
p_scatter2 <- FeatureScatter(
  seurat,
  feature1 = "nCount_RNA",
  feature2 = "percent.rb",
  raster = FALSE
) + RotatedAxis()

print(p_scatter2)

# 7.3 保存散点图
ggsave(
  filename = file.path(scRNA_default_dir, paste0(scRNA_name_scatter2, ".pdf")),
  plot = p_scatter2,
  width = scRNA_w_scatter2,
  height = scRNA_h_scatter2,
  device = "pdf"
)


# 8. nCount_RNA vs percent.mt 散点图模块

# 8.1 散点图参数
scRNA_w_scatter3 <- 7
scRNA_h_scatter3 <- 7
scRNA_name_scatter3 <- "5_nCount_vs_percent_mt"

# 8.2 生成散点图
p_scatter3 <- FeatureScatter(
  seurat,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt",
  raster = FALSE
) + RotatedAxis()

print(p_scatter3)

# 8.3 保存散点图
ggsave(
  filename = file.path(scRNA_default_dir, paste0(scRNA_name_scatter3, ".pdf")),
  plot = p_scatter3,
  width = scRNA_w_scatter3,
  height = scRNA_h_scatter3,
  device = "pdf"
)


# 9. 完成提示
cat("\nStep1 完成：seurat 已过滤并写回当前 R 环境。\n")
cat("全部图已保存到目录：", scRNA_default_dir, "\n")