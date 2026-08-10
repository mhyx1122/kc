# 1. 加载必要的 R 包

suppressPackageStartupMessages({
  library(Seurat)
})


# 2. 检查 Seurat 对象

if (!exists("seurat", envir = .GlobalEnv)) stop("全局环境中没有找到 seurat 对象。")
if (!inherits(seurat, "Seurat")) stop("全局环境中的 seurat 不是 Seurat 对象。")
if (!"RNA" %in% names(seurat@assays)) stop("seurat 对象中不存在 RNA assay。")

DefaultAssay(seurat) <- "RNA"


# 3. 设置分析参数

gene_name <- "CD274"
plot_colors <- c("grey", "red")
plot_width <- 6.5
plot_height <- 5
output_dir <- "单基因分组"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (!(gene_name %in% rownames(seurat[["RNA"]]))) stop(paste0("目标基因 ", gene_name, " 不存在于 RNA assay 中。"))


# 4. 获取目标基因表达量并按照中位数分组

gene_expr <- FetchData(seurat, vars = gene_name, assay = "RNA", layer = "data")[, 1]
median_expr <- median(gene_expr)

seurat$cellType <- ifelse(gene_expr > median_expr, paste0(gene_name, "_High"), paste0(gene_name, "_Low"))
seurat$cellType <- factor(seurat$cellType, levels = c(paste0(gene_name, "_Low"), paste0(gene_name, "_High")))

print(table(seurat$cellType))

Idents(seurat) <- seurat$cellType


# 5. 绘制并保存 UMAP 分组图

if (!"umap" %in% names(seurat@reductions)) stop("seurat 对象中不存在 umap 降维结果。")

pdf(file.path(output_dir, "单基因中位数对所有细胞分组(umap).pdf"), width = plot_width, height = plot_height)
print(DimPlot(seurat, group.by = "cellType", reduction = "umap", cols = plot_colors))
dev.off()


# 6. 绘制并保存 tSNE 分组图

if (!"tsne" %in% names(seurat@reductions)) stop("seurat 对象中不存在 tsne 降维结果。")

pdf(file.path(output_dir, "单基因中位数对所有细胞分组(tsne).pdf"), width = plot_width, height = plot_height)
print(DimPlot(seurat, group.by = "cellType", reduction = "tsne", cols = plot_colors))
dev.off()