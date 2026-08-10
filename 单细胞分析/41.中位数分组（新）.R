# 1. 加载必要的 R 包

suppressPackageStartupMessages({
  library(Seurat)
})


# 2. 检查 Seurat 对象

if (!exists("seurat", envir = .GlobalEnv)) stop("全局环境中没有找到 seurat 对象。")
if (!inherits(seurat, "Seurat")) stop("全局环境中的 seurat 不是 Seurat 对象。")
if (!"cellType" %in% colnames(seurat@meta.data)) stop("seurat@meta.data 中不存在 cellType 列。")
if (!"RNA" %in% names(seurat@assays)) stop("seurat 对象中不存在 RNA assay。")
if (!"data" %in% Layers(seurat[["RNA"]])) stop("RNA assay 中不存在 data 层，请先运行 NormalizeData。")


# 3. 设置分析参数

gene_name <- "CD274"
cell_name <- "TCells"

high_group_name <- "T_cell_gene_positive"
low_group_name <- "T_cell_gene_negative"

plot_colors <- c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#7E6148FF", "#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF", "#7AA6DCFF", "#003C67FF", "#8F7700FF", "#3B3B3BFF", "#A73030FF", "#374E55FF", "#DF8F44FF", "#00A1D5FF", "#B24745FF", "#79AF97FF", "#6A6599FF", "#80796BFF", "#BC3C29FF", "#0072B5FF", "#E18727FF", "#20854EFF", "#7876B1FF", "#6F99ADFF", "#FFDC91FF", "#EE4C97FF")

plot_width <- 8
plot_height <- 6
output_dir <- "特定细胞群的单基因中位数分组"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (!(gene_name %in% rownames(seurat[["RNA"]]))) stop(paste0("目标基因 ", gene_name, " 不存在于 RNA assay 中。"))
if (!(cell_name %in% unique(as.character(seurat$cellType)))) stop(paste0("cellType 列中不存在指定细胞类型：", cell_name))
if (high_group_name == low_group_name) stop("高表达组名称和低表达组名称不能相同。")

other_celltypes <- setdiff(unique(as.character(seurat$cellType)), cell_name)

if (high_group_name %in% other_celltypes) stop(paste0("高表达组名称 ", high_group_name, " 与原有细胞类型名称重复。"))
if (low_group_name %in% other_celltypes) stop(paste0("低表达组名称 ", low_group_name, " 与原有细胞类型名称重复。"))


# 4. 提取指定细胞群的细胞条形码

original_celltype <- as.character(seurat$cellType)
target_cells <- rownames(seurat@meta.data)[original_celltype == cell_name]

if (length(target_cells) == 0) stop(paste0("没有筛选到 cellType 为 ", cell_name, " 的细胞。"))


# 5. 提取目标基因表达量并按照目标细胞群内部中位数分组

gene_expr_data <- FetchData(seurat, vars = gene_name, cells = target_cells, assay = "RNA", layer = "data")
gene_expr <- gene_expr_data[, 1]
median_expr <- median(gene_expr, na.rm = TRUE)

gene_group <- ifelse(gene_expr > median_expr, high_group_name, low_group_name)

cat("目标细胞群：", cell_name, "\n")
cat("目标基因：", gene_name, "\n")
cat("目标细胞群内部表达量中位数：", median_expr, "\n")

print(table(gene_group))

if (!any(gene_group == high_group_name)) warning("高表达组中没有细胞，可能是目标基因在该细胞群中的表达量全部相同。")
if (!any(gene_group == low_group_name)) warning("低表达组中没有细胞，可能是目标基因在该细胞群中的表达量分布异常。")


# 6. 更新指定细胞群的 cellType

seurat@meta.data$cellType <- original_celltype
target_positions <- match(names(gene_expr), rownames(seurat@meta.data))
seurat@meta.data$cellType[target_positions] <- gene_group

original_levels <- unique(original_celltype)
new_levels <- unlist(lapply(original_levels, function(x) if (x == cell_name) c(low_group_name, high_group_name) else x))

seurat@meta.data$cellType <- factor(seurat@meta.data$cellType, levels = unique(new_levels))

Idents(seurat) <- "cellType"

print(table(seurat$cellType))


# 7. 设置绘图颜色

celltype_levels <- levels(seurat$cellType)

if (length(plot_colors) < length(celltype_levels)) stop(paste0("颜色数量不足。当前 cellType 有 ", length(celltype_levels), " 组，但只提供了 ", length(plot_colors), " 个颜色。"))

plot_colors_use <- plot_colors[seq_along(celltype_levels)]
names(plot_colors_use) <- celltype_levels


# 8. 绘制并保存 UMAP 分组图

if (!"umap" %in% names(seurat@reductions)) stop("seurat 对象中不存在 umap 降维结果。")

umap_file <- file.path(output_dir, paste0(cell_name, "_", gene_name, "_中位数分组(umap).pdf"))

pdf(umap_file, width = plot_width, height = plot_height)
print(DimPlot(seurat, group.by = "cellType", reduction = "umap", cols = plot_colors_use))
dev.off()


# 9. 绘制并保存 tSNE 分组图

if (!"tsne" %in% names(seurat@reductions)) stop("seurat 对象中不存在 tsne 降维结果。")

tsne_file <- file.path(output_dir, paste0(cell_name, "_", gene_name, "_中位数分组(tsne).pdf"))

pdf(tsne_file, width = plot_width, height = plot_height)
print(DimPlot(seurat, group.by = "cellType", reduction = "tsne", cols = plot_colors_use))
dev.off()