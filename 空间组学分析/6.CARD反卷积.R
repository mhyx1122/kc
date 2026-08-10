suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(qs)
  library(CARD)
})

# 1. 构建CARD对象并进行空间反卷积

# 1.1 参数设置

seurat_file <- "单细胞Seurat对象.qs"
out_dir <- "反卷积细胞注释结果"

# 1.2 创建输出文件夹

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# 1.3 读取空间转录组对象

if (!exists("Spatial_Data", envir = .GlobalEnv)) {
  stop("全局环境中不存在 Spatial_Data，请先准备空间对象。")
}

srt_spatial <- get("Spatial_Data", envir = .GlobalEnv)

if (!"cellType" %in% colnames(srt_spatial@meta.data)) {
  stop("Spatial_Data@meta.data 中不存在 cellType 列。")
}

Idents(srt_spatial) <- srt_spatial$cellType

# 1.4 读取单细胞Seurat对象

if (!file.exists(seurat_file)) {
  stop(paste0("未找到单细胞Seurat文件：", seurat_file))
}

sc_obj <- qread(seurat_file)

required_meta_columns <- c("orig.ident", "cellType")
missing_meta_columns <- setdiff(required_meta_columns, colnames(sc_obj@meta.data))

if (length(missing_meta_columns) > 0) {
  stop(
    paste0(
      "单细胞对象的meta.data中缺少以下列：",
      paste(missing_meta_columns, collapse = ", ")
    )
  )
}

# 1.5 提取空间表达矩阵和空间坐标

image_names <- names(srt_spatial@images)

if (length(image_names) == 0) {
  stop("Spatial_Data 中未找到空间图像信息。")
}

image_use <- image_names[1]

spatial_count <- GetAssayData(
  srt_spatial,
  assay = "Spatial",
  layer = "counts"
)

spatial_loca <- GetTissueCoordinates(
  srt_spatial,
  image = image_use
)

spatial_location <- spatial_loca[, 1:2, drop = FALSE]

# 1.6 提取单细胞表达矩阵和元数据

sc_count <- GetAssayData(
  sc_obj,
  assay = "RNA",
  layer = "counts"
)

sc_meta <- sc_obj@meta.data
sc_meta$cellID <- rownames(sc_meta)
sc_meta <- sc_meta[, c("cellID", "orig.ident", "cellType"), drop = FALSE]
sc_meta$CB <- sc_meta$cellID
rownames(sc_meta) <- sc_meta$CB
sc_meta$CB <- NULL

# 1.7 构建CARD对象

CARD_obj <- createCARDObject(
  sc_count = sc_count,
  sc_meta = sc_meta,
  spatial_count = spatial_count,
  spatial_location = spatial_location,
  ct.varname = "cellType",
  ct.select = unique(sc_meta$cellType),
  sample.varname = "orig.ident"
)

# 1.8 运行CARD反卷积

CARD_obj <- suppressWarnings(
  CARD_deconvolution(
    CARD_object = CARD_obj
  )
)

# 1.9 记录本步骤参数

param_step1_text <- paste0(
  "Step 1：CARD反卷积对象构建\n",
  "运行流程说明：\n",
  "1. 读取单细胞Seurat对象和全局Spatial_Data对象。\n",
  "2. 分别提取空间counts、空间坐标、单细胞counts和单细胞meta信息。\n",
  "3. 使用CARD::createCARDObject()构建CARD对象。\n",
  "4. 使用CARD::CARD_deconvolution()进行空间反卷积。\n\n",
  "本次运行参数：\n",
  "- 单细胞文件：", seurat_file, "\n",
  "- 输出目录：", out_dir, "\n",
  "- 空间图像：", image_use, "\n",
  "- 空间spot数量：", nrow(CARD_obj@Proportion_CARD), "\n",
  "- 细胞类型数量：", ncol(CARD_obj@Proportion_CARD), "\n"
)

# 2. 绘制细胞比例饼图

# 2.1 参数设置

Pieradius <- 12

colorsAll <- c(
  "#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF",
  "#91D1C2FF", "#7E6148FF", "#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF",
  "#7AA6DCFF", "#003C67FF", "#8F7700FF", "#3B3B3BFF", "#A73030FF", "#374E55FF",
  "#DF8F44FF", "#00A1D5FF", "#B24745FF", "#79AF97FF", "#6A6599FF", "#80796BFF",
  "#BC3C29FF", "#0072B5FF", "#E18727FF", "#20854EFF", "#7876B1FF", "#6F99ADFF",
  "#FFDC91FF", "#EE4C97FF"
)

w_pie <- 10
h_pie <- 10
name_pie <- "1.细胞比例饼图"

# 2.2 绘制细胞比例饼图

pie_plot <- CARD.visualize.pie(
  proportion = CARD_obj@Proportion_CARD,
  spatial_location = CARD_obj@spatial_location,
  radius = Pieradius,
  colors = colorsAll
)

print(pie_plot)

# 2.3 保存细胞比例饼图

ggsave(
  filename = file.path(out_dir, paste0(name_pie, ".pdf")),
  plot = pie_plot,
  width = w_pie,
  height = h_pie,
  device = "pdf"
)

# 2.4 记录本步骤参数

param_step2_text <- paste0(
  "Step 2：细胞比例饼图\n",
  "运行流程说明：\n",
  "1. 使用CARD.visualize.pie()根据反卷积结果绘制每个空间点的细胞组成饼图。\n",
  "2. 使用radius和colors控制圆半径与颜色。\n\n",
  "本次运行参数：\n",
  "- Pieradius：", Pieradius, "\n",
  "- colors：", paste(colorsAll, collapse = ","), "\n",
  "- 保存宽高：", w_pie, " × ", h_pie, "英寸\n",
  "- 文件名：", name_pie, ".pdf\n"
)

# 3. 绘制特定细胞比例图

# 3.1 参数设置

ct_visualize <- c("TCells", "Endothelial.Cells", "B.Cells")
colors_prop <- c("lightblue", "lightyellow", "red")
NumCols_prop <- 2
pt_size <- 1.6

w_prop <- 10
h_prop <- 10
name_prop <- "2.特定细胞比例图"

# 3.2 检查细胞类型

missing_ct <- setdiff(
  ct_visualize,
  colnames(CARD_obj@Proportion_CARD)
)

if (length(missing_ct) > 0) {
  stop(
    paste0(
      "以下细胞类型不在CARD反卷积结果中：",
      paste(missing_ct, collapse = ", ")
    )
  )
}

# 3.3 绘制特定细胞比例图

prop_plot <- CARD.visualize.prop(
  proportion = CARD_obj@Proportion_CARD,
  spatial_location = CARD_obj@spatial_location,
  ct.visualize = ct_visualize,
  colors = colors_prop,
  NumCols = NumCols_prop,
  pointSize = pt_size
)

print(prop_plot)

# 3.4 保存特定细胞比例图

ggsave(
  filename = file.path(out_dir, paste0(name_prop, ".pdf")),
  plot = prop_plot,
  width = w_prop,
  height = h_prop,
  device = "pdf"
)

# 3.5 记录本步骤参数

param_step3_text <- paste0(
  "Step 3：特定细胞比例图\n",
  "运行流程说明：\n",
  "1. 使用CARD.visualize.prop()对指定细胞类型进行空间比例可视化。\n",
  "2. 使用ct.visualize、colors、NumCols和pointSize控制图形显示。\n\n",
  "本次运行参数：\n",
  "- ct.visualize：", paste(ct_visualize, collapse = ", "), "\n",
  "- colors：", paste(colors_prop, collapse = ","), "\n",
  "- NumCols：", NumCols_prop, "\n",
  "- pointSize：", pt_size, "\n",
  "- 保存宽高：", w_prop, " × ", h_prop, "英寸\n",
  "- 文件名：", name_prop, ".pdf\n"
)

# 4. 绘制特定基因表达图

# 4.1 参数设置

gene_list <- c("CD8A", "TOP2A", "CDK1")
colors_gene <- c("lightblue", "lightyellow", "red")
NumCols_gene <- 2

w_gene <- 10
h_gene <- 10
name_gene <- "3.特定基因表达图"

# 4.2 检查基因

missing_genes <- setdiff(
  gene_list,
  rownames(CARD_obj@spatial_countMat)
)

if (length(missing_genes) > 0) {
  stop(
    paste0(
      "以下基因不在空间表达矩阵中：",
      paste(missing_genes, collapse = ", ")
    )
  )
}

# 4.3 绘制特定基因表达图

gene_plot <- CARD.visualize.gene(
  spatial_expression = CARD_obj@spatial_countMat,
  spatial_location = CARD_obj@spatial_location,
  gene.visualize = gene_list,
  colors = colors_gene,
  NumCols = NumCols_gene
)

print(gene_plot)

# 4.4 保存特定基因表达图

ggsave(
  filename = file.path(out_dir, paste0(name_gene, ".pdf")),
  plot = gene_plot,
  width = w_gene,
  height = h_gene,
  device = "pdf"
)

# 4.5 记录本步骤参数

param_step4_text <- paste0(
  "Step 4：特定基因表达图\n",
  "运行流程说明：\n",
  "1. 使用CARD.visualize.gene()对指定基因在空间中的表达进行展示。\n",
  "2. 使用gene.visualize、colors和NumCols控制图形显示。\n\n",
  "本次运行参数：\n",
  "- gene.visualize：", paste(gene_list, collapse = ", "), "\n",
  "- colors：", paste(colors_gene, collapse = ","), "\n",
  "- NumCols：", NumCols_gene, "\n",
  "- 保存宽高：", w_gene, " × ", h_gene, "英寸\n",
  "- 文件名：", name_gene, ".pdf\n"
)

# 5. 绘制双细胞类型对比图

# 5.1 参数设置

ct2_visualize <- c("TCells", "B.Cells")
colors_ct2_1 <- c("lightblue", "lightyellow", "red")
colors_ct2_2 <- c("lightblue", "lightyellow", "black")

w_ct2 <- 10
h_ct2 <- 10
name_ct2 <- "4.双细胞类型对比图"

# 5.2 检查细胞类型

if (length(ct2_visualize) != 2) {
  stop("双细胞类型对比图需要且只能指定2个细胞类型。")
}

missing_ct2 <- setdiff(
  ct2_visualize,
  colnames(CARD_obj@Proportion_CARD)
)

if (length(missing_ct2) > 0) {
  stop(
    paste0(
      "以下细胞类型不在CARD反卷积结果中：",
      paste(missing_ct2, collapse = ", ")
    )
  )
}

# 5.3 绘制双细胞类型对比图

ct2_plot <- CARD.visualize.prop.2CT(
  proportion = CARD_obj@Proportion_CARD,
  spatial_location = CARD_obj@spatial_location,
  ct2.visualize = ct2_visualize,
  colors = list(colors_ct2_1, colors_ct2_2)
)

print(ct2_plot)

# 5.4 保存双细胞类型对比图

ggsave(
  filename = file.path(out_dir, paste0(name_ct2, ".pdf")),
  plot = ct2_plot,
  width = w_ct2,
  height = h_ct2,
  device = "pdf"
)

# 5.5 记录本步骤参数

param_step5_text <- paste0(
  "Step 5：双细胞类型对比图\n",
  "运行流程说明：\n",
  "1. 使用CARD.visualize.prop.2CT()对两个指定细胞类型进行空间比例对比。\n",
  "2. 使用两组颜色向量分别控制两个细胞类型的色彩映射。\n\n",
  "本次运行参数：\n",
  "- ct2.visualize：", paste(ct2_visualize, collapse = ", "), "\n",
  "- colors_ct2_1：", paste(colors_ct2_1, collapse = ","), "\n",
  "- colors_ct2_2：", paste(colors_ct2_2, collapse = ","), "\n",
  "- 保存宽高：", w_ct2, " × ", h_ct2, "英寸\n",
  "- 文件名：", name_ct2, ".pdf\n"
)

# 6. 绘制各细胞类型间相关性图

# 6.1 参数设置

w_cor <- 10
h_cor <- 10
name_cor <- "5.各细胞间相关性图"

# 6.2 绘制相关性图

cor_plot <- CARD.visualize.Cor(
  CARD_obj@Proportion_CARD,
  colors = NULL
)

print(cor_plot)

# 6.3 保存相关性图

ggsave(
  filename = file.path(out_dir, paste0(name_cor, ".pdf")),
  plot = cor_plot,
  width = w_cor,
  height = h_cor,
  device = "pdf"
)

# 6.4 记录本步骤参数

param_step6_text <- paste0(
  "Step 6：各细胞间相关性图\n",
  "运行流程说明：\n",
  "1. 使用CARD.visualize.Cor()对反卷积得到的细胞比例矩阵进行相关性可视化。\n",
  "2. 结果反映不同细胞类型在空间中的共分布和相关模式。\n\n",
  "本次运行参数：\n",
  "- 保存宽高：", w_cor, " × ", h_cor, "英寸\n",
  "- 文件名：", name_cor, ".pdf\n"
)

# 7. 生成分析概览并保存参数记录

# 7.1 参数设置

name_params <- "CARD_Spatial_parameters"

# 7.2 生成分析概览

summary_text <- paste0(
  "CARD反卷积分析已完成。\n",
  "空间spot数量：", nrow(CARD_obj@Proportion_CARD), "\n",
  "细胞类型数量：", ncol(CARD_obj@Proportion_CARD), "\n",
  "特定细胞比例图：", paste(ct_visualize, collapse = ", "), "\n",
  "特定基因表达图：", paste(gene_list, collapse = ", "), "\n",
  "双细胞类型对比图：", paste(ct2_visualize, collapse = ", "), "\n"
)

# 7.3 合并参数记录

parameter_summary <- paste(
  param_step1_text,
  param_step2_text,
  param_step3_text,
  param_step4_text,
  param_step5_text,
  param_step6_text,
  sep = "\n\n------------------------------\n\n"
)

# 7.4 保存参数记录

writeLines(
  parameter_summary,
  con = file.path(out_dir, paste0(name_params, ".txt"))
)