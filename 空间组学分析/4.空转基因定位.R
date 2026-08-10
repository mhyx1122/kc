suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(grid)
})

# 1. 读取对象并检查基因

# 1.1 参数设置

out_dir <- "特定基因的展示"

genes_plot <- "AIF1,C1QA,C1QC"

# 1.2 创建输出文件夹

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# 1.3 检查并读取Spatial_Data对象

if (!exists("Spatial_Data", envir = .GlobalEnv)) {
  stop("全局环境中不存在 Spatial_Data，请先完成前面的分析步骤。")
}

srt <- get("Spatial_Data", envir = .GlobalEnv)

if (!inherits(srt, "Seurat")) {
  stop("Spatial_Data 不是有效的Seurat对象。")
}

# 1.4 解析需要绘制的基因

genes_use <- unlist(strsplit(genes_plot, ","))
genes_use <- trimws(genes_use)
genes_use <- genes_use[genes_use != ""]
genes_use <- unique(genes_use)

if (length(genes_use) == 0) {
  stop("请输入至少一个基因。")
}

# 1.5 检查基因是否存在

missing_genes <- genes_use[!genes_use %in% rownames(srt)]

if (length(missing_genes) > 0) {
  stop(
    paste0(
      "以下基因在Spatial_Data中不存在：",
      paste(missing_genes, collapse = ", ")
    )
  )
}

# 1.6 记录基因载入参数

param_step1_text <- paste0(
  "Step 1：基因列表载入与检查\n",
  "运行流程说明：\n",
  "1. 读取基因字符串，并按英文逗号拆分为基因向量。\n",
  "2. 从全局环境读取Spatial_Data对象。\n",
  "3. 检查每个基因是否存在于Spatial_Data对象中。\n\n",
  "本次运行参数：\n",
  "- 基因列表：", paste(genes_use, collapse = ", "), "\n",
  "- 基因数量：", length(genes_use), "\n"
)

# 2. 绘制UMAP基因分布图

# 2.1 参数设置

selected_gene_umap <- genes_use[1]

w_umap_gene <- 12
h_umap_gene <- 12

pt_size <- 0.5
umap_min_cutoff <- 0
umap_max_cutoff <- NA

umap_color_low <- "grey"
umap_color_high <- "#FF0000"

name_umap_gene <- "UMAP_FeaturePlot"
name_all_umap_prefix <- "UMAP_FeaturePlot"

# 2.2 检查UMAP结果和当前基因

if (!"umap" %in% Reductions(srt)) {
  stop("Spatial_Data中不存在UMAP降维结果。")
}

if (!selected_gene_umap %in% genes_use) {
  stop(
    paste0(
      "selected_gene_umap不在基因列表中，可选基因为：",
      paste(genes_use, collapse = ", ")
    )
  )
}

# 2.3 绘制当前基因的UMAP分布图

umap_gene_plot <- FeaturePlot(
  srt,
  features = selected_gene_umap,
  reduction = "umap",
  cols = c(umap_color_low, umap_color_high),
  min.cutoff = umap_min_cutoff,
  max.cutoff = umap_max_cutoff,
  pt.size = pt_size
)

print(umap_gene_plot)

# 2.4 保存当前基因的UMAP分布图

ggsave(
  filename = file.path(
    out_dir,
    paste0(name_umap_gene, "_", selected_gene_umap, ".pdf")
  ),
  plot = umap_gene_plot,
  width = w_umap_gene,
  height = h_umap_gene,
  device = "pdf"
)

# 2.5 批量导出全部基因的UMAP分布图

for (gene in genes_use) {
  current_umap_plot <- FeaturePlot(
    srt,
    features = gene,
    reduction = "umap",
    cols = c(umap_color_low, umap_color_high),
    min.cutoff = umap_min_cutoff,
    max.cutoff = umap_max_cutoff,
    pt.size = pt_size
  )
  
  ggsave(
    filename = file.path(
      out_dir,
      paste0(name_all_umap_prefix, "_", gene, ".pdf")
    ),
    plot = current_umap_plot,
    width = w_umap_gene,
    height = h_umap_gene,
    device = "pdf"
  )
}

# 2.6 记录UMAP基因分布图参数

param_step2_text <- paste0(
  "Step 2：UMAP基因分布图\n",
  "运行流程说明：\n",
  "1. 使用Seurat::FeaturePlot(reduction = umap)绘制指定基因的UMAP分布图。\n",
  "2. 使用low、high颜色、min.cutoff、max.cutoff和pt.size控制图形样式。\n",
  "3. 保存当前选中基因图，并批量导出全部基因的UMAP图。\n\n",
  "本次运行参数：\n",
  "- 当前基因：", selected_gene_umap, "\n",
  "- low：", umap_color_low, "\n",
  "- high：", umap_color_high, "\n",
  "- min.cutoff：", umap_min_cutoff, "\n",
  "- max.cutoff：", ifelse(is.na(umap_max_cutoff), "NA", umap_max_cutoff), "\n",
  "- pt.size：", pt_size, "\n",
  "- 保存宽高：", w_umap_gene, " × ", h_umap_gene, "英寸\n",
  "- 当前图文件名前缀：", name_umap_gene, "\n",
  "- 批量导出文件名前缀：", name_all_umap_prefix, "\n"
)

# 3. 绘制原版Spatial基因分布图

# 3.1 参数设置

selected_gene_spatial <- genes_use[1]

w_spatial_gene <- 12
h_spatial_gene <- 12

spatial_pt_size_factor <- 1.6
image_alpha <- 1
image_scale <- "lowres"

spatial_min_cutoff <- 0
spatial_max_cutoff <- NA

name_spatial_gene <- "Spatial_FeaturePlot"
name_all_spatial_prefix <- "Spatial_FeaturePlot"

# 3.2 检查空间图像和当前基因

if (length(Images(srt)) == 0) {
  stop("Spatial_Data中不存在空间图像信息，无法绘制SpatialFeaturePlot。")
}

if (!image_scale %in% c("lowres", "hires")) {
  stop("image_scale只能设置为lowres或hires。")
}

if (!selected_gene_spatial %in% genes_use) {
  stop(
    paste0(
      "selected_gene_spatial不在基因列表中，可选基因为：",
      paste(genes_use, collapse = ", ")
    )
  )
}

# 3.3 绘制当前基因的原版Spatial分布图

spatial_gene_plot <- SpatialFeaturePlot(
  srt,
  features = selected_gene_spatial,
  min.cutoff = spatial_min_cutoff,
  max.cutoff = spatial_max_cutoff,
  pt.size.factor = spatial_pt_size_factor,
  image.alpha = image_alpha,
  image.scale = image_scale
)

print(spatial_gene_plot)

# 3.4 保存当前基因的原版Spatial分布图

ggsave(
  filename = file.path(
    out_dir,
    paste0(name_spatial_gene, "_", selected_gene_spatial, ".pdf")
  ),
  plot = spatial_gene_plot,
  width = w_spatial_gene,
  height = h_spatial_gene,
  device = "pdf"
)

# 3.5 批量导出全部基因的原版Spatial分布图

for (gene in genes_use) {
  current_spatial_plot <- SpatialFeaturePlot(
    srt,
    features = gene,
    min.cutoff = spatial_min_cutoff,
    max.cutoff = spatial_max_cutoff,
    pt.size.factor = spatial_pt_size_factor,
    image.alpha = image_alpha,
    image.scale = image_scale
  )
  
  ggsave(
    filename = file.path(
      out_dir,
      paste0(name_all_spatial_prefix, "_", gene, ".pdf")
    ),
    plot = current_spatial_plot,
    width = w_spatial_gene,
    height = h_spatial_gene,
    device = "pdf"
  )
}

# 3.6 记录原版Spatial基因分布图参数

param_step3_text <- paste0(
  "Step 3：Spatial基因分布图（原版）\n",
  "运行流程说明：\n",
  "1. 使用Seurat::SpatialFeaturePlot()绘制指定基因的空间表达分布图。\n",
  "2. 使用min.cutoff、max.cutoff、pt.size.factor、image.alpha和image.scale控制图形样式。\n",
  "3. 保存当前选中基因图，并批量导出全部基因的Spatial图。\n\n",
  "本次运行参数：\n",
  "- 当前基因：", selected_gene_spatial, "\n",
  "- min.cutoff：", spatial_min_cutoff, "\n",
  "- max.cutoff：", ifelse(is.na(spatial_max_cutoff), "NA", spatial_max_cutoff), "\n",
  "- pt.size.factor：", spatial_pt_size_factor, "\n",
  "- image.alpha：", image_alpha, "\n",
  "- image.scale：", image_scale, "\n",
  "- 保存宽高：", w_spatial_gene, " × ", h_spatial_gene, "英寸\n",
  "- 当前图文件名前缀：", name_spatial_gene, "\n",
  "- 批量导出文件名前缀：", name_all_spatial_prefix, "\n"
)

# 4. 绘制自定义颜色Spatial基因分布图

# 4.1 参数设置

selected_gene_spatial_custom <- genes_use[1]

w_spatial_gene_custom <- 12
h_spatial_gene_custom <- 12

spatial_pt_size_factor_custom <- 1.6
image_alpha_custom <- 1
image_scale_custom <- "lowres"

spatial_min_cutoff_custom <- 0
spatial_max_cutoff_custom <- NA

spatial_color_low <- "grey"
spatial_color_high <- "#FF0000"

legend_title_size_custom <- 12
legend_text_size_custom <- 10
legend_key_size_custom <- 0.6

name_spatial_gene_custom <- "Spatial_FeaturePlot_CustomColor"
name_all_spatial_custom_prefix <- "Spatial_FeaturePlot_CustomColor"

# 4.2 检查参数和当前基因

if (!image_scale_custom %in% c("lowres", "hires")) {
  stop("image_scale_custom只能设置为lowres或hires。")
}

if (!selected_gene_spatial_custom %in% genes_use) {
  stop(
    paste0(
      "selected_gene_spatial_custom不在基因列表中，可选基因为：",
      paste(genes_use, collapse = ", ")
    )
  )
}

# 4.3 绘制当前基因的自定义颜色Spatial分布图

spatial_gene_custom_plot <- SpatialFeaturePlot(
  srt,
  features = selected_gene_spatial_custom,
  min.cutoff = spatial_min_cutoff_custom,
  max.cutoff = spatial_max_cutoff_custom,
  pt.size.factor = spatial_pt_size_factor_custom,
  image.alpha = image_alpha_custom,
  image.scale = image_scale_custom
) +
  scale_fill_gradient(
    low = spatial_color_low,
    high = spatial_color_high
  ) +
  theme(
    legend.title = element_text(size = legend_title_size_custom),
    legend.text = element_text(size = legend_text_size_custom),
    legend.key.size = unit(legend_key_size_custom, "cm")
  )

print(spatial_gene_custom_plot)

# 4.4 保存当前基因的自定义颜色Spatial分布图

ggsave(
  filename = file.path(
    out_dir,
    paste0(
      name_spatial_gene_custom,
      "_",
      selected_gene_spatial_custom,
      ".pdf"
    )
  ),
  plot = spatial_gene_custom_plot,
  width = w_spatial_gene_custom,
  height = h_spatial_gene_custom,
  device = "pdf"
)

# 4.5 批量导出全部基因的自定义颜色Spatial分布图

for (gene in genes_use) {
  current_spatial_custom_plot <- SpatialFeaturePlot(
    srt,
    features = gene,
    min.cutoff = spatial_min_cutoff_custom,
    max.cutoff = spatial_max_cutoff_custom,
    pt.size.factor = spatial_pt_size_factor_custom,
    image.alpha = image_alpha_custom,
    image.scale = image_scale_custom
  ) +
    scale_fill_gradient(
      low = spatial_color_low,
      high = spatial_color_high
    ) +
    theme(
      legend.title = element_text(size = legend_title_size_custom),
      legend.text = element_text(size = legend_text_size_custom),
      legend.key.size = unit(legend_key_size_custom, "cm")
    )
  
  ggsave(
    filename = file.path(
      out_dir,
      paste0(name_all_spatial_custom_prefix, "_", gene, ".pdf")
    ),
    plot = current_spatial_custom_plot,
    width = w_spatial_gene_custom,
    height = h_spatial_gene_custom,
    device = "pdf"
  )
}

# 4.6 记录自定义颜色Spatial基因分布图参数

param_step4_text <- paste0(
  "Step 4：Spatial基因分布图（自定义颜色版）\n",
  "运行流程说明：\n",
  "1. 使用Seurat::SpatialFeaturePlot()生成空间表达图。\n",
  "2. 使用scale_fill_gradient()设置表达量渐变颜色。\n",
  "3. 使用theme()调整图例标题、图例文字和图例色块大小。\n",
  "4. 保存当前选中基因图，并批量导出全部基因的自定义颜色Spatial图。\n\n",
  "本次运行参数：\n",
  "- 当前基因：", selected_gene_spatial_custom, "\n",
  "- low：", spatial_color_low, "\n",
  "- high：", spatial_color_high, "\n",
  "- min.cutoff：", spatial_min_cutoff_custom, "\n",
  "- max.cutoff：",
  ifelse(is.na(spatial_max_cutoff_custom), "NA", spatial_max_cutoff_custom), "\n",
  "- pt.size.factor：", spatial_pt_size_factor_custom, "\n",
  "- image.alpha：", image_alpha_custom, "\n",
  "- image.scale：", image_scale_custom, "\n",
  "- legend.title.size：", legend_title_size_custom, "\n",
  "- legend.text.size：", legend_text_size_custom, "\n",
  "- legend.key.size(cm)：", legend_key_size_custom, "\n",
  "- 保存宽高：",
  w_spatial_gene_custom, " × ", h_spatial_gene_custom, "英寸\n",
  "- 当前图文件名前缀：", name_spatial_gene_custom, "\n",
  "- 批量导出文件名前缀：", name_all_spatial_custom_prefix, "\n"
)

# 5. 保存参数记录

# 5.1 参数设置

name_params <- "GenePlot_Spatial_parameters"

# 5.2 合并参数记录

parameter_summary <- paste(
  param_step1_text,
  param_step2_text,
  param_step3_text,
  param_step4_text,
  sep = "\n\n------------------------------\n\n"
)

# 5.3 保存参数记录

writeLines(
  parameter_summary,
  con = file.path(out_dir, paste0(name_params, ".txt"))
)