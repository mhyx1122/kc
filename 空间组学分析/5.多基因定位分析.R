suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

# 1. 设置输出文件夹

out_dir <- "基因共定位可视化"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# 2. 读取空间转录组对象并解析基因

# 2.1 参数设置

genes_input <- "BCL2,CD8A"

# 表达量大于该阈值时记为表达
expr_cutoff <- 0

# 2.2 读取Spatial_Data对象

if (!exists("Spatial_Data", envir = .GlobalEnv)) {
  stop("全局环境中不存在 Spatial_Data，请先完成前面的分析步骤。")
}

srt <- get("Spatial_Data", envir = .GlobalEnv)

if (!inherits(srt, "Seurat")) {
  stop("Spatial_Data 不是有效的Seurat对象。")
}

# 2.3 解析输入基因

gene_names <- unlist(strsplit(genes_input, ","))
gene_names <- trimws(gene_names)
gene_names <- gene_names[gene_names != ""]
gene_names <- unique(gene_names)

if (length(gene_names) < 2) {
  stop("至少需要输入2个基因。")
}

# 3. 计算多基因表达组合

# 3.1 提取基因表达数据

expr_df <- FetchData(
  srt,
  vars = gene_names
)

missing_genes <- setdiff(
  gene_names,
  colnames(expr_df)
)

if (length(missing_genes) > 0) {
  stop(
    paste0(
      "以下基因未能成功提取：",
      paste(missing_genes, collapse = ", ")
    )
  )
}

# 3.2 根据表达阈值进行二值化

expr_binary <- as.data.frame(
  expr_df > expr_cutoff
)

expr_binary[] <- lapply(
  expr_binary,
  as.integer
)

# 3.3 将每个基因的二值表达结果写入Seurat对象

for (gene in gene_names) {
  srt[[paste0(gene, "_sf")]] <- expr_binary[[gene]]
}

# 3.4 根据每个spot或cell的表达情况生成组合标签

combo_labels <- apply(
  expr_binary,
  1,
  function(binary_values) {
    positive_genes <- gene_names[which(binary_values == 1)]
    
    if (length(positive_genes) == 0) {
      "None"
    } else if (length(positive_genes) == 1) {
      positive_genes
    } else {
      paste(positive_genes, collapse = "+")
    }
  }
)

# 3.5 设置组合标签顺序

combo_levels <- unique(
  c(
    "None",
    sort(
      unique(
        combo_labels[combo_labels != "None"]
      )
    )
  )
)

srt$coexp_multi <- factor(
  combo_labels,
  levels = combo_levels
)

# 3.6 统计每种表达组合的数量

combo_table <- as.data.frame(
  table(srt$coexp_multi),
  stringsAsFactors = FALSE
)

colnames(combo_table) <- c(
  "Combination",
  "Count"
)

combo_table <- combo_table[
  combo_table$Count > 0,
  ,
  drop = FALSE
]

# 3.7 生成分析摘要

summary_text <- paste0(
  "多基因共定位分析已完成。\n",
  "分析基因：", paste(gene_names, collapse = ", "), "\n",
  "基因数量：", length(gene_names), "\n",
  "表达阈值：", expr_cutoff, "\n",
  "有效组合数量：", nrow(combo_table), "\n",
  "总spot/cell数：", ncol(srt), "\n"
)

# 3.8 更新全局环境中的Spatial_Data

assign(
  "Spatial_Data",
  srt,
  envir = .GlobalEnv
)

# 3.9 记录共定位组合计算参数

param_step1_text <- paste0(
  "Step 1：多基因共定位组合计算\n",
  "运行流程说明：\n",
  "1. 从全局环境读取Spatial_Data对象。\n",
  "2. 使用Seurat::FetchData()提取输入基因的表达数据。\n",
  "3. 按表达阈值将每个基因二值化，表达记为1，不表达记为0。\n",
  "4. 对每个spot/cell汇总多基因表达组合，生成coexp_multi分组。\n",
  "5. 将每个基因的二值表达结果和coexp_multi写入Spatial_Data对象。\n",
  "6. 统计每种共定位组合的数量。\n\n",
  "本次运行参数：\n",
  "- 基因列表：", paste(gene_names, collapse = ", "), "\n",
  "- 基因数量：", length(gene_names), "\n",
  "- 表达阈值：", expr_cutoff, "\n",
  "- 组合数量：", nrow(combo_table), "\n"
)

# 4. 绘制多基因共定位空间图

# 4.1 参数设置

spatial_width <- 8
spatial_height <- 6

spatial_label <- FALSE
spatial_label_size <- 5
spatial_pt_size <- 1

spatial_colors <- c(
  "#BDBDBD", "#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF",
  "#8491B4FF", "#91D1C2FF", "#7E6148FF", "#0073C2FF", "#EFC000FF", "#868686FF",
  "#CD534CFF", "#7AA6DCFF", "#003C67FF", "#8F7700FF", "#3B3B3BFF", "#A73030FF",
  "#374E55FF", "#DF8F44FF", "#00A1D5FF", "#B24745FF", "#79AF97FF", "#6A6599FF",
  "#80796BFF", "#BC3C29FF", "#0072B5FF", "#E18727FF", "#20854EFF", "#7876B1FF",
  "#6F99ADFF", "#FFDC91FF", "#EE4C97FF"
)

name_spatial <- "1.多基因共定位空间图"

# 4.2 检查空间图像信息

if (length(Images(srt)) == 0) {
  stop("Spatial_Data中不存在空间图像信息，无法绘制SpatialPlot。")
}

# 4.3 检查颜色数量

combination_levels <- levels(srt$coexp_multi)

if (length(spatial_colors) < length(combination_levels)) {
  stop(
    paste0(
      "颜色数量不足以覆盖所有共定位组合。",
      "当前需要 ",
      length(combination_levels),
      " 个颜色，但只提供了 ",
      length(spatial_colors),
      " 个颜色。"
    )
  )
}

combination_colors <- setNames(
  spatial_colors[seq_along(combination_levels)],
  combination_levels
)

# 4.4 绘制多基因共定位空间图

spatial_plot <- suppressWarnings(
  SpatialPlot(
    srt,
    group.by = "coexp_multi",
    label = spatial_label,
    label.size = spatial_label_size,
    cols = combination_colors,
    pt.size.factor = spatial_pt_size
  )
)

print(spatial_plot)

# 4.5 保存多基因共定位空间图

ggsave(
  filename = file.path(
    out_dir,
    paste0(name_spatial, ".pdf")
  ),
  plot = spatial_plot,
  width = spatial_width,
  height = spatial_height,
  device = "pdf"
)

# 4.6 记录空间图参数

param_step2_text <- paste0(
  "Step 2：共定位空间定位图\n",
  "运行流程说明：\n",
  "1. 使用Seurat::SpatialPlot(group.by = coexp_multi)绘制多基因共定位空间分布图。\n",
  "2. 根据不同共定位组合进行分类着色。\n",
  "3. 使用label、label.size、pt.size.factor和颜色集合控制图形样式。\n\n",
  "本次运行参数：\n",
  "- label：", spatial_label, "\n",
  "- label.size：", spatial_label_size, "\n",
  "- pt.size.factor：", spatial_pt_size, "\n",
  "- colors：", paste(spatial_colors, collapse = ","), "\n",
  "- 保存宽高：", spatial_width, " × ", spatial_height, "英寸\n",
  "- 文件名：", name_spatial, ".pdf\n"
)

# 5. 保存参数记录

# 5.1 参数设置

name_params <- "Gene_CoLocalization_parameters"

# 5.2 合并参数记录

parameter_summary <- paste(
  param_step1_text,
  param_step2_text,
  sep = "\n\n------------------------------\n\n"
)

# 5.3 保存参数记录

writeLines(
  parameter_summary,
  con = file.path(
    out_dir,
    paste0(name_params, ".txt")
  )
)