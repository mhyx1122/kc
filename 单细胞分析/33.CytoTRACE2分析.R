library(Seurat)
library(ggplot2)
library(CytoTRACE2)

# 1. 检查全局环境中是否存在 seurat 对象
if (!exists("seurat", envir = .GlobalEnv)) {
  stop("全局环境中没有 seurat 对象，请先加载。")
}

srt <- get("seurat", envir = .GlobalEnv)
md <- srt@meta.data


# 2. CytoTRACE2 分析参数模块

# 2.1 CytoTRACE2 参数
species <- "human"
slot_type <- "counts"
annotation_column <- "cellType"

# 2.2 检查 meta.data 中是否存在指定注释列
if (!annotation_column %in% colnames(md)) {
  stop(paste0("meta.data 中不存在列：", annotation_column))
}


# 3. 输出目录模块

# 3.1 输出目录参数
out_dir <- "CytoTRACE2分析"

# 3.2 创建输出目录
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}


# 4. 运行 CytoTRACE2 分析模块

# 4.1 运行 CytoTRACE2
cytotrace2_result <- cytotrace2(
  srt,
  species = species,
  is_seurat = TRUE,
  slot_type = slot_type
)


# 5. 构建细胞注释信息模块

# 5.1 提取细胞注释信息
annotation <- data.frame(
  cellType = md[[annotation_column]],
  row.names = rownames(md)
)


# 6. 导出 CytoTRACE2 结果表格模块

# 6.1 提取 CytoTRACE2 结果中的核心列
cytotrace2_meta_export <- data.frame(
  cell_id = rownames(cytotrace2_result@meta.data),
  cytotrace2_result@meta.data[, c(
    "CytoTRACE2_Score",
    "CytoTRACE2_Potency",
    "CytoTRACE2_Relative",
    "preKNN_CytoTRACE2_Score",
    "preKNN_CytoTRACE2_Potency"
  ), drop = FALSE],
  row.names = NULL,
  check.names = FALSE
)

# 6.2 保存 CytoTRACE2 结果表格
write.csv(
  cytotrace2_meta_export,
  file = file.path(out_dir, "CytoTRACE2_export.csv"),
  row.names = FALSE
)


# 7. 生成 CytoTRACE2 图像模块

# 7.1 生成所有 CytoTRACE2 图像
plots <- plotData(
  cytotrace2_result = cytotrace2_result,
  annotation = annotation,
  expression_data = srt,
  is_seurat = TRUE
)


# 8. 保存 CytoTRACE2 UMAP 图模块

# 8.1 CytoTRACE2 UMAP 图参数
w_cytotrace2_umap <- 8
h_cytotrace2_umap <- 8
name_cytotrace2_umap <- "1.CytoTRACE2_UMAP"

# 8.2 提取 CytoTRACE2 UMAP 图
p_cytotrace2_umap <- plots[["CytoTRACE2_UMAP"]]

if (is.null(p_cytotrace2_umap)) {
  stop("未找到图对象：CytoTRACE2_UMAP")
}

print(p_cytotrace2_umap)

# 8.3 保存 CytoTRACE2 UMAP 图
ggsave(
  filename = file.path(out_dir, paste0(name_cytotrace2_umap, ".pdf")),
  plot = p_cytotrace2_umap,
  width = w_cytotrace2_umap,
  height = h_cytotrace2_umap,
  device = "pdf"
)


# 9. 保存 CytoTRACE2 Potency UMAP 图模块

# 9.1 CytoTRACE2 Potency UMAP 图参数
w_cytotrace2_potency_umap <- 8
h_cytotrace2_potency_umap <- 8
name_cytotrace2_potency_umap <- "2.CytoTRACE2_Potency_UMAP"

# 9.2 提取 CytoTRACE2 Potency UMAP 图
p_cytotrace2_potency_umap <- plots[["CytoTRACE2_Potency_UMAP"]]

if (is.null(p_cytotrace2_potency_umap)) {
  stop("未找到图对象：CytoTRACE2_Potency_UMAP")
}

print(p_cytotrace2_potency_umap)

# 9.3 保存 CytoTRACE2 Potency UMAP 图
ggsave(
  filename = file.path(out_dir, paste0(name_cytotrace2_potency_umap, ".pdf")),
  plot = p_cytotrace2_potency_umap,
  width = w_cytotrace2_potency_umap,
  height = h_cytotrace2_potency_umap,
  device = "pdf"
)


# 10. 保存 CytoTRACE2 Relative UMAP 图模块

# 10.1 CytoTRACE2 Relative UMAP 图参数
w_cytotrace2_relative_umap <- 8
h_cytotrace2_relative_umap <- 8
name_cytotrace2_relative_umap <- "3.CytoTRACE2_Relative_UMAP"

# 10.2 提取 CytoTRACE2 Relative UMAP 图
p_cytotrace2_relative_umap <- plots[["CytoTRACE2_Relative_UMAP"]]

if (is.null(p_cytotrace2_relative_umap)) {
  stop("未找到图对象：CytoTRACE2_Relative_UMAP")
}

print(p_cytotrace2_relative_umap)

# 10.3 保存 CytoTRACE2 Relative UMAP 图
ggsave(
  filename = file.path(out_dir, paste0(name_cytotrace2_relative_umap, ".pdf")),
  plot = p_cytotrace2_relative_umap,
  width = w_cytotrace2_relative_umap,
  height = h_cytotrace2_relative_umap,
  device = "pdf"
)


# 11. 保存 Phenotype UMAP 图模块

# 11.1 Phenotype UMAP 图参数
w_phenotype_umap <- 8
h_phenotype_umap <- 8
name_phenotype_umap <- "4.Phenotype_UMAP"

# 11.2 提取 Phenotype UMAP 图
p_phenotype_umap <- plots[["Phenotype_UMAP"]]

if (is.null(p_phenotype_umap)) {
  stop("未找到图对象：Phenotype_UMAP")
}

print(p_phenotype_umap)

# 11.3 保存 Phenotype UMAP 图
ggsave(
  filename = file.path(out_dir, paste0(name_phenotype_umap, ".pdf")),
  plot = p_phenotype_umap,
  width = w_phenotype_umap,
  height = h_phenotype_umap,
  device = "pdf"
)


# 12. 保存 CytoTRACE2 Boxplot by Pheno 图模块

# 12.1 CytoTRACE2 Boxplot by Pheno 图参数
w_boxplot_bypheno <- 8
h_boxplot_bypheno <- 8
name_boxplot_bypheno <- "5.CytoTRACE2_Boxplot_byPheno"

# 12.2 提取 CytoTRACE2 Boxplot by Pheno 图
p_boxplot_bypheno <- plots[["CytoTRACE2_Boxplot_byPheno"]]

if (is.null(p_boxplot_bypheno)) {
  stop("未找到图对象：CytoTRACE2_Boxplot_byPheno")
}

print(p_boxplot_bypheno)

# 12.3 保存 CytoTRACE2 Boxplot by Pheno 图
ggsave(
  filename = file.path(out_dir, paste0(name_boxplot_bypheno, ".pdf")),
  plot = p_boxplot_bypheno,
  width = w_boxplot_bypheno,
  height = h_boxplot_bypheno,
  device = "pdf"
)


# 13. 保存参数记录模块

# 13.1 参数记录文件参数
name_params <- "CytoTRACE2_parameters"

# 13.2 生成参数记录文本
param_text <- paste0(
  "本次分析参数总结：\n",
  "- 数据来源：全局 seurat 对象\n",
  "- species：", species, "\n",
  "- slot_type：", slot_type, "\n",
  "- annotation 列：", annotation_column, "\n",
  "- 输出目录：", out_dir, "\n",
  "写作提示词（自行组装语言，或借助AI组装）：\n",
  "1.基于 CytoTRACE2 对单细胞数据进行发育潜能推断，评估不同细胞的分化状态。\n",
  "2.结合 CytoTRACE2 分值、Potency 分类和 Relative score 的 UMAP 展示细胞状态在低维空间中的分布特征。\n",
  "3.进一步结合表型注释与箱线图比较不同细胞群体之间的潜能差异。"
)

# 13.3 保存参数记录
writeLines(
  param_text,
  con = file.path(out_dir, paste0(name_params, ".txt"))
)


# 14. 完成提示
cat("\nCytoTRACE2 分析完成。\n")
cat("结果表格已保存：", file.path(out_dir, "CytoTRACE2_export.csv"), "\n")
cat("图片和参数文件已保存到目录：", out_dir, "\n")