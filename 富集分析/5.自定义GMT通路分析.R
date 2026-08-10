library(ggplot2)
library(clusterProfiler)
library(enrichplot)

# 1. 设置输入文件和输出文件夹
folder_path <- "自定义的GSEA分析"
limma_file <- "GSEA输入表格.csv"
gmt_file <- "自定义基因集.gmt"
num_zdy <- 15
Pvava <- "pvalue"
geneSetID <- "1,2,3"
color_GSEA <- "#FFCC33,#333399,#FF0033"

# 2. 创建结果保存文件夹
if (!dir.exists(folder_path)) {
  dir.create(folder_path)
}

# 3. 读取输入文件
limma_data <- read.csv(
  limma_file,
  header = TRUE,
  check.names = FALSE,
  row.names = NULL
)
gmt <- read.gmt(gmt_file)

# 4. 构建 GSEA 排序向量
logFC_vector <- setNames(
  limma_data$logFC,
  limma_data$gene_symbol
)
logFC_vector <- sort(logFC_vector, decreasing = TRUE)

# 5. 处理 GSEA 主图展示的 geneSetID 参数
gene_set_id <- unlist(strsplit(geneSetID, ","))
if (all(grepl("^[0-9]+$", gene_set_id))) {
  gene_set_id <- as.numeric(gene_set_id)
} else {
  gene_set_id <- gene_set_id
}

# 6. 处理 GSEA 主图颜色参数
gsea_colors <- unlist(strsplit(color_GSEA, ","))

# 7. 运行 GSEA 分析
ges_result <- suppressWarnings(
  GSEA(
    logFC_vector,
    TERM2GENE = gmt,
    pvalueCutoff = 1,
    eps = 0
  )
)
ges_AB <- as.data.frame(ges_result)

# 8. 绘制并保存 GSEA 山峦图
p <- ridgeplot(
  ges_result,
  showCategory = as.numeric(num_zdy),
  fill = Pvava,
  decreasing = TRUE
) +
  theme(
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8)
  )
ggsave(
  filename = file.path(folder_path, "1.GSEA-MSigDB山峦图.pdf"),
  plot = p,
  width = 15,
  height = 15,
  units = "cm"
)
message("\nGSEA-MSigDB山峦图绘制完成，并保存在文件夹中")

# 9. 绘制并保存 GSEA 气泡图
pp <- dotplot(
  ges_result,
  showCategory = num_zdy,
  font.size = 8,
  color = Pvava
)
ggsave(
  filename = file.path(folder_path, "2.GSEA-MSigDB气泡图.pdf"),
  plot = pp,
  width = 15,
  height = 15,
  units = "cm"
)
message("\nGSEA-MSigDB气泡图绘制完成，并保存在文件夹中")

# 10. 绘制并保存 GSEA 主图
p5 <- gseaplot2(
  ges_result,
  geneSetID = gene_set_id,
  color = gsea_colors,
  pvalue_table = TRUE,
  ES_geom = "line"
)
ggsave(
  filename = file.path(folder_path, "3.GSEA-MSigDB主图.pdf"),
  plot = p5,
  width = 30,
  height = 30,
  units = "cm"
)
message("\nGSEA-MSigDB主图绘制完成，并保存在文件夹中")

# 11. 保存 GSEA 结果表
write.csv(
  ges_AB,
  file = file.path(folder_path, "GSEA.result-自定义.csv"),
  row.names = FALSE
)

# 12. 结束提示
message("\n全部分析完成。")