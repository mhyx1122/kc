library(ggplot2)
library(clusterProfiler)
library(enrichplot)

# 1. 设置输入文件路径
input_file <- "差异分析结果.csv"
gmt_file <- "自定义基因集示例格式.csv"

# 2. 设置输出文件夹
output_folder <- "自定义的GSEA分析"

# 3. 设置展示通路个数
show_category_num <- 15

# 4. 设置 P 值展示方式
# 可选："p.adjust" 或 "pvalue"
pvalue_display <- "pvalue"

# 5. 设置 GSEA 主图展示的基因集
# 可以写数字编号，例如 "1,2,3"
# 也可以写具体基因集 ID，例如 "HALLMARK_TNFA_SIGNALING_VIA_NFKB,HALLMARK_IL6_JAK_STAT3_SIGNALING"
gene_set_id_input <- "1,2,3"

# 6. 设置 GSEA 主图颜色
gsea_colors_input <- "#FFCC33,#333399,#FF0033"

# 7. 创建输出文件夹
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

# 8. 解析 GSEA 主图展示的基因集参数
gene_set_id <- unlist(strsplit(gene_set_id_input, ","))

if (all(grepl("^[0-9]+$", gene_set_id))) {
  gene_set_id <- as.numeric(gene_set_id)
} else {
  gene_set_id <- gene_set_id
}

# 9. 解析 GSEA 主图颜色
gsea_colors <- unlist(strsplit(gsea_colors_input, ","))

# 10. 读取 GSEA 输入表格
limma_data <- read.csv(
  input_file,
  header = TRUE,
  check.names = FALSE,
  row.names = NULL
)

# 11. 读取自定义基因集文件
gmt <- read.csv(
  gmt_file,
  header = TRUE,
  check.names = FALSE,
  row.names = NULL
)

# 12. 构建 GSEA 输入向量
logFC_vector <- setNames(limma_data$logFC, limma_data$gene_symbol)
logFC_vector <- sort(logFC_vector, decreasing = TRUE)

# 13. 运行自定义 GSEA 分析
ges_result <- suppressWarnings(
  GSEA(
    logFC_vector,
    TERM2GENE = gmt,
    pvalueCutoff = 1,
    eps = 0
  )
)

# 14. 提取 GSEA 结果表
ges_AB <- as.data.frame(ges_result)

# 15. 绘制并保存 GSEA 山峦图
p <- ridgeplot(
  ges_result,
  showCategory = as.numeric(show_category_num),
  fill = pvalue_display,
  decreasing = TRUE
) +
  theme(
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8)
  )

ggsave(
  filename = file.path(output_folder, "1.GSEA-MSigDB山峦图.pdf"),
  plot = p,
  width = 15,
  height = 15,
  units = "cm"
)

message("\nGSEA-MSigDB山峦图绘制完成，并保存在文件夹中")

# 16. 绘制并保存 GSEA 气泡图
pp <- dotplot(
  ges_result,
  showCategory = show_category_num,
  font.size = 8,
  color = pvalue_display
)

ggsave(
  filename = file.path(output_folder, "2.GSEA-MSigDB气泡图.pdf"),
  plot = pp,
  width = 15,
  height = 15,
  units = "cm"
)

message("\nGSEA-MSigDB气泡图绘制完成，并保存在文件夹中")

# 17. 绘制并保存 GSEA 主图
p5 <- gseaplot2(
  ges_result,
  geneSetID = gene_set_id,
  color = gsea_colors,
  pvalue_table = TRUE,
  ES_geom = "line"
)

ggsave(
  filename = file.path(output_folder, "3.GSEA-MSigDB主图.pdf"),
  plot = p5,
  width = 30,
  height = 30,
  units = "cm"
)

message("\nGSEA-MSigDB主图绘制完成，并保存在文件夹中")

# 18. 保存 GSEA 结果表

write.csv(
  ges_AB,
  file = file.path(output_folder, "GSEA.result-自定义.csv"),
  row.names = FALSE
)

# 19. 输出完成提示
cat("自定义 GSEA 分析完成，结果已保存到文件夹：", output_folder, "\n")
