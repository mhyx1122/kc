library(tibble)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pathview)

# 1. 设置输入参数
input_file <- "差异基因表格.csv"
showPath <- "hsa04936"
limit_gene <- 1
output_dir <- "通路图展示"
GO_min <- org.Hs.eg.db
KEGG_min <- "hsa"

# 2. 读取差异基因表格
diffGene <- read.csv(
  input_file,
  header = TRUE,
  stringsAsFactors = FALSE
)

# 3. 设置颜色区间
user_limit <- list(
  gene = limit_gene,
  cpd = limit_gene
)

# 4. 统一列名
colnames(diffGene)[colnames(diffGene) == "gene_symbol"] <- "gene"
colnames(diffGene)[colnames(diffGene) == "logFC"] <- "logFC"

# 5. 基因Symbol转换为ENTREZID
entrezID <- suppressWarnings(
  suppressMessages(
    bitr(
      diffGene$gene,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = GO_min
    )
  )
)

# 6. 合并差异基因表和ENTREZID
diffGene_merged <- merge(
  diffGene,
  entrezID,
  by.x = "gene",
  by.y = "SYMBOL"
)

# 7. 构建pathview需要的geneFC向量
geneFC <- deframe(
  diffGene_merged[, c("ENTREZID", "logFC")]
)

# 8. 创建输出文件夹
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# 9. 生成KEGG通路图，并直接保存到指定文件夹
pathView <- pathview(
  pathway.id = showPath,
  gene.data = geneFC,
  cpd.data = NULL,
  species = KEGG_min,
  limit = user_limit,
  out.suffix = "pathview",
  same.layer = FALSE,
  keys.align = "y",
  kegg.dir = output_dir
)

# 10. 输出完成提示
message("通路图已生成，保存在文件夹: ", output_dir)