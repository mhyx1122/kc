# 1. 设置输入文件路径
bulk_geneexp_file <- "表达矩阵文件.csv"
bulk_gene_name <- "基因名称文件.csv"

# 2. 读取表达矩阵文件和基因名称文件
rt1 <- read.csv(bulk_geneexp_file, header=TRUE, row.names=1, check.names=FALSE)
GENEname <- read.csv(bulk_gene_name, header=TRUE)

# 3. 获取要提取的基因名
Genes_name <- GENEname[, 1]
cat("1.接下来提取以下基因的表达量：\n")
print(Genes_name)

# 4. 提取表达数据
exp <- rt1[Genes_name, , drop=FALSE]
exp <- as.data.frame(exp)

# 5. 保存提取后的数据到 CSV 文件
write.csv(exp, file = "基因名提取表达量结果/exp_data.csv", row.names = TRUE)