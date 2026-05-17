library(dplyr)

# 1. 先统计 grch38 注释表中的 RNA 类型数量
anno_grch38_biotype_table <- table(anno_grch38$biotype)

# 2. 查看 grch38 基因组中 RNA 类型统计表
print("grch38基因组中RNA类型统计表：")
print(as.data.frame(anno_grch38_biotype_table))

# 3. 设置输入文件路径
counts_geneexp_file <- "表达矩阵文件.csv"

# 4. 设置要提取的 RNA 类型
bulkqc_gene_type <- "protein_coding"

# 5. 如果结果文件夹不存在，就创建
if (!dir.exists("RNA类型提取结果")) {
  dir.create("RNA类型提取结果")
}

# 6. 读取表达矩阵文件
idcounts <- read.csv(
  file = counts_geneexp_file,
  sep = ",",
  header = TRUE,
  check.names = FALSE
)

# 7. 提取注释表
Map <- anno_grch38
Map <- Map[, c(1, 5, 10)]

# 8. 输出提示
print("如果报错了，请先检查你的数据中到底含不含多种类型的RNA")

# 9. 提取基因名和基因类型两列
B <- data.frame(Map[, c(2, 3)])

# 10. 保留指定 RNA 类型的行
B_protein_coding <- filter(B, B[[2]] == bulkqc_gene_type)

# 11. 删除基因名为空或缺失的行
B_protein_coding_cleaned <- B_protein_coding %>%
  filter(!is.na(B_protein_coding[[1]]) & B_protein_coding[[1]] != "")

# 12. 提取要保留的基因名向量
genes_to_keep <- B_protein_coding_cleaned[[1]]

# 13. 在表达矩阵中保留这些基因
idcounts_filtered <- filter(idcounts, idcounts[[1]] %in% genes_to_keep)

# 14. 生成输出文件名
output_file_name <- paste0("RNA类型提取结果/", bulkqc_gene_type, "_gene_data.csv")

# 15. 保存结果文件
write.csv(idcounts_filtered, output_file_name, row.names = FALSE)
