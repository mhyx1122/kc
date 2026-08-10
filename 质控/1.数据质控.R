library(dplyr)

# 1. 设置输入文件路径
# 请把下面两个路径改成你自己的实际文件路径
gse_data_file <- "表达矩阵文件.csv"
gpl_data_file <- "基因注释文件.csv"

# 2. 读取输入文件
cat("数据处理中...\n")
cat("分析进度：1.读取数据文件...\n")

gse_data <- read.csv(gse_data_file, stringsAsFactors = FALSE)
gpl_data <- read.csv(gpl_data_file, stringsAsFactors = FALSE)

# 3. 创建输出目录
if (!dir.exists("注释基因名结果")) {
  dir.create("注释基因名结果")
}

# 4. 去除前后空字符串
cat("分析进度：2.去除前后空字符串...\n")

gpl_data <- gpl_data %>%
  mutate(across(everything(), ~ trimws(.)))

# 5. 修改列名
colnames(gse_data)[1] <- "ID"
colnames(gpl_data)[1] <- "ProbeID"
colnames(gpl_data)[2] <- "GENE_SYMBOL"

# 6. 创建映射关系
cat("分析进度：3.创建映射关系用以注释...\n")
mapping <- gpl_data[, c("ProbeID", "GENE_SYMBOL")]

# 7. 合并数据
cat("分析进度：4.合并数据...\n")

merged_data <- merge(gse_data, mapping, by.x = "ID", by.y = "ProbeID")

# 8. 用 GENE_SYMBOL 替换 ID 列
final_data <- merged_data
final_data$ID <- final_data$GENE_SYMBOL


# 9. 删除多余列
final_data$GENE_SYMBOL <- NULL


# 10. 保存结果
write.csv(final_data, "注释基因名结果/基因注释后的表格.csv", row.names = FALSE)

cat("分析进度：5.结果文件已经保存到文件夹中了...\n")
cat("此处不需要写作\n")
cat("数据处理完成。\n")
