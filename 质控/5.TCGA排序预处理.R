# 1. 设置输入文件路径
TCGA_group_file <- "选择的TCGA表达矩阵文件.csv"

# 2. 读取TCGA表达矩阵文件
data_group <- read.csv(TCGA_group_file, header=TRUE, row.names=1, check.names=FALSE)

# 3. 获取列名
col_names <- colnames(data_group)

# 4. 定义排序函数
sort_function <- function(col_name) {
  # 提取第14和第15个字符并转换为数字
  num <- as.numeric(substr(col_name, 14, 15))
  return(num < 10)
}

# 5. 按排序函数的结果对列名进行排序
sorted_col_names <- c(col_names[sapply(col_names, sort_function)], col_names[!sapply(col_names, sort_function)])

# 6. 重新排序数据框
data_group1 <- data_group[, sorted_col_names]

# 7. 保存为csv文件
write.csv(data_group1, file = "TCGA分组排序结果/after_group_TCGA.csv", row.names = T)