library(edgeR)

# 1. 设置输入文件路径
file_path <- "表达矩阵文件.csv"

# 2. 设置参数
bili0 <- 0.5
zuidiz <- 0.1
threshold <- 0.1

# 3. 读取表达矩阵
data_DEGA <- read.csv(
  file_path,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

# 4. 如果结果文件夹不存在，就创建
if (!dir.exists("去除低表达基因结果")) {
  dir.create("去除低表达基因结果")
}

# 5. 去除低表达基因（Count）
data10 <- data_DEGA
data10 <- data10[apply(cpm(data10), 1, sum) > zuidiz, ]
data10 <- round(as.matrix(data10))

write.csv(
  data10,
  file = file.path("去除低表达基因结果", "去除低表达基因（Count）.csv"),
  row.names = TRUE
)

message("低表达基因（Count）过滤完成！")

# 6. 去除低表达基因（非Count）
gene_means <- rowMeans(data_DEGA)
filtered_data <- data_DEGA[gene_means > threshold, ]

write.csv(
  filtered_data,
  file = file.path("去除低表达基因结果", "去除低表达基因（非Count）.csv"),
  row.names = TRUE
)

message("低表达基因（非Count）过滤完成！")

# 7. 去除高0值比例的基因
expMatrix <- data_DEGA
zero_counts <- rowSums(expMatrix == 0)
threshold_zero <- ncol(expMatrix) * bili0
expMatrix <- expMatrix[zero_counts <= threshold_zero, ]

write.csv(
  expMatrix,
  file = file.path("去除低表达基因结果", "去除高0值比例.csv"),
  row.names = TRUE
)

message("高0值比例基因过滤完成！")

# 8. 输出完成提示
cat("全部处理完成，结果文件已保存到“去除低表达基因结果”文件夹中。\n")