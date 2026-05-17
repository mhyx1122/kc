library(dplyr)

# 设置输入文件路径
# 请把这里改成你自己的文件路径
bulk_unique_data_file <- "基因注释后的表格.csv"

# 选择去重方法
# 1 = 相加取平均值
# 2 = 平均值最大的一行
# 3 = 只保留第一个出现的
Remove_duplicate_methods <- 2

cat("数据处理中...\n")
cat("分析进度：1.读取基因注释后的数据...\n")

# 读取基因注释后的表格文件
datqc2 <- read.csv(bulk_unique_data_file, stringsAsFactors = FALSE)

# 如果结果文件夹不存在，就创建
if (!dir.exists("去除重复行结果")) {
  dir.create("去除重复行结果")
}

# 第一列改名为 ID
colnames(datqc2)[1] <- "ID"

cat("分析进度：2.空白行名删除\n")

# 删除第一列为空白的行
data1 <- datqc2[!is.na(datqc2$ID) & trimws(datqc2$ID) != "", ]

# 根据 Remove_duplicate_methods 的值执行相应处理
if (Remove_duplicate_methods == 1) {
  cat("分析进度：3.现在进行的是:方法1（相加取平均值）\n")
  
  # 在原始数据中添加一个索引列
  data1$index <- seq_along(data1$ID)
  
  # 按 ID 分组，对各列求平均值
  average_data <- data1 %>%
    group_by(ID) %>%
    summarise(
      index = min(index),
      across(-index, ~ mean(.x, na.rm = TRUE))
    ) %>%
    arrange(index)
  
  # 删除辅助索引列
  average_data$index <- NULL
  
  # 保存结果
  write.csv(average_data, file = "去除重复行结果/average__unique.csv", row.names = FALSE)
  
  result_data <- average_data
  
} else if (Remove_duplicate_methods == 2) {
  cat("分析进度：3.现在进行的是:方法2（平均值最大的一行）\n")
  
  # 计算每行平均表达值
  data1$MeanExpression <- rowMeans(data1[, -1, drop = FALSE])
  # 添加索引列
  data1$index <- seq_along(data1$ID)
  
  # 每个 ID 保留平均表达值最大的一行
  max_mean_data <- data1 %>%
    group_by(ID) %>%
    slice_max(order_by = MeanExpression, n = 1) %>%
    ungroup() %>%
    arrange(index)
  
  # 删除辅助列
  max_mean_data$index <- NULL
  max_mean_data$MeanExpression <- NULL
  
  # 再去一次重
  max_mean_data <- max_mean_data[!duplicated(max_mean_data[, 1]), ]
  
  # 保存结果
  write.csv(max_mean_data, file = "去除重复行结果/maxmean__unique.csv", row.names = FALSE)
  
  result_data <- max_mean_data
  
} else if (Remove_duplicate_methods == 3) {
  cat("分析进度：3.现在进行的是:方法3(只保留第一个出现的)\n")
  
  # 只保留第一个出现的重复项
  data_unique <- data1[!duplicated(data1[, 1]), ]
  
  # 保存结果
  write.csv(data_unique, file = "去除重复行结果/first_unique.csv", row.names = FALSE)
  
  result_data <- data_unique
  
} else {
  stop("Remove_duplicate_methods 只能设置为 1、2 或 3。")
}

cat("分析完成，结果文件已经保存到“去除重复行结果”文件夹中。\n")