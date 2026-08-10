# 1. 设置输入文件路径
geneGroup_data_file <- "表达矩阵文件.csv"

# 2. 设置按哪个基因进行分组
group_by112 <- "AADAT"

# 3. 如果结果文件夹不存在，就创建
if (!dir.exists("中位数分组结果")) {
  dir.create("中位数分组结果")
}

# 4. 读取表达矩阵文件
data_Select_group <- read.csv(
  geneGroup_data_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

data <- data_Select_group

# 5. 开始按中位数分组
prompt <- "按中位数进行分组，开始"
print(prompt)

# 6. 计算目标基因在所有样本中的中位数
median_value <- median(as.numeric(data[group_by112, ]), na.rm = TRUE)

cat("该基因中位数为：")
print(median_value)

# 7. 根据中位数给样本列名加后缀
new_colnames <- ifelse(
  data[group_by112, ] > median_value,
  paste(colnames(data), "_high", sep = ""),
  paste(colnames(data), "_low", sep = "")
)

# 8. 更新列名
colnames(data) <- new_colnames

print("分组后的列名")
print(colnames(data))

# 9. 将 high 组排前面，low 组排后面
data_gp_byMedian <- data[, order(-grepl("_high$", names(data)), grepl("_low$", names(data)))]

print("对列名进行排序")
print(colnames(data_gp_byMedian))

# 10. 保存结果文件
write.csv(
  data_gp_byMedian,
  file = "中位数分组结果/new_group_by_median.csv",
  row.names = TRUE
)

# 11. 输出完成提示
prompt <- "结束，分组后的文件已下载，文件名为new_group_by_median"
prompt1 <- "如需其它分组方法，请加作者微信cgxr410提供思路"
print(prompt)
print(prompt1)