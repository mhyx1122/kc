library(limma)

# 1. 设置输入文件路径
file_path <- "芯片表达矩阵.csv"

# 2. 读取数据文件
data <- read.csv(
  file_path,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

expMatrix <- data

# 3. 如果结果文件夹不存在，就创建
if (!dir.exists("芯片数据标准化")) {
  dir.create("芯片数据标准化")
}

# 4. 先画 Log2 转换前的箱线图
pdf("芯片数据标准化/1.log2转化前箱线图.pdf", width = 9, height = 6)
par(mar = c(7, 7, 2, 2))
boxplot(expMatrix, col = rainbow(ncol(expMatrix)), notch = TRUE, las = 2)
dev.off()

print("已完成 Log2 转换前的箱线图生成")

# 5. 判断是否需要做 Log2 转换
ex <- expMatrix
qx <- as.numeric(quantile(ex, c(0.00, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))
LogC <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)

# 6. 如果需要，则进行 Log2 转换并生成箱线图
if (LogC) {
  expMatrix <- log2(ex + 1)
  
  pdf("芯片数据标准化/2.log2转化后箱线图.pdf", width = 9, height = 6)
  par(mar = c(7, 7, 2, 2))
  boxplot(expMatrix, notch = TRUE, col = rainbow(ncol(expMatrix)), las = 2)
  dev.off()
  
  write.csv(
    expMatrix,
    file = "芯片数据标准化/microarray_genes_log2.csv",
    row.names = TRUE
  )
  
  print("log2 transform finished")
} else {
  print("log2 transform not needed")
}

# 7. 去除样本间系统差异
expMatrix <- normalizeBetweenArrays(expMatrix)

# 8. 画标准化后的箱线图
pdf("芯片数据标准化/2.log2转化后且normalized的箱线图.pdf", width = 9, height = 6)
par(mar = c(7, 7, 2, 2))
boxplot(expMatrix, notch = TRUE, col = rainbow(ncol(expMatrix)), las = 2)
dev.off()

# 9. 保存标准化后的表达矩阵
write.csv(
  expMatrix,
  file = "芯片数据标准化/microarray_genes_log2_normalized.csv",
  row.names = TRUE
)

# 10. 输出完成提示
print("Normalization finished")
cat("芯片数据标准化处理完成。\n")