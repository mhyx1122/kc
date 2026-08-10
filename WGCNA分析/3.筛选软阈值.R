library(WGCNA)

# 1. 设置输入参数
h <- 0.85
pit <- 15
pit_width <- 12
pit_hight <- 6
nThreads <- 1

# 2. 设置WGCNA多线程
if (is.numeric(nThreads) && nThreads >= 2) {
  enableWGCNAThreads(nThreads = nThreads)
} else {
  message("木火医学：v：cgxr410")
}

# 3. 设置软阈值筛选范围
powers <- c(1:pit)

# 4. 计算软阈值
sft <- pickSoftThreshold(
  datExpr0,
  powerVector = powers,
  verbose = 5
)

# 5. 绘制并保存软阈值筛选图
pdf(
  file = file.path(output_dir,"4_independence.pdf"),
  width = pit_width,
  height = pit_hight
)
par(mfrow = c(1,2))
cex1 <- 0.9
plot(
  sft$fitIndices[,1],
  -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
  xlab = "Soft Threshold (power)",
  ylab = "Scale Free Topology Model Fit,signed R^2",
  type = "n",
  main = paste("Scale independence")
)
text(
  sft$fitIndices[,1],
  -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
  labels = powers,
  cex = cex1,
  col = "red"
)
abline(h = h,col = "red")
plot(
  sft$fitIndices[,1],
  sft$fitIndices[,5],
  xlab = "Soft Threshold (power)",
  ylab = "Mean Connectivity",
  type = "n",
  main = paste("Mean connectivity")
)
text(
  sft$fitIndices[,1],
  sft$fitIndices[,5],
  labels = powers,
  cex = cex1,
  col = "red"
)
dev.off()

# 6. 查看软阈值筛选结果
print(sft)

# 7. 提取最佳软阈值
softPower <- sft$powerEstimate
result <- sft$fitIndices[,2] * (-sign(sft$fitIndices[,3]))
softPower_row <- sft$powerEstimate
softPower_value <- result[softPower_row]

# 8. 判断自动筛选结果是否可用
if (is.na(softPower)) {
  print("自动筛选软阈值失败，请增大 pit 参数后再次尝试，或手动筛选 0.8 以上的点")
} else if (softPower_value > 0.8) {
  print(paste("自动筛选软阈值有效，R2 值为正数:",softPower_value,"，softPower =",softPower))
} else if (softPower_value < 0) {
  print(paste("自动筛选软阈值方法失效，R2 为负数:",softPower_value))
  print("自动筛选的结果不可用，请手动寻找 R2 值大于 0.8 以上的软阈值。")
} else {
  print(paste("筛选的 R2 值为:",softPower_value,"小于 0.8，请手动调整，softPower =",softPower))
}

# 9. 输出完成提示
message("数据处理完成。软阈值为: ",softPower)