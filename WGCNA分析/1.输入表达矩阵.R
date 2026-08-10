library(WGCNA)

# 1. 设置输入参数
input_file <- "表达矩阵.csv"
is_TCGA <- FALSE
IFlog <- TRUE
mean123 <- 0
MAD <- 3000
output_dir <- "4.WGCNA筛选"
method <- "average"
cutli <- 1
cex <- 1
tree_width <- 10
tree_hight <- 10
BianJu <- c(1,3,2,1)
cex.axis <- 1
cex.main <- 1
cex.lab <- 1

# 2. 读取表达矩阵
rt <- read.csv(
  input_file,
  header = TRUE,
  sep = ",",
  check.names = FALSE,
  row.names = 1
)

# 3. 如果是TCGA数据，则截断样本名
if (is_TCGA) {
  colnames(rt) <- substr(colnames(rt),1,12)
}

# 4. 去除重复样本列
duplicated(colnames(rt))
rt <- rt[,!duplicated(colnames(rt))]

# 5. 转置表达矩阵，使行为样本、列为基因
datExpr0 <- as.data.frame(t(rt))

# 6. 检查样本和基因质量
gsg <- goodSamplesGenes(datExpr0,verbose = 3)
if (!gsg$allOK) {
  if (sum(!gsg$goodGenes) > 0) {
    print(paste("Removing genes:",paste(names(datExpr0)[!gsg$goodGenes],collapse = ", ")))
  }
  if (sum(!gsg$goodSamples) > 0) {
    print(paste("Removing samples:",paste(rownames(datExpr0)[!gsg$goodSamples],collapse = ", ")))
  }
  datExpr0 <- datExpr0[gsg$goodSamples,gsg$goodGenes]
}

# 7. 判断并执行log2转换
if (IFlog) {
  ex <- datExpr0
  qx <- as.numeric(quantile(ex,c(0.00,0.25,0.5,0.75,0.99,1.0),na.rm = TRUE))
  LogC <- (qx[5] > 100) ||
    (qx[6] - qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
  if (LogC) {
    datExpr0 <- log2(ex + 1)
    print("log2 transform finished")
  } else {
    print("log2 transform not needed")
  }
} else {
  print("log2 transform skipped")
}

# 8. 创建输出文件夹
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# 9. 转置数据并计算MAD值
datExpr <- as.data.frame(t(datExpr0))
data_mad <- apply(datExpr,1,mad)
datExpr0 <- datExpr[order(data_mad,decreasing = TRUE)[1:MAD],]
datExpr0 <- as.data.frame(t(datExpr0))

# 10. 计算并过滤表达量低的基因
means <- apply(datExpr0,2,mean)
datExpr0 <- datExpr0[,means > mean123]

# 11. 导出过滤后的表达矩阵
filtered_fpkm <- t(datExpr0)
filtered_fpkm <- data.frame(sample = rownames(filtered_fpkm),filtered_fpkm)
rownames(filtered_fpkm) <- NULL
filtered_fpkm_path <- file.path(output_dir,"rt_filter.xls")
write.table(
  filtered_fpkm,
  file = filtered_fpkm_path,
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE,
  sep = "\t"
)

# 12. 样本聚类
sampleTree <- hclust(dist(datExpr0),method = method)

# 13. 保存原始样本聚类图
pdf(
  file = file.path(output_dir,"1_samplenet.pdf"),
  width = tree_width,
  height = tree_hight
)
par(mar = BianJu)
plot(
  sampleTree,
  main = "Sample clustering to detect outliers",
  sub = "",
  xlab = "",
  ylab = "",
  cex.axis = cex.axis,
  cex.main = cex.main,
  cex = cex
)
dev.off()

# 14. 计算树状图高度和裁切位置
heights <- sampleTree$height
hlevel <- quantile(heights,cutli)
message("总树高为：")
print(heights)
message("从树高多少开始裁切：")
print(hlevel)

# 15. 保存带裁切线的样本聚类图
pdf(
  file = file.path(output_dir,"2_samplenet_cut.pdf"),
  width = tree_width,
  height = tree_hight
)
par(mar = BianJu)
plot(
  sampleTree,
  main = "Sample clustering to detect outliers",
  sub = "",
  xlab = "",
  ylab = "",
  cex.axis = cex.axis,
  cex.main = cex.main,
  cex = cex
)
abline(h = hlevel,col = "red")
dev.off()

# 16. 保留非离群样本
clust <- cutreeStatic(sampleTree,cutHeight = hlevel,minSize = 15)
clust_table <- table(clust)
print(clust_table)
max_clust <- as.numeric(names(which.max(clust_table)))
keepSamples <- clust == max_clust
datExpr0 <- datExpr0[keepSamples,]

# 17. 输出完成提示
message("所有输出文件已保存到：",output_dir)