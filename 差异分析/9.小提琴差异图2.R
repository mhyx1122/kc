library(ggplot2)
library(reshape2)

# 1. 设置输入文件路径
file_path <- "基因表达CSV文件.csv"
gene_file_path <- "感兴趣基因CSV文件.csv"

# 2. 设置分组参数
Group1 <- 8
Group2 <- 8
GroupNAME1 <- "T"
GroupNAME2 <- "N"

# 3. 设置分组颜色
color_Group1 <- "#1f78b4"
color_Group2 <- "#e31a1c"

# 4. 设置图像大小
Figwidth <- 9
Figheight <- 5

# 5. 读取数据
rt <- read.csv(file_path, header = TRUE, sep = ",", check.names = FALSE, row.names = 1)

# 6. 判断是否需要 log2 转换
expMatrix <- rt
ex <- expMatrix
qx <- as.numeric(quantile(ex, c(0.00, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))
LogC <- (qx[6] > 100) || (qx[6] - qx[1] > 50)

if (LogC) {
  expMatrix <- log2(ex + 1)
  print("log2 transform finished")
} else {
  print("log2 transform not needed")
}

# 7. 转置表达矩阵并添加分组信息
rt <- t(expMatrix)
group <- factor(c(rep(GroupNAME1, Group1), rep(GroupNAME2, Group2)))
rt <- data.frame(Type = group, rt, check.names = FALSE)

# 8. 读取基因名称文件
geneData <- read.csv(gene_file_path, header = FALSE, stringsAsFactors = FALSE)
geneNAME <- as.character(geneData[, 1])

# 9. 提取感兴趣基因的数据
rt_subset <- rt[, c("Type", geneNAME)]

# 10. 获取 x 轴名称
x <- colnames(rt_subset)[1]

# 11. 差异分析并添加显著性标记
geneSig <- c("")

for (gene in colnames(rt_subset)[2:ncol(rt_subset)]) {
  rt_subset1 <- rt_subset[, c(gene, "Type")]
  colnames(rt_subset1) <- c("expression", "Type")
  
  p <- if (length(levels(factor(rt_subset1$Type))) > 2) {
    kruskal.test(expression ~ Type, data = rt_subset1)$p.value
  } else {
    wilcox.test(expression ~ Type, data = rt_subset1)$p.value
  }
  
  Sig <- ifelse(p < 0.001, "***",
                ifelse(p < 0.01, "**",
                       ifelse(p < 0.05, "*", "")))
  geneSig <- c(geneSig, Sig)
}

colnames(rt_subset) <- paste0(colnames(rt_subset), geneSig)

# 12. 数据转换成长格式
data <- melt(rt_subset, id.vars = c("Type"))
colnames(data) <- c("Type", "Gene", "Expression")

# 13. 设置自定义颜色
custom_colors <- setNames(c(color_Group1, color_Group2), c(GroupNAME1, GroupNAME2))

# 14. 创建输出目录
if (!dir.exists("1.差异分析/指定基因的可视化")) {
  dir.create("1.差异分析/指定基因的可视化", recursive = TRUE)
}

# 15. 生成图形对象
p1 <- ggplot(data, aes(x = Type, y = Expression, fill = Type)) +
  geom_violin(alpha = 0.6, color = "black", linewidth = 0.5) +
  geom_boxplot(width = 0.2, position = position_dodge(0.9), outlier.shape = NA, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = custom_colors) +
  facet_wrap(~Gene, nrow = 1) +
  labs(x = x, y = "Gene Expression-log2(x+1)", fill = "Group") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 16. 输出 PDF 文件
pdf(file = "1.差异分析/指定基因的可视化/vioplot.pdf", width = Figwidth, height = Figheight)
print(p1)
dev.off()

# 17. 输出完成提示
cat("数据处理完成。图像已保存为 PDF 文件：1.差异分析/指定基因的可视化/vioplot.pdf\n")