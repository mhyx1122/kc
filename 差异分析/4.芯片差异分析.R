library(dplyr)
library(ggplot2)
library(limma)
library(tinyarray)
library(patchwork)
library(future)
library(future.apply)
library(grid)

# 1. 设置输入文件路径
exprflie <- "基因表达数据文件.csv"

# 2. 设置分组样本数量
num_T <- 8
num_N <- 8

# 3. 设置阈值参数
手动设置阈值 <- FALSE
cuf <- 1

# 4. 设置 P 值类型
P值选择 <- "P.Value"

# 5. 设置是否标准化
是否标准化 <- TRUE

# 6. 读取表达矩阵
data_DEGA <- read.csv(exprflie, header = TRUE, row.names = 1)

# 7. 输出文件维度
cat("数据文件行数:", nrow(data_DEGA), "列数:", ncol(data_DEGA), "\n")

# 8. 生成分组信息
group <- factor(c(rep("T", num_T), rep("N", num_N)))
data_DEGA_with_group <- data.frame(id = colnames(data_DEGA), group = group)
write.csv(data_DEGA_with_group, "data_with_group.csv", row.names = FALSE)

# 9. 创建结果目录
if (!dir.exists("1.差异分析/microarray")) {
  dir.create("1.差异分析/microarray", recursive = TRUE)
}

if (P值选择 == "P.Value") {
  adjust <- FALSE
} else if (P值选择 == "adj.P.Val") {
  adjust <- TRUE
}

data <- data_DEGA
expMatrix <- data

# 10. 去除含有大量0值的样本
zero_counts <- rowSums(expMatrix == 0)
threshold <- ncol(expMatrix) * 0.5
expMatrix <- expMatrix[zero_counts <= threshold, ]

# 11. log2转化箱线图
pdf("1.差异分析/microarray/1.log2转化前箱线图.pdf", width = 9, height = 6)
par(mar = c(7, 7, 2, 2))
boxplot(expMatrix, col = group, notch = TRUE, las = 2)
dev.off()

ex <- expMatrix
qx <- as.numeric(quantile(ex, c(0.00, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))
LogC <- (qx[6] > 100) || (qx[6] - qx[1] > 50)

if (LogC) {
  expMatrix <- log2(ex + 1)
  pdf("1.差异分析/microarray/2.log2转化后箱线图.pdf", width = 9, height = 6)
  par(mar = c(7, 7, 2, 2))
  boxplot(expMatrix, notch = TRUE, col = group, las = 2)
  dev.off()
  write.csv(expMatrix, file = "1.差异分析/microarray/microarray_genes_log2.csv", row.names = TRUE)
  print("log2 transform finished")
} else {
  print("log2 transform not needed")
}

print("对于芯片数据的处理，大家通常默认为先进行lo2转换，再进行normalizeBetweenArrays标准化")

if (是否标准化 == TRUE) {
  expMatrix <- normalizeBetweenArrays(expMatrix)
  pdf("1.差异分析/microarray/2.log2转化后且normalized的箱线图.pdf", width = 9, height = 6)
  par(mar = c(7, 7, 2, 2))
  boxplot(expMatrix, notch = TRUE, col = group, las = 2)
  dev.off()
  write.csv(expMatrix, file = "1.差异分析/microarray/microarray_genes_log2_normalized.csv", row.names = TRUE)
}

# 12. limma 差异分析
data10 <- expMatrix
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
rownames(design) <- colnames(data10)

constrasts <- paste(rev(levels(group)), collapse = "-")
cont.matrix <- makeContrasts(contrasts = constrasts, levels = design)

fit <- lmFit(data10, design)
fit1 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit1)
limma_DEG <- topTable(fit2, coef = constrasts, n = Inf, adjust.method = "BH")
limma_DEG <- na.omit(limma_DEG)

logFC_cutoff <- with(limma_DEG, mean(abs(logFC)) + 2 * sd(abs(logFC)))
if (手动设置阈值) {
  logFC_cutoff <- cuf
}

message("limma:动态log2FC阈值为：|logFC| > [mean(|logFC|) + 2sd(|logFC|)]")
message("本次分析采用的阈值为：")
print(logFC_cutoff)

k1 <- (limma_DEG[[P值选择]] < 0.05) & (limma_DEG$logFC < -logFC_cutoff)
k2 <- (limma_DEG[[P值选择]] < 0.05) & (limma_DEG$logFC > logFC_cutoff)
limma_DEG$change <- ifelse(k1, "DOWN", ifelse(k2, "UP", "NOT"))
print(table(limma_DEG$change))
write.csv(limma_DEG, file = "1.差异分析/microarray/limma差异分析结果.csv", row.names = TRUE)

limma_log2FC_microarray <- data.frame(gene_symbol = rownames(limma_DEG), change = limma_DEG$change, logFC = limma_DEG$logFC)
write.csv(limma_log2FC_microarray, file = "1.差异分析/microarray/limma_log2FC.csv", row.names = FALSE)

limma_log2FC_fpkm_filtered <- limma_log2FC_microarray %>%
  filter(change != "NOT")
write.csv(limma_log2FC_fpkm_filtered, file = "1.差异分析/microarray/limma_上下调基因.csv", row.names = FALSE)

# 13. 非参数秩和检验
count_norm <- as.data.frame(expMatrix)

plan(multisession)

pvalues <- future_lapply(1:nrow(count_norm), function(i) {
  data <- cbind.data.frame(gene = as.numeric(t(count_norm[i, ])), group)
  p <- suppressWarnings(wilcox.test(gene ~ group, data)$p.value)
  return(p)
})

plan(sequential)

fdr <- p.adjust(unlist(pvalues), method = "BH")
bonferroni <- p.adjust(unlist(pvalues), method = "bonferroni")
log2FC <- rowMeans(count_norm[, group == "T"]) - rowMeans(count_norm[, group == "N"])

RankTest <- data.frame(logFC = log2FC, P.Value = unlist(pvalues), adj.P.Val = fdr, Bonferroni = bonferroni)
rownames(RankTest) <- rownames(count_norm)
RankTest <- na.omit(RankTest)
write.csv(RankTest, file = "1.差异分析/microarray/非参数秩和检验结果.csv")

logFC_cutoff1 <- with(RankTest, mean(abs(logFC)) + 2 * sd(abs(logFC)))
if (手动设置阈值) {
  logFC_cutoff1 <- cuf
}

message("Wilcox:动态log2FC阈值为：|logFC| > [mean(|logFC|) + 2sd(|logFC|)]")
message("本次分析采用的阈值为：")
print(logFC_cutoff1)

RTk1 <- (RankTest[[P值选择]] < 0.05) & (RankTest$logFC < -logFC_cutoff1)
RTk2 <- (RankTest[[P值选择]] < 0.05) & (RankTest$logFC > logFC_cutoff1)
RankTest$change <- ifelse(RTk1, "DOWN", ifelse(RTk2, "UP", "NOT"))
print(table(RankTest$change))

RankTest_log2FC_microarray <- data.frame(gene_symbol = rownames(RankTest), change = RankTest$change, logFC = RankTest$logFC)
write.csv(RankTest_log2FC_microarray, file = "1.差异分析/microarray/Wilco_Test_log2FC.csv", row.names = FALSE)

limma_log2FC_fpkm_filtered <- RankTest_log2FC_microarray %>%
  filter(change != "NOT")
write.csv(limma_log2FC_fpkm_filtered, file = "1.差异分析/microarray/Wilco_Test_上下调基因.csv", row.names = FALSE)

# 14. 设置作图参数
bottom <- 2
left <- 2
top <- 2
right <- 2

down_color <- "#2874C5"
no_change_color <- "grey"
up_color <- "#f87669"

heatmap_genecolor_low <- "#2fa1dd"
heatmap_genecolor_mid <- "white"
heatmap_genecolor_high <- "#f87669"

heatmap_samplecolor_T <- "#2fa1dd"
heatmap_samplecolor_N <- "#f87669"

gene_tree <- TRUE
sample_tree <- FALSE

pdf_width <- 15
pdf_height <- 15
venn_pdf_width <- 12
venn_pdf_height <- 6
heatmap_cutoff <- 2

volcano_color <- c(down_color, no_change_color, up_color)
heatmap_genecolor <- (grDevices::colorRampPalette(c(heatmap_genecolor_low, heatmap_genecolor_mid, heatmap_genecolor_high)))(100)
heatmap_samplecolor <- c(heatmap_samplecolor_T, heatmap_samplecolor_N)

# 15. limma 热图和火山图
cg1 <- rownames(limma_DEG)[limma_DEG$change != "NOT"]
h1 <- draw_heatmap(data10[cg1, ], group, n_cutoff = heatmap_cutoff, cluster_rows = gene_tree, cluster_cols = sample_tree, clustering_method = "mcquitty", color = heatmap_genecolor, color_an = heatmap_samplecolor)
h1 <- h1 + theme(plot.margin = unit(c(0, 1, 0, 0), "cm"))

p1 <- draw_heatmap(data10[cg1, ], group, n_cutoff = heatmap_cutoff, cluster_rows = gene_tree, cluster_cols = sample_tree, clustering_method = "mcquitty", color = heatmap_genecolor, color_an = heatmap_samplecolor, legend = TRUE, annotation_legend = TRUE)
ggsave("1.差异分析/microarray/3.limma单独热图.pdf", plot = p1, width = pdf_width, height = pdf_height, units = "cm")

v1 <- draw_volcano(limma_DEG, pkg = 4, logFC_cutoff = logFC_cutoff, color = volcano_color, adjust = adjust)
ggsave("1.差异分析/microarray/4.limma单独火山图.pdf", plot = v1, width = pdf_height, height = pdf_height, units = "cm")

pdf("1.差异分析/microarray/5.组合作图.pdf", width = 10, height = 5)
par(mar = c(bottom, left, top, right))
plot(v1 + h1)
dev.off()

# 16. Wilcoxon 热图和火山图
cg1 <- rownames(RankTest)[RankTest$change != "NOT"]
h1 <- draw_heatmap(data10[cg1, ], group, n_cutoff = heatmap_cutoff, cluster_rows = gene_tree, cluster_cols = sample_tree, clustering_method = "mcquitty", color = heatmap_genecolor, color_an = heatmap_samplecolor)
h1 <- h1 + theme(plot.margin = unit(c(0, 1, 0, 0), "cm"))

p6 <- draw_heatmap(data10[cg1, ], group, n_cutoff = heatmap_cutoff, cluster_rows = gene_tree, cluster_cols = sample_tree, clustering_method = "mcquitty", color = heatmap_genecolor, color_an = heatmap_samplecolor, legend = TRUE, annotation_legend = TRUE)
ggsave("1.差异分析/microarray/6.Wilcox单独热图.pdf", plot = p6, width = pdf_height, height = pdf_height, units = "cm")

v1 <- draw_volcano(RankTest, pkg = 4, lab = "RankTest", xlab.package = FALSE, logFC_cutoff = logFC_cutoff1, color = volcano_color, adjust = adjust)
ggsave("1.差异分析/microarray/7.Wilcox单独火山图.pdf", plot = v1, width = pdf_height, height = pdf_height, units = "cm")

pdf("1.差异分析/microarray/8.Wilcox组合作图.pdf", width = 10, height = 5)
par(mar = c(bottom, left, top, right))
plot(v1 + h1)
dev.off()