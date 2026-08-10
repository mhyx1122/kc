library(dplyr)
library(stringr)
library(ggplot2)
library(DESeq2)
library(edgeR)
library(limma)
library(VennDiagram)
library(tinyarray)
library(patchwork)
library(future.apply)

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

# 5. 读取表达矩阵
data_DEGA <- read.csv(exprflie, header = TRUE, row.names = 1)

# 6. 输出文件维度
cat("数据文件行数:", nrow(data_DEGA), "列数:", ncol(data_DEGA), "\n")

# 7. 生成分组信息
group <- factor(c(rep("T", num_T), rep("N", num_N)))
data_DEGA_with_group <- data.frame(id = colnames(data_DEGA), group = group)
write.csv(data_DEGA_with_group, "data_with_group.csv", row.names = FALSE)

# 8. 创建结果目录
if (!dir.exists("1.差异分析/Count")) {
  dir.create("1.差异分析/Count", recursive = TRUE)
}

# 9. 开始 Count 差异分析
if (P值选择 == "P.Value") {
  adjust <- FALSE
} else if (P值选择 == "adj.P.Val") {
  adjust <- TRUE
}

data <- data_DEGA
data10 <- data
data10 <- data[apply(cpm(data), 1, sum) > 0, ]
data10 <- round(as.matrix(data10))

# 10. DESeq2 差异分析
exprSet <- data10
group_list <- group
(colData <- data.frame(row.names = colnames(exprSet), group_list = group_list))

dds <- DESeqDataSetFromMatrix(
  countData = exprSet,
  colData = colData,
  design = ~ group_list
)
dds <- DESeq(dds)
res <- results(dds, contrast = c("group_list", "T", "N"))
resOrdered <- res[order(res$padj), ]
DESeq2_DEG <- as.data.frame(resOrdered)
DESeq2_DEG <- na.omit(DESeq2_DEG)

logFC_cutoff1 <- with(DESeq2_DEG, mean(abs(log2FoldChange)) + 2 * sd(abs(log2FoldChange)))
if (手动设置阈值) {
  logFC_cutoff1 <- cuf
}

DESeq2_adjust <- data.frame(
  gene_symbol = rownames(DESeq2_DEG),
  P.Value = DESeq2_DEG$pvalue,
  adj.P.Val = DESeq2_DEG$padj
)
row.names(DESeq2_adjust) <- DESeq2_adjust$gene_symbol
DESeq2_adjust$gene_symbol <- NULL

k1 <- (DESeq2_adjust[[P值选择]] < 0.05) & (DESeq2_DEG$log2FoldChange < -logFC_cutoff1)
k2 <- (DESeq2_adjust[[P值选择]] < 0.05) & (DESeq2_DEG$log2FoldChange > logFC_cutoff1)
DESeq2_DEG$change <- ifelse(k1, "DOWN", ifelse(k2, "UP", "NOT"))
print(table(DESeq2_DEG$change))

colnames(DESeq2_DEG)[2] <- "logFC"
colnames(DESeq2_DEG)[5] <- "P.Value"
colnames(DESeq2_DEG)[6] <- "adj.P.Val"

write.csv(DESeq2_DEG, file = "1.差异分析/Count/DESeq2_DEG.csv", row.names = TRUE)

DESeq2_log2FC <- data.frame(
  gene_symbol = rownames(DESeq2_DEG),
  change = DESeq2_DEG$change,
  logFC = DESeq2_DEG$logFC
)
write.csv(DESeq2_log2FC, file = "1.差异分析/Count/DESeq2_log2FC.csv", row.names = FALSE)

DESeq2_log2FC_filtered <- DESeq2_log2FC %>%
  filter(change != "NOT")
write.csv(DESeq2_log2FC_filtered, file = "1.差异分析/Count/DESeq2_上下调基因.csv", row.names = FALSE)

# 11. edgeR 差异分析
dge <- DGEList(counts = data10, group = group)
dge$samples$lib.size <- colSums(dge$counts)
dge <- calcNormFactors(dge)
design <- model.matrix(~0 + group)
rownames(design) <- colnames(dge)
colnames(design) <- levels(group)

dge <- estimateGLMCommonDisp(dge, design)
dge <- estimateGLMTrendedDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)

bcv <- dge$common.dispersion
fit <- glmFit(dge, design)
fit2 <- glmLRT(fit, contrast = c(-1, 1))
edgeR_DEG <- topTags(fit2, n = nrow(data10))
edgeR_DEG <- as.data.frame(edgeR_DEG)

logFC_cutoff2 <- with(edgeR_DEG, mean(abs(logFC)) + 2 * sd(abs(logFC)))
if (手动设置阈值) {
  logFC_cutoff2 <- cuf
}

edgeR_adjust <- data.frame(
  gene_symbol = rownames(edgeR_DEG),
  P.Value = edgeR_DEG$PValue,
  adj.P.Val = edgeR_DEG$FDR
)
row.names(edgeR_adjust) <- edgeR_adjust$gene_symbol
edgeR_adjust$gene_symbol <- NULL

k1 <- (edgeR_adjust[[P值选择]] < 0.05) & (edgeR_DEG$logFC < -logFC_cutoff2)
k2 <- (edgeR_adjust[[P值选择]] < 0.05) & (edgeR_DEG$logFC > logFC_cutoff2)
edgeR_DEG$change <- ifelse(k1, "DOWN", ifelse(k2, "UP", "NOT"))
print(table(edgeR_DEG$change))

colnames(edgeR_DEG)[4] <- "P.Value"
colnames(edgeR_DEG)[5] <- "adj.P.Val"
write.csv(edgeR_DEG, file = "1.差异分析/Count/edgeR_DEG.csv", row.names = TRUE)

edgeR_log2FC <- data.frame(
  gene_symbol = rownames(edgeR_DEG),
  change = edgeR_DEG$change,
  logFC = edgeR_DEG$logFC
)
write.csv(edgeR_log2FC, file = "1.差异分析/Count/edgeR_log2FC.csv", row.names = FALSE)

edgeR_log2FC_filtered <- edgeR_log2FC %>%
  filter(change != "NOT")
write.csv(edgeR_log2FC_filtered, file = "1.差异分析/Count/edgeR_上下调基因.csv", row.names = FALSE)

# 12. limma-voom 差异分析
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
rownames(design) <- colnames(data10)
dge <- calcNormFactors(dge)
v <- voom(dge, design, normalize = "quantile")
constrasts <- paste(rev(levels(group)), collapse = "-")
cont.matrix <- makeContrasts(contrasts = constrasts, levels = design)

fit <- lmFit(v, design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
limma_voom_DEG <- topTable(fit2, coef = constrasts, n = Inf)
limma_voom_DEG <- na.omit(limma_voom_DEG)

logFC_cutoff3 <- with(limma_voom_DEG, mean(abs(logFC)) + 2 * sd(abs(logFC)))
if (手动设置阈值) {
  logFC_cutoff3 <- cuf
}

k1 <- (limma_voom_DEG[[P值选择]] < 0.05) & (limma_voom_DEG$logFC < -logFC_cutoff3)
k2 <- (limma_voom_DEG[[P值选择]] < 0.05) & (limma_voom_DEG$logFC > logFC_cutoff3)
limma_voom_DEG$change <- ifelse(k1, "DOWN", ifelse(k2, "UP", "NOT"))
print(table(limma_voom_DEG$change))
write.csv(limma_voom_DEG, file = "1.差异分析/Count/limma_voom_DEG.csv", row.names = TRUE)

limma_log2FC <- data.frame(
  gene_symbol = rownames(limma_voom_DEG),
  change = limma_voom_DEG$change,
  logFC = limma_voom_DEG$logFC
)
write.csv(limma_log2FC, file = "1.差异分析/Count/limma_log2FC.csv", row.names = FALSE)

limma_log2FC_filtered <- limma_log2FC %>%
  filter(change != "NOT")
write.csv(limma_log2FC_filtered, file = "1.差异分析/Count/limma_voom_上下调基因.csv", row.names = FALSE)

print("动态log2FC阈值为：|logFC| > [mean(|logFC|) + 2sd(|logFC|)]")
message("本次分析采用的阈值为：")
message("DEseq2")
print(logFC_cutoff1)
message("EdgeR")
print(logFC_cutoff2)
message("limma-voom")
print(logFC_cutoff3)

# 13. 设置作图参数
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

# 14. 三大 R 包差异基因对比
A <- data.frame(
  DEseq2 = as.integer(table(DESeq2_DEG$change)),
  edgeR = as.integer(table(edgeR_DEG$change)),
  limma_voom = as.integer(table(limma_voom_DEG$change)),
  row.names = c("DOWN", "NOT", "UP")
)
print(A)

# 15. 上调、下调基因分别画维恩图
UP <- function(df) {
  rownames(df)[df$change == "UP"]
}
DOWN <- function(df) {
  rownames(df)[df$change == "DOWN"]
}

up <- intersect(intersect(UP(DESeq2_DEG), UP(edgeR_DEG)), UP(limma_voom_DEG))
down <- intersect(intersect(DOWN(DESeq2_DEG), DOWN(edgeR_DEG)), DOWN(limma_voom_DEG))

up_genes <- list(
  DEseq2 = UP(DESeq2_DEG),
  edgeR = UP(edgeR_DEG),
  limma = UP(limma_voom_DEG)
)
down_genes <- list(
  DEseq2 = DOWN(DESeq2_DEG),
  edgeR = DOWN(edgeR_DEG),
  limma = DOWN(limma_voom_DEG)
)

up.plot <- draw_venn(up_genes, lwd = 2, "UPgene")
down.plot <- draw_venn(down_genes, lwd = 2, "DOWNgene")

veen_count <- (up.plot + down.plot) + plot_layout(guides = "collect")
intersection_genes <- intersect(up, down)
up_intersection_genes <- intersect(intersect(up_genes$DEseq2, up_genes$edgeR), up_genes$limma)
down_intersection_genes <- intersect(intersect(down_genes$DEseq2, down_genes$edgeR), down_genes$limma)

down_intersection_df <- data.frame(gene_symbol = down_intersection_genes)
write.table(down_intersection_df, file = "1.差异分析/Count/down_intersection_genes.csv", sep = ",", row.names = FALSE)

up_intersection_df <- data.frame(gene_symbol = up_intersection_genes)
write.table(up_intersection_df, file = "1.差异分析/Count/up_intersection_genes.csv", sep = ",", row.names = FALSE)

up_intersection_df$change <- "UP"
down_intersection_df$change <- "DOWN"
combined_intersection_df <- bind_rows(up_intersection_df, down_intersection_df)

# 16. 热图
cg1 <- rownames(DESeq2_DEG)[DESeq2_DEG$change != "NOT"]
cg2 <- rownames(edgeR_DEG)[edgeR_DEG$change != "NOT"]
cg3 <- rownames(limma_voom_DEG)[limma_voom_DEG$change != "NOT"]

p1 <- draw_heatmap(data10[cg1, ], group, n_cutoff = heatmap_cutoff, cluster_rows = gene_tree, cluster_cols = sample_tree, clustering_method = "mcquitty", color = heatmap_genecolor, color_an = heatmap_samplecolor, legend = TRUE, annotation_legend = TRUE)
ggsave("1.差异分析/Count/1.DESeq2单独热图.pdf", plot = p1, width = pdf_width, height = pdf_height, units = "cm")

p2 <- draw_heatmap(data10[cg2, ], group, n_cutoff = heatmap_cutoff, cluster_rows = gene_tree, cluster_cols = sample_tree, clustering_method = "mcquitty", color = heatmap_genecolor, color_an = heatmap_samplecolor, legend = TRUE, annotation_legend = TRUE)
ggsave("1.差异分析/Count/2.edgeR单独热图.pdf", plot = p2, width = pdf_width, height = pdf_height, units = "cm")

p3 <- draw_heatmap(data10[cg3, ], group, n_cutoff = heatmap_cutoff, cluster_rows = gene_tree, cluster_cols = sample_tree, clustering_method = "mcquitty", color = heatmap_genecolor, color_an = heatmap_samplecolor, legend = TRUE, annotation_legend = TRUE)
ggsave("1.差异分析/Count/3.limma_voom单独热图.pdf", plot = p3, width = pdf_width, height = pdf_height, units = "cm")

# 17. 火山图
p4 <- draw_volcano(DESeq2_DEG, pkg = 1, logFC_cutoff = logFC_cutoff1, color = volcano_color, adjust = adjust)
ggsave("1.差异分析/Count/4.DESeq2单独火山图.pdf", plot = p4, width = pdf_width, height = pdf_height, units = "cm")

p5 <- draw_volcano(edgeR_DEG, pkg = 2, logFC_cutoff = logFC_cutoff2, color = volcano_color, adjust = adjust)
ggsave("1.差异分析/Count/5.edgeR单独火山图.pdf", plot = p5, width = pdf_width, height = pdf_height, units = "cm")

p6 <- draw_volcano(limma_voom_DEG, pkg = 3, logFC_cutoff = logFC_cutoff3, color = volcano_color, adjust = adjust)
ggsave("1.差异分析/Count/6.limma_voom单独火山图.pdf", plot = p6, width = pdf_width, height = pdf_height, units = "cm")

h1 <- draw_heatmap(data10[cg1, ], group, n_cutoff = heatmap_cutoff, cluster_rows = gene_tree, cluster_cols = sample_tree, clustering_method = "mcquitty", color = heatmap_genecolor, color_an = heatmap_samplecolor)
h2 <- draw_heatmap(data10[cg2, ], group, n_cutoff = heatmap_cutoff, cluster_rows = gene_tree, cluster_cols = sample_tree, clustering_method = "mcquitty", color = heatmap_genecolor, color_an = heatmap_samplecolor)
h3 <- draw_heatmap(data10[cg3, ], group, n_cutoff = heatmap_cutoff, cluster_rows = gene_tree, cluster_cols = sample_tree, clustering_method = "mcquitty", color = heatmap_genecolor, color_an = heatmap_samplecolor)

v1 <- draw_volcano(DESeq2_DEG, pkg = 1, logFC_cutoff = logFC_cutoff1, color = volcano_color, adjust = adjust)
v2 <- draw_volcano(edgeR_DEG, pkg = 2, logFC_cutoff = logFC_cutoff2, color = volcano_color, adjust = adjust)
v3 <- draw_volcano(limma_voom_DEG, pkg = 3, logFC_cutoff = logFC_cutoff3, color = volcano_color, adjust = adjust)

p <- ((h1 + h2 + h3) / (v1 + v2 + v3)) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")

ggsave(
  filename = "1.差异分析/Count/7.heatmap_volcano.png",
  plot = p,
  width = 8, height = 6, units = "in",
  dpi = 600, device = "png"
)

# 18. 维恩图拼图
pdf("1.差异分析/Count/8.veen.pdf", width = venn_pdf_width, height = venn_pdf_height)
plot((up.plot + down.plot) + plot_layout(guides = "collect"))
dev.off()

message("作图完毕，更多生信课程请关注：b站-生信科学家，V：cgxr410")