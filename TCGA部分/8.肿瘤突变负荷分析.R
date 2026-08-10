# 1. 加载必要 R 包

library(TCGAbiolinks)
library(maftools)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggExtra)
library(limma)


# 2. 设置 TCGA TMB 下载与计算参数

project_id <- "TCGA-COAD"

capture_size <- 30

folder_name <- "TMB肿瘤突变负荷分析"

plot_width <- 10

plot_height <- 8


# 3. 设置表达量与 TMB 相关性分析参数

exp_file <- "after_group_TCGA.csv"

tmb_file <- file.path(
  folder_name,
  paste0("TMBscore_", project_id, "_处理后结果.csv")
)

TBMgene_name <- "PAF1"

cor_method <- "pearson"
# 可选：
# cor_method <- "spearman"
# cor_method <- "kendall"

folder_TMB2 <- "TMB肿瘤突变负荷分析"

topngene <- 10

cor_plot_width <- 5

cor_plot_height <- 5

x_color <- "orange"

y_color <- "blue"


# 4. 创建结果保存文件夹

if (!dir.exists(folder_name)) {
  dir.create(folder_name, recursive = TRUE)
}

if (!dir.exists(folder_TMB2)) {
  dir.create(folder_TMB2, recursive = TRUE)
}


# 5. 从 TCGA 下载 MAF 突变数据

query <- GDCquery(
  project = project_id,
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  access = "open"
)

GDCdownload(query)

maf_data <- GDCprepare(query)

maf <- read.maf(
  maf = maf_data
)


# 6. 计算 TMB

tmb_result <- tmb(
  maf = maf,
  captureSize = capture_size
)

head(tmb_result)


# 7. 保存 TMB 原始结果

csv_file <- file.path(
  folder_name,
  paste0("TMB_results_", project_id, "_原始数据.csv")
)

write.csv(
  tmb_result,
  file = csv_file,
  row.names = TRUE
)


# 8. 提取 TMBscore 并保存

TMBscore <- tmb_result %>%
  dplyr::select(
    Tumor_Sample_Barcode,
    total_perMB
  )

csv_file_tmbscore <- file.path(
  folder_name,
  paste0("TMBscore_", project_id, "_处理后结果.csv")
)

write.csv(
  TMBscore,
  file = csv_file_tmbscore,
  row.names = FALSE
)


# 9. 保存 MAF Summary 和 TMB 分布图

pdf_file <- file.path(
  folder_name,
  paste0("TMB_Plots_", project_id, ".pdf")
)

pdf(
  pdf_file,
  width = plot_width,
  height = plot_height
)

plotmafSummary(
  maf = maf,
  addStat = "median",
  dashboard = TRUE
)

hist(
  tmb_result$total_perMB,
  main = "TMB Distribution",
  xlab = "TMB (mutations/Mb)",
  breaks = 20
)

dev.off()


# 10. 生成突变基因 summary 数据

mutation_summary <- getGeneSummary(maf)

total_samples <- length(
  maf@clinical.data$Tumor_Sample_Barcode
)

mutation_summary$total_sample <- total_samples

mutation_summary$Ratio <- mutation_summary$MutatedSamples / mutation_summary$total_sample

csv_summary_file <- file.path(
  folder_TMB2,
  paste0("瀑布图配套数据_", project_id, ".csv")
)

write.csv(
  mutation_summary,
  file = csv_summary_file,
  row.names = FALSE
)


# 11. 绘制突变瀑布图

oncoplot_pdf_file <- file.path(
  folder_TMB2,
  paste0("突变瀑布图_", project_id, ".pdf")
)

pdf(
  oncoplot_pdf_file,
  width = cor_plot_width,
  height = cor_plot_height
)

oncoplot(
  maf = maf,
  top = topngene
)

dev.off()


# 12. 读取表达矩阵和 TMB 文件

expTMBData <- read.csv(
  exp_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

tmbTMBData <- read.csv(
  tmb_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

rownames(tmbTMBData) <- gsub(
  "-",
  ".",
  rownames(tmbTMBData)
)


# 13. 提取目标基因表达量

rt <- expTMBData

rt <- log2(rt + 1)

rt <- t(rt)

rt_data <- rt[, TBMgene_name, drop = FALSE]

rt_data <- as.matrix(rt_data)

rownames(rt_data) <- substr(
  rownames(rt_data),
  1,
  12
)

data_clean <- avereps(rt_data)


# 14. 整理 TMB 数据

tmb_data <- tmbTMBData

tmb_data <- t(tmb_data)

tmb_data <- t(tmb_data)

colnames(tmb_data)[1] <- "TMB"

rownames(tmb_data) <- substr(
  rownames(tmb_data),
  1,
  12
)

tmb_clean <- avereps(tmb_data)


# 15. 合并表达量和 TMB 数据

sameSample <- intersect(
  row.names(data_clean),
  row.names(tmb_clean)
)

data_clean <- data_clean[
  sameSample,
  ,
  drop = FALSE
]

tmb_clean <- tmb_clean[
  sameSample,
  ,
  drop = FALSE
]

rt1 <- cbind(
  data_clean,
  tmb_clean
)

output_file <- file.path(
  folder_TMB2,
  "肿瘤突变负荷合并后的数据.csv"
)

write.csv(
  rt1,
  file = output_file,
  row.names = TRUE
)


# 16. 目标基因表达量与 TMB 相关性分析

gene <- colnames(rt1)[1]

x <- as.numeric(
  rt1[, gene]
)

y <- log2(
  as.numeric(rt1[, "TMB"]) + 1
)

df1 <- as.data.frame(
  cbind(
    x,
    y
  )
)

corT <- cor.test(
  x,
  y,
  method = cor_method
)

print(corT)


# 17. 绘制相关性散点图和边际密度图

p1 <- ggplot(
  df1,
  aes(
    x,
    y
  )
) +
  xlab(
    paste0(
      gene,
      " expression"
    )
  ) +
  ylab(
    "Tumor mutation burden"
  ) +
  geom_point() +
  geom_smooth(
    method = "lm",
    formula = y ~ x
  ) +
  theme_bw() +
  stat_cor(
    method = cor_method,
    aes(
      x = x,
      y = y
    )
  )

p2 <- ggMarginal(
  p1,
  type = "density",
  xparams = list(
    fill = x_color
  ),
  yparams = list(
    fill = y_color
  )
)


# 18. 保存相关性图

output_pdf <- file.path(
  folder_TMB2,
  "cor.pdf"
)

pdf(
  file = output_pdf,
  width = cor_plot_width,
  height = cor_plot_height
)

print(p2)

dev.off()


# 19. 输出结果路径

cat("TMB 原始结果已保存到：", csv_file, "\n")

cat("TMBscore 结果已保存到：", csv_file_tmbscore, "\n")

cat("TMB 汇总图已保存到：", pdf_file, "\n")

cat("突变 summary 数据已保存到：", csv_summary_file, "\n")

cat("突变瀑布图已保存到：", oncoplot_pdf_file, "\n")

cat("表达量与 TMB 合并数据已保存到：", output_file, "\n")

cat("相关性图已保存到：", output_pdf, "\n")