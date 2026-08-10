# 1. 加载R包

suppressPackageStartupMessages({
  library(ieugwasr)
  library(gwasglue)
  library(coloc)
  library(locuscomparer)
})


# 2. 设置两套GWAS数据和共定位区域参数

# 第一套GWAS数据ID
QTLdata1 <- "ieu-a-300"

# 第二套GWAS数据ID
GWASdata2 <- "ieu-a-7"

# 最显著SNP上下游区域范围，单位为bp
SNP_range <- 90000

if (!is.numeric(SNP_range) || length(SNP_range) != 1 || is.na(SNP_range) || SNP_range <= 0) {
  stop("SNP_range必须是一个大于0的数值。")
}


# 3. 创建结果文件夹

folder_name <- paste0(QTLdata1, "__和__", GWASdata2)
dir.create(folder_name, showWarnings = FALSE, recursive = TRUE)


# 4. 提取第一套GWAS数据中最显著的SNP

top <- suppressWarnings(ieugwasr::tophits(QTLdata1))

if (is.null(top) || nrow(top) == 0) {
  stop("没有从第一套GWAS数据中提取到显著SNP。")
}

required_top_cols <- c("chr", "position", "p", "rsid")
missing_top_cols <- setdiff(required_top_cols, colnames(top))

if (length(missing_top_cols) > 0) {
  stop("tophits结果缺少以下必要列：", paste(missing_top_cols, collapse = "、"))
}

top$p <- suppressWarnings(as.numeric(top$p))
top$position <- suppressWarnings(as.numeric(top$position))
top <- top[!is.na(top$p) & !is.na(top$position) & !is.na(top$chr) & !is.na(top$rsid), , drop = FALSE]

if (nrow(top) == 0) {
  stop("去除缺失值后，没有可用于确定共定位区域的显著SNP。")
}

top <- top[order(top$p), , drop = FALSE]

lead_chr <- as.character(top$chr[1])
lead_position <- as.numeric(top$position[1])
lead_snp <- as.character(top$rsid[1])

region_start <- max(1, lead_position - SNP_range)
region_end <- lead_position + SNP_range
chrpos <- paste0(lead_chr, ":", region_start, "-", region_end)

message("提取到的最显著SNP为：", lead_snp)
message("该SNP位于", lead_chr, "号染色体，基因组位置为：", lead_position)
message("最终选择的共定位区域为：", chrpos)
message("请根据实际研究需要判断该共定位区域是否需要调整。")


# 5. 设置第一套GWAS数据类型和样本量

# 可设置为"quant"或"cc"
GWAS1_size <- "quant"

# 仅当GWAS1_size为"quant"时使用
N_colo1 <- 173082

if (!GWAS1_size %in% c("quant", "cc")) {
  stop("GWAS1_size只能设置为'quant'或'cc'。")
}

if (GWAS1_size == "quant" && (!is.numeric(N_colo1) || length(N_colo1) != 1 || is.na(N_colo1) || N_colo1 <= 0)) {
  stop("第一套GWAS为连续型数据时，必须设置有效的N_colo1。")
}


# 6. 设置第二套GWAS数据类型和样本量

# 可设置为"quant"或"cc"
GWAS2_size <- "cc"

# 仅当GWAS2_size为"quant"时使用；病例对照数据可设为NULL
N_colo2 <- NULL

if (!GWAS2_size %in% c("quant", "cc")) {
  stop("GWAS2_size只能设置为'quant'或'cc'。")
}

if (GWAS2_size == "quant" && (!is.numeric(N_colo2) || length(N_colo2) != 1 || is.na(N_colo2) || N_colo2 <= 0)) {
  stop("第二套GWAS为连续型数据时，必须设置有效的N_colo2。")
}


# 7. 提取并协调两套GWAS区域数据

colodata_GG <- suppressWarnings(
  gwasglue::ieugwasr_to_coloc(
    id1 = QTLdata1,
    id2 = GWASdata2,
    chrompos = chrpos,
    type1 = GWAS1_size,
    type2 = GWAS2_size
  )
)

if (is.null(colodata_GG) || !all(c("dataset1", "dataset2") %in% names(colodata_GG))) {
  stop("没有获得完整的共定位数据，结果中缺少dataset1或dataset2。")
}

dataset1 <- colodata_GG$dataset1
dataset2 <- colodata_GG$dataset2

required_dataset_cols <- c("snp", "pos", "beta", "varbeta")

missing_dataset1_cols <- setdiff(required_dataset_cols, names(dataset1))
missing_dataset2_cols <- setdiff(required_dataset_cols, names(dataset2))

if (length(missing_dataset1_cols) > 0) {
  stop("第一套共定位数据缺少以下必要字段：", paste(missing_dataset1_cols, collapse = "、"))
}

if (length(missing_dataset2_cols) > 0) {
  stop("第二套共定位数据缺少以下必要字段：", paste(missing_dataset2_cols, collapse = "、"))
}

if (length(dataset1$snp) == 0 || length(dataset2$snp) == 0) {
  stop("指定区域内没有可用于共定位分析的共同SNP。")
}

if (length(dataset1$snp) != length(dataset2$snp)) {
  stop("两套共定位数据的SNP数量不一致，无法继续分析。")
}

if (!all(as.character(dataset1$snp) == as.character(dataset2$snp))) {
  stop("两套共定位数据的SNP顺序不一致，无法继续分析。")
}

message("共提取到", length(dataset1$snp), "个可用于共定位分析的共同SNP。")


# 8. 准备第一套GWAS的coloc输入数据

gwas1data <- list(
  snp = as.character(dataset1$snp),
  position = as.numeric(dataset1$pos),
  beta = as.numeric(dataset1$beta),
  varbeta = as.numeric(dataset1$varbeta),
  type = GWAS1_size
)

if (GWAS1_size == "quant") {
  if (is.null(dataset1$MAF) || length(dataset1$MAF) == 0) {
    stop("第一套连续型GWAS数据缺少MAF，无法运行coloc.abf。")
  }
  
  gwas1data$MAF <- as.numeric(dataset1$MAF)
  gwas1data$N <- as.numeric(N_colo1)
}


# 9. 准备第二套GWAS的coloc输入数据

gwas2data <- list(
  snp = as.character(dataset2$snp),
  position = as.numeric(dataset2$pos),
  beta = as.numeric(dataset2$beta),
  varbeta = as.numeric(dataset2$varbeta),
  type = GWAS2_size
)

if (GWAS2_size == "quant") {
  if (is.null(dataset2$MAF) || length(dataset2$MAF) == 0) {
    stop("第二套连续型GWAS数据缺少MAF，无法运行coloc.abf。")
  }
  
  gwas2data$MAF <- as.numeric(dataset2$MAF)
  gwas2data$N <- as.numeric(N_colo2)
}


# 10. 检查共定位分析所需数值

if (anyNA(gwas1data$position) || anyNA(gwas1data$beta) || anyNA(gwas1data$varbeta)) {
  stop("第一套GWAS数据的position、beta或varbeta中存在缺失值。")
}

if (anyNA(gwas2data$position) || anyNA(gwas2data$beta) || anyNA(gwas2data$varbeta)) {
  stop("第二套GWAS数据的position、beta或varbeta中存在缺失值。")
}

if (any(gwas1data$varbeta <= 0) || any(gwas2data$varbeta <= 0)) {
  stop("varbeta中存在小于或等于0的数值，无法运行coloc.abf。")
}

if (GWAS1_size == "quant" && (anyNA(gwas1data$MAF) || any(gwas1data$MAF <= 0 | gwas1data$MAF >= 1))) {
  stop("第一套连续型GWAS数据的MAF存在缺失值或不在0至1之间。")
}

if (GWAS2_size == "quant" && (anyNA(gwas2data$MAF) || any(gwas2data$MAF <= 0 | gwas2data$MAF >= 1))) {
  stop("第二套连续型GWAS数据的MAF存在缺失值或不在0至1之间。")
}


# 11. 运行共定位分析

result <- suppressWarnings(
  coloc::coloc.abf(
    dataset1 = gwas1data,
    dataset2 = gwas2data
  )
)

if (is.null(result$summary) || is.null(result$results)) {
  stop("coloc.abf没有返回完整的共定位结果。")
}

if (!"SNP.PP.H4" %in% colnames(result$results)) {
  stop("共定位结果中没有SNP.PP.H4列。")
}


# 12. 整理共定位结果

# SNP共享因果变异后验概率筛选阈值
SNP_PP_H4_threshold <- 0.75

need_result <- result$results[order(result$results$SNP.PP.H4, decreasing = TRUE, na.last = TRUE), , drop = FALSE]

need_result_sig <- need_result[
  !is.na(need_result$SNP.PP.H4) & need_result$SNP.PP.H4 > SNP_PP_H4_threshold,
  ,
  drop = FALSE
]


# 13. 保存共定位结果

write.csv(result$summary, file.path(folder_name, "coloc_summary_result.csv"))
write.csv(need_result, file.path(folder_name, "coloc_SNP_result.csv"), row.names = FALSE)
write.csv(need_result_sig, file.path(folder_name, "coloc_SNP_result_sig.csv"), row.names = FALSE)

print(result$summary)

if ("PP.H4.abf" %in% names(result$summary)) {
  message("共定位分析的PP.H4为：", signif(result$summary["PP.H4.abf"], 4))
}

message("SNP.PP.H4大于", SNP_PP_H4_threshold, "的SNP数量为：", nrow(need_result_sig))


# 14. 准备LocusCompare绘图数据

if (is.null(dataset1$pvalues) || is.null(dataset2$pvalues)) {
  stop("共定位数据中缺少pvalues，无法绘制LocusCompare图。")
}

gwas_1 <- data.frame(
  rsid = as.character(dataset1$snp),
  pval = suppressWarnings(as.numeric(dataset1$pvalues)),
  stringsAsFactors = FALSE
)

gwas_2 <- data.frame(
  rsid = as.character(dataset2$snp),
  pval = suppressWarnings(as.numeric(dataset2$pvalues)),
  stringsAsFactors = FALSE
)

gwas_1 <- gwas_1[
  !is.na(gwas_1$rsid) & gwas_1$rsid != "" &
    !is.na(gwas_1$pval) & gwas_1$pval > 0 & gwas_1$pval <= 1,
  ,
  drop = FALSE
]

gwas_2 <- gwas_2[
  !is.na(gwas_2$rsid) & gwas_2$rsid != "" &
    !is.na(gwas_2$pval) & gwas_2$pval > 0 & gwas_2$pval <= 1,
  ,
  drop = FALSE
]

if (nrow(gwas_1) == 0 || nrow(gwas_2) == 0) {
  stop("整理后没有可用于LocusCompare绘图的SNP数据。")
}


# 15. 设置LocusCompare图参数

title1 <- "GWAS1"
title2 <- "GWAS2"

# 图片尺寸，单位为英寸
width_coloc <- 9.6
height_coloc <- 6


# 16. 绘制LocusCompare图

locuscompare_plot <- locuscomparer::locuscompare(
  in_fn1 = gwas_1,
  in_fn2 = gwas_2,
  title1 = title1,
  title2 = title2
)


# 17. 保存LocusCompare图

output_plot_file <- file.path(folder_name, "coloc_locuscompare_plot.pdf")

pdf(output_plot_file, width = width_coloc, height = height_coloc)
print(locuscompare_plot)
dev.off()


# 18. 保留主要共定位结果

result_coloc <- need_result

message("共定位分析完成，所有结果已保存至：", folder_name)
