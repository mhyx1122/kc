suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
})

# 1、读取GWAS本地数据

# GWAS数据文件路径，支持txt或tsv文件
gwas_file <- "GWAS数据.txt"

# 文件分隔符
file_sep <- "\t"

if (!file.exists(gwas_file)) {
  stop("没有找到GWAS数据文件：", gwas_file)
}

gwas_data <- data.table::fread(
  file = gwas_file,
  sep = file_sep,
  data.table = TRUE
)


# 2、根据P-value阈值筛选数据

# P-value列名
pval_col <- "p_value"

# P-value筛选阈值
pvalue_threshold <- 5e-8

if (!pval_col %in% colnames(gwas_data)) {
  stop("GWAS数据中没有找到P-value列：", pval_col)
}

# 去除P-value为NaN的行
df1 <- gwas_data[!is.nan(gwas_data[[pval_col]]), ]

# 筛选P-value小于设定阈值的数据
df1 <- subset(df1, df1[[pval_col]] < pvalue_threshold)


# 3、设置表型名称

# 表型名称
phenotype <- "BRCA"

df1$Phenotype <- phenotype
df1 <- as.data.frame(df1)


# 4、设置GWAS数据对应的列名

# SNP列名
snp_col <- "SNP"

# Beta值列名
beta_col <- "beta"

# 标准误差列名
se_col <- "standard_error"

# 次等位基因频率列名
eaf_col <- "eaf"

# 效应等位基因列名
effect_allele_col <- "effect_allele"

# 其他等位基因列名
other_allele_col <- "other_allele"

# 样本量列名
samplesize_col <- "N"

required_columns <- c(
  snp_col,
  beta_col,
  se_col,
  eaf_col,
  effect_allele_col,
  other_allele_col,
  pval_col,
  samplesize_col
)

missing_columns <- setdiff(required_columns, colnames(df1))

if (length(missing_columns) > 0) {
  stop(
    "GWAS数据中缺少以下列：",
    paste(missing_columns, collapse = "、")
  )
}


# 5、格式化为孟德尔随机化暴露数据

exposure_data <- TwoSampleMR::format_data(
  dat = df1,
  type = "exposure",
  phenotype_col = "Phenotype",
  snp_col = snp_col,
  beta_col = beta_col,
  se_col = se_col,
  eaf_col = eaf_col,
  effect_allele_col = effect_allele_col,
  other_allele_col = other_allele_col,
  pval_col = pval_col,
  samplesize_col = samplesize_col
)

exposure_data$id.exposure <- phenotype


# 6、查看格式化后的数据

print(head(exposure_data))


# 7、保存处理结果

# 结果保存文件夹
folder_path <- "exposure_data"

if (!dir.exists(folder_path)) {
  dir.create(folder_path, recursive = TRUE)
}

output_file <- file.path(
  folder_path,
  paste0("gwas_exposure_data_", phenotype, ".csv")
)

write.csv(
  exposure_data,
  file = output_file,
  row.names = FALSE
)

message("GWAS暴露数据处理完成。")
message("筛选后的SNP数量：", nrow(exposure_data))
message("结果保存位置：", normalizePath(output_file, winslash = "/", mustWork = FALSE))