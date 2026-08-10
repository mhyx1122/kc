suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
})

# 1、读取GWAS结局数据

# GWAS数据文件路径
gwas_file <- "GWAS数据.txt"

# 文件分隔符
file_sep <- "\t"

if (!file.exists(gwas_file)) {
  stop("没有找到GWAS数据文件：", gwas_file)
}

gwas_data_out <- data.table::fread(
  file = gwas_file,
  sep = file_sep,
  data.table = TRUE
)


# 2、设置结局表型名称

# 表型名称
phenotype <- "BRCA"

gwas_data_out$Phenotype <- phenotype
gwas_data_out1 <- as.data.frame(gwas_data_out)


# 3、设置GWAS数据列名

# SNP列名
snp_col <- "SNP"

# Beta值列名
beta_col <- "beta"

# 标准误差列名
se_col <- "standard_error"

# 效应等位基因频率列名
eaf_col <- "eaf"

# 效应等位基因列名
effect_allele_col <- "effect_allele"

# 其他等位基因列名
other_allele_col <- "other_allele"

# P-value列名
pval_col <- "p_value"

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

missing_columns <- setdiff(required_columns, colnames(gwas_data_out1))

if (length(missing_columns) > 0) {
  stop("GWAS数据中缺少以下列：", paste(missing_columns, collapse = "、"))
}


# 4、格式化为孟德尔随机化结局数据

outcome_data <- TwoSampleMR::format_data(
  dat = gwas_data_out1,
  type = "outcome",
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

outcome_data$outcome <- phenotype
outcome_data$id.outcome <- phenotype


# 5、查看格式化后的数据

print(head(outcome_data))


# 6、保存处理结果

# 结果保存文件夹
folder_path <- "outcome_data"

if (!dir.exists(folder_path)) {
  dir.create(folder_path, recursive = TRUE)
}

output_file <- file.path(
  folder_path,
  paste0("outcome_data_", phenotype, ".csv")
)

write.csv(
  outcome_data,
  file = output_file,
  row.names = FALSE
)

message("GWAS结局数据处理完成。")
message("结局数据SNP数量：", nrow(outcome_data))
message("结果保存位置：", normalizePath(output_file, winslash = "/", mustWork = FALSE))