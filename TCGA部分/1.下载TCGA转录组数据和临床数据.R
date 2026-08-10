# 1. 加载必要 R 包
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(stringr)

# 2. 设置 TCGA 项目 ID
TumorName <- "TCGA-LIHC"

# 3. 设置表达数据输出文件夹
folder_name <- "exp"

if (!dir.exists(folder_name)) {
  dir.create(folder_name)
}

# 4. 下载 TCGA 表达数据
expquery <- GDCquery(
  project = TumorName,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(
  expquery,
  files.per.chunk = 40
)

saveRDS(
  expquery,
  paste0(TumorName, ".rds")
)

expr123 <- GDCprepare(expquery)

print("下一步提取：counts/TPM/FPKM数据")

# 5. 提取 counts 数据并保存
counts <- as.data.frame(assay(expr123))

data <- as.data.frame(rowRanges(expr123))

expr123_count <- cbind(
  gene_type = data$gene_type,
  gene_name = data$gene_name,
  counts
)

write.csv(
  expr123_count,
  file = paste0(folder_name, "/all_gene_counts.csv"),
  row.names = TRUE
)

# 6. 提取 TPM 数据并保存
TPM <- as.data.frame(
  assay(expr123, i = "tpm_unstrand")
)

expr123_TPM <- cbind(
  gene_type = data$gene_type,
  gene_name = data$gene_name,
  TPM
)

write.csv(
  expr123_TPM,
  file = paste0(folder_name, "/all_gene_TPM.csv"),
  row.names = TRUE
)

# 7. 提取 FPKM 数据并保存
FPKM <- as.data.frame(
  assay(expr123, i = "fpkm_unstrand")
)

expr123_FPKM <- cbind(
  gene_type = data$gene_type,
  gene_name = data$gene_name,
  FPKM
)

write.csv(
  expr123_FPKM,
  file = paste0(folder_name, "/all_gene_FPKM.csv"),
  row.names = TRUE
)

# 8. 计算每种 gene_type 的出现次数并保存
gene_type_counts <- data %>%
  group_by(gene_type) %>%
  summarise(count = n()) %>%
  arrange(desc(count))

write.csv(
  gene_type_counts,
  file = paste0(folder_name, "/gene_type_counts.csv"),
  row.names = FALSE
)

# 9. 读取 gene_type 统计文件
gene_type_counts <- read.csv(
  "exp/gene_type_counts.csv",
  header = TRUE,
  sep = ",",
  check.names = FALSE
)

# 10. 设置需要筛选的 GeneType
GeneType <- c("protein_coding")

# 11. 读取 counts / TPM / FPKM 文件
expr_count <- read.csv(
  "exp/all_gene_counts.csv",
  header = TRUE,
  row.names = 1,
  stringsAsFactors = FALSE
)

expr_TPM <- read.csv(
  "exp/all_gene_TPM.csv",
  header = TRUE,
  row.names = 1,
  stringsAsFactors = FALSE
)

expr_FPKM <- read.csv(
  "exp/all_gene_FPKM.csv",
  header = TRUE,
  row.names = 1,
  stringsAsFactors = FALSE
)

# 12. 根据 GeneType 筛选数据
mRNA_count <- filter(
  expr_count,
  gene_type %in% GeneType
)

mRNA_tpm <- filter(
  expr_TPM,
  gene_type %in% GeneType
)

mRNA_fpkm <- filter(
  expr_FPKM,
  gene_type %in% GeneType
)

# 13. 保存筛选后的数据
output_count_file <- paste0(
  "exp/mRNA_count_",
  paste(GeneType, collapse = "_"),
  ".csv"
)

output_tpm_file <- paste0(
  "exp/mRNA_tpm_",
  paste(GeneType, collapse = "_"),
  ".csv"
)

output_fpkm_file <- paste0(
  "exp/mRNA_fpkm_",
  paste(GeneType, collapse = "_"),
  ".csv"
)

write.csv(
  mRNA_count,
  file = output_count_file,
  row.names = TRUE
)

write.csv(
  mRNA_tpm,
  file = output_tpm_file,
  row.names = TRUE
)

write.csv(
  mRNA_fpkm,
  file = output_fpkm_file,
  row.names = TRUE
)

# 14. 设置表达矩阵去重输入文件
bulk_unique_data_file <- "exp/mRNA_tpm_protein_coding.csv"

# 15. 设置去重方法
# 1 = 相加取平均值
# 2 = 平均值最大的一行
# 3 = 只保留第一个出现的
Remove_duplicate_methods <- 1

# 16. 创建表达矩阵去重输出目录
file_name <- tools::file_path_sans_ext(
  basename(bulk_unique_data_file)
)

output_dir <- file.path(
  "1.数据准备",
  file_name
)

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 17. 定义列排序函数
sort_columns <- function(data_group) {
  col_names <- colnames(data_group)
  
  sort_function <- function(col_name) {
    num <- as.numeric(substr(col_name, 14, 15))
    return(num < 10)
  }
  
  sorted_col_names <- c(
    col_names[sapply(col_names, sort_function)],
    col_names[!sapply(col_names, sort_function)]
  )
  
  return(sorted_col_names)
}

# 18. 读取需要去重的表达矩阵
data <- read.csv(bulk_unique_data_file)

data$X <- NULL
data$gene_type <- NULL

colnames(data)[1] <- "ID"

# 19. 删除第一列为空白的行
data1 <- data[data$ID != "", ]

# 20. 去重方法 1：相加取平均值
if (Remove_duplicate_methods == 1) {
  data1$index <- seq_along(data1$ID)
  
  average_data <- data1 %>%
    group_by(ID) %>%
    summarise(across(everything(), mean)) %>%
    arrange(index)
  
  average_data$index <- NULL
  
  write.csv(
    average_data,
    file = file.path(output_dir, "average__unique.csv"),
    row.names = FALSE
  )
  
  data_group <- read.csv(
    file.path(output_dir, "average__unique.csv"),
    header = TRUE,
    row.names = 1,
    check.names = FALSE
  )
  
  sorted_col_names <- sort_columns(data_group)
  data_group1 <- data_group[, sorted_col_names]
  
  write.csv(
    data_group1,
    file = file.path(output_dir, "after_group_TCGA.csv"),
    row.names = TRUE
  )
  
  message("数据预处理完成：1.average__unique.csv（去除重复基因名后的数据）")
  message("数据预处理完成：2.after_group_TCGA.csv（分组排序后的数据）")
}

# 21. 去重方法 2：平均值最大的一行
if (Remove_duplicate_methods == 2) {
  data1$MeanExpression <- rowMeans(data1[, -1])
  data1$index <- seq_along(data1$ID)
  
  max_mean_data <- data1 %>%
    group_by(ID) %>%
    slice_max(order_by = MeanExpression, n = 1) %>%
    ungroup() %>%
    arrange(index)
  
  max_mean_data$index <- NULL
  max_mean_data$MeanExpression <- NULL
  
  write.csv(
    max_mean_data,
    file = file.path(output_dir, "maxmean__unique.csv"),
    row.names = FALSE
  )
  
  data_group <- read.csv(
    file.path(output_dir, "maxmean__unique.csv"),
    header = TRUE,
    row.names = 1,
    check.names = FALSE
  )
  
  sorted_col_names <- sort_columns(data_group)
  data_group1 <- data_group[, sorted_col_names]
  
  write.csv(
    data_group1,
    file = file.path(output_dir, "after_group_TCGA.csv"),
    row.names = TRUE
  )
  
  message("数据预处理完成：1.maxmean__unique.csv（去除重复基因名后的数据）")
  message("数据预处理完成：2.after_group_TCGA.csv（分组排序后的数据）")
}

# 22. 去重方法 3：只保留第一个出现的
if (Remove_duplicate_methods == 3) {
  data_unique <- data1[!duplicated(data1[, 1]), ]
  
  write.csv(
    data_unique,
    file = file.path(output_dir, "first_unique.csv"),
    row.names = FALSE
  )
  
  data_group <- read.csv(
    file.path(output_dir, "first_unique.csv"),
    header = TRUE,
    row.names = 1,
    check.names = FALSE
  )
  
  sorted_col_names <- sort_columns(data_group)
  data_group1 <- data_group[, sorted_col_names]
  
  write.csv(
    data_group1,
    file = file.path(output_dir, "after_group_TCGA.csv"),
    row.names = TRUE
  )
  
  message("数据预处理完成：1.first_unique.csv（去除重复基因名后的数据）")
  message("数据预处理完成：2.after_group_TCGA.csv（分组排序后的数据）")
}

# 23. 设置临床数据文件
clinical_data_file <- "clinical.tsv"

# 24. 创建临床数据输出文件夹
if (!dir.exists("1.数据准备")) {
  dir.create("1.数据准备")
}

# 25. 读取临床数据
clinical_data <- read.csv(
  clinical_data_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# 26. 提取关键临床字段
clin_time <- clinical_data %>%
  dplyr::select(
    cases.submitter_id,
    demographic.vital_status,
    demographic.days_to_death,
    diagnoses.days_to_last_follow_up,
    demographic.age_at_index,
    demographic.gender,
    diagnoses.ajcc_pathologic_t,
    diagnoses.ajcc_pathologic_m,
    diagnoses.ajcc_pathologic_n,
    diagnoses.ajcc_pathologic_stage
  ) %>%
  dplyr::filter(!duplicated(cases.submitter_id))

# 27. 构建原始生存数据
clin_merge <- clin_time %>%
  dplyr::mutate(
    futime = case_when(
      demographic.vital_status == "Alive" ~ diagnoses.days_to_last_follow_up,
      demographic.vital_status == "Dead" ~ demographic.days_to_death
    )
  ) %>%
  dplyr::mutate(
    fustat = case_when(
      demographic.vital_status == "Alive" ~ 0,
      demographic.vital_status == "Dead" ~ 1
    )
  ) %>%
  dplyr::mutate(
    demographic.gender = case_when(
      demographic.gender == "female" ~ 0,
      demographic.gender == "male" ~ 1
    )
  )

write.csv(
  clin_merge,
  file = "1.数据准备/1.原始生存数据.csv",
  row.names = FALSE
)

# 28. 整理有序变量生存数据
survival_data <- data.frame(
  id = clin_merge$cases.submitter_id,
  fustat = clin_merge$fustat,
  futime = clin_merge$futime,
  age = clin_merge$demographic.age_at_index,
  gender = clin_merge$demographic.gender,
  stage = clin_merge$diagnoses.ajcc_pathologic_stage,
  T = clin_merge$diagnoses.ajcc_pathologic_t,
  N = clin_merge$diagnoses.ajcc_pathologic_n,
  M = clin_merge$diagnoses.ajcc_pathologic_m
)

merged_data <- survival_data %>%
  mutate(
    stage = case_when(
      stage %in% c("Stage I", "Stage IA", "Stage IB", "Stage IC") ~ "1",
      stage %in% c("Stage II", "Stage IIA", "Stage IIB", "Stage IIC") ~ "2",
      stage %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC") ~ "3",
      stage %in% c("Stage IV", "Stage IVA", "Stage IVB", "Stage IVC") ~ "4",
      TRUE ~ stage
    )
  ) %>%
  mutate(
    T = str_extract(T, "[0-9]+"),
    N = str_extract(N, "[0-9]+"),
    M = str_extract(M, "[0-9]+")
  )

merged_data$id <- gsub("-", ".", merged_data$id)
merged_data[merged_data == "'--"] <- NA

write.csv(
  merged_data,
  file = "1.数据准备/2.处理好的生存数据(有序变量).csv",
  row.names = FALSE
)

# 29. 整理因子型生存数据
data_clean <- merged_data

data_clean <- data_clean %>%
  mutate(
    stage = as.character(stage),
    stage = case_when(
      stage %in% c("1") ~ "StageI",
      stage %in% c("2") ~ "StageII",
      stage %in% c("3") ~ "StageIII",
      stage %in% c("4") ~ "StageIV",
      TRUE ~ stage
    )
  )

data_clean <- data_clean %>%
  mutate(
    gender = as.character(gender),
    gender = case_when(
      gender %in% c("0") ~ "FEMALE",
      gender %in% c("1") ~ "MALE",
      TRUE ~ gender
    )
  )

data_clean <- data_clean %>%
  mutate(
    across(
      c(T, M, N),
      ~ if_else(
        is.na(.),
        as.character(.),
        paste(cur_column(), ., sep = "")
      )
    )
  )

data_clean[data_clean == "'--"] <- NA

write.csv(
  data_clean,
  file = "1.数据准备/2.处理好的生存数据(因子型).csv",
  row.names = FALSE
)

# 30. 输出完成提示
message("TCGA 表达数据下载完成。")
message("counts / TPM / FPKM 已保存到 exp 文件夹。")
message("gene_type 统计已保存：exp/gene_type_counts.csv")
message("GeneType 筛选完成。")
message("表达矩阵去重完成。")
message("临床数据整理完成。")
message("表达数据保存目录：", folder_name)
message("表达矩阵预处理目录：", output_dir)
message("临床数据保存目录：1.数据准备")