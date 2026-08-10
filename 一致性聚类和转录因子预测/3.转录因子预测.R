# 1. 加载必要 R 包
library(TFTF)
library(UpSetR)
library(dplyr)
library(tidyr)
library(ggsci)

# 2. 设置通用数据库参数
datasets <- c(
  "hTFtarget",
  "KnockTF",
  "FIMO_JASPAR",
  "PWMEnrich_JASPAR",
  "ENCODE",
  "CHEA",
  "TRRUST",
  "GTRD",
  "ChIP_Atlas"
)

cor_DB <- c("TCGA", "GTEx")

TCGA_tissue <- "COAD"
GTEx_tissue <- "Colon"

cor_cutoff <- 0.3
FIMO_score <- 10
PWMEnrich_p <- 0.1
cut_log2FC <- 1
down_only <- TRUE

folder_path <- "TF_Target_Results"

if (!dir.exists(folder_path)) {
  dir.create(folder_path, recursive = TRUE)
}

# 3. 定义结果整理函数
pad_results_to_df <- function(result_object) {
  if (is.null(result_object)) {
    stop("预测结果为空。")
  }
  
  if (is.null(result_object$results)) {
    stop("预测结果中没有 results 字段。")
  }
  
  max_length <- max(sapply(result_object$results, length))
  
  result_padded <- lapply(result_object$results, function(x) {
    length(x) <- max_length
    return(x)
  })
  
  result_df <- as.data.frame(result_padded)
  
  return(result_df)
}

# 4. 模块一：根据靶基因预测可能的转录因子
target_gene_symbol <- "GAPDH"

TF_Gene <- TFTF::predict_TF(
  datasets = datasets,
  target = target_gene_symbol,
  TCGA_tissue = TCGA_tissue,
  GTEx_tissue = GTEx_tissue,
  cor_DB = cor_DB,
  cor_cutoff = cor_cutoff,
  FIMO.score = FIMO_score,
  PWMEnrich.p = PWMEnrich_p,
  cut.log2FC = cut_log2FC,
  down.only = down_only,
  app = FALSE
)

TF_GeneNAME <- pad_results_to_df(TF_Gene)

tf_save_path <- file.path(
  folder_path,
  paste0("预测的转录因子基因（", target_gene_symbol, "）.csv")
)

write.csv(
  TF_GeneNAME,
  tf_save_path,
  row.names = FALSE
)

message("靶基因对应转录因子预测完成：", tf_save_path)

# 5. 模块二：根据转录因子预测可能的靶基因
tf_symbol <- "STAT3"

target_gene <- TFTF::predict_target(
  datasets = datasets,
  tf = tf_symbol,
  TCGA_tissue = TCGA_tissue,
  GTEx_tissue = GTEx_tissue,
  cor_DB = cor_DB,
  cor_cutoff = cor_cutoff,
  FIMO.score = FIMO_score,
  PWMEnrich.p = PWMEnrich_p,
  cut.log2FC = cut_log2FC,
  down.only = down_only,
  app = FALSE
)

Target_GeneNAME <- pad_results_to_df(target_gene)

target_save_path <- file.path(
  folder_path,
  paste0("预测的靶基因（", tf_symbol, "）.csv")
)

write.csv(
  Target_GeneNAME,
  target_save_path,
  row.names = FALSE
)

message("转录因子对应靶基因预测完成：", target_save_path)

# 6. 模块三：设置 UpSet 图输入和参数
upset_input_file <- target_save_path

output_dir <- "TF_Target_Results"

if (output_dir == "") {
  output_dir <- "5.基因取交集"
}

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

color_point <- "#e76c54"

upset_pdf_width <- 8
upset_pdf_height <- 6

# 7. 读取 UpSet 输入文件
mcc_data <- read.csv(
  upset_input_file,
  header = TRUE,
  check.names = FALSE
)

# 8. 整理成长格式数据
all_elements <- data.frame(
  Element = character(),
  Hub = character(),
  stringsAsFactors = FALSE
)

for (hub in names(mcc_data)) {
  column_data <- na.omit(mcc_data[[hub]])
  column_data <- column_data[column_data != ""]
  
  temp <- data.frame(
    Element = column_data,
    Hub = hub,
    stringsAsFactors = FALSE
  )
  
  all_elements <- rbind(all_elements, temp)
}

all_elements <- na.omit(all_elements)
all_elements <- unique(all_elements)

# 9. 转换为 UpSetR 需要的 0/1 矩阵
all_elements_upset <- all_elements %>%
  mutate(incidence = 1) %>%
  tidyr::spread(Hub, incidence, fill = 0)

# 10. 生成交集颜色
interCsecGtions <- apply(
  all_elements_upset[, -1, drop = FALSE],
  1,
  function(x) paste(which(x == 1), collapse = ",")
)

unique_intersections <- unique(interCsecGtions)

randomColor123 <- function() {
  paste0(
    "#",
    paste0(
      sample(c(0:9, letters[1:6]), 6, replace = TRUE),
      collapse = ""
    )
  )
}

intersection_colors <- replicate(
  length(unique_intersections),
  randomColor123()
)

# 11. 设置集合颜色
set_colors <- c(
  ggsci::pal_npg()(9),
  ggsci::pal_jco()(9),
  ggsci::pal_jama()(7),
  ggsci::pal_nejm()(8)
)

set_colors <- rep(
  set_colors,
  length.out = length(names(mcc_data))
)

# 12. 保存 UpSet 图
upset_pdf_path <- file.path(
  output_dir,
  "upset_plot_colored_elements.pdf"
)

upset_plot <- UpSetR::upset(
  all_elements_upset,
  sets = names(mcc_data),
  order.by = "degree",
  matrix.color = color_point,
  main.bar.color = intersection_colors,
  sets.bar.color = set_colors,
  sets.x.label = "Set Size",
  point.size = 3,
  line.size = 1,
  mb.ratio = c(0.6, 0.4),
  shade.color = "gray80",
  shade.alpha = 0.3,
  matrix.dot.alpha = 0.6,
  show.numbers = "yes",
  number.angles = 0,
  group.by = "degree",
  text.scale = 1.2,
  set_size.angles = 0
)

pdf(
  upset_pdf_path,
  width = upset_pdf_width,
  height = upset_pdf_height
)

print(upset_plot)

dev.off()

# 13. 保存 UpSet 矩阵结果
write.csv(
  all_elements_upset,
  file.path(output_dir, "all_elements_upset.csv"),
  row.names = FALSE
)

int <- all_elements_upset
rownames(int) <- int[[1]]
int$row_sums <- rowSums(int[, -1, drop = FALSE])

write.csv(
  int,
  file.path(output_dir, "all__upset.csv"),
  row.names = FALSE
)

# 14. 输出完成提示
message("全部分析完成。")
message("预测转录因子结果：", tf_save_path)
message("预测靶基因结果：", target_save_path)
message("UpSet 图：", upset_pdf_path)
message("UpSet 矩阵：", file.path(output_dir, "all_elements_upset.csv"))