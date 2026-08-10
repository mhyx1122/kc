# 1. 加载必要 R 包
library(survival)
library(survminer)
library(gridExtra)

# 2. 设置输入文件
survival_file <- "TCGA-COAD.survival.csv"
cluster_file <- "cluster_output(K=4).csv"

# 3. 设置输出参数
output_dir <- "ConsensusCluster"
cluster_pdf_width <- 8
cluster_pdf_height <- 8

# 4. 设置是否为 TCGA 数据
is_tcga_data <- TRUE

# 5. 创建输出文件夹
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# 6. 读取生存数据和聚类结果
datacluster1 <- read.csv(
  survival_file,
  header = TRUE,
  check.names = FALSE
)

datacluster2 <- read.csv(
  cluster_file,
  header = TRUE,
  check.names = FALSE
)

# 7. 统一第一列列名为 ID
colnames(datacluster1)[1] <- "ID"
colnames(datacluster2)[1] <- "ID"

# 8. 检查必要列
if (!all(c("ID", "OS.time", "OS") %in% colnames(datacluster1))) {
  stop("生存数据必须包含 ID、OS.time、OS 三列。")
}

if (!all(c("ID", "Cluster") %in% colnames(datacluster2))) {
  stop("聚类结果必须包含 ID、Cluster 两列。")
}

# 9. 去除重复样本
datacluster1 <- datacluster1[!duplicated(datacluster1$ID), , drop = FALSE]
datacluster2 <- datacluster2[!duplicated(datacluster2$ID), , drop = FALSE]

# 10. 按数据类型处理样本 ID
if (is_tcga_data) {
  datacluster1$ID <- substr(datacluster1$ID, 1, 12)
  datacluster2$ID <- substr(datacluster2$ID, 1, 12)
  
  datacluster1$ID <- gsub("-", ".", datacluster1$ID)
  datacluster2$ID <- gsub("-", ".", datacluster2$ID)
}

# 11. 合并生存数据和聚类结果
merged_data <- merge(
  datacluster1,
  datacluster2,
  by = "ID",
  all = TRUE
)

cleaned_data <- na.omit(merged_data)

if (nrow(cleaned_data) == 0) {
  stop("生存数据和聚类结果没有成功匹配到样本，请检查两个文件的样本 ID。")
}

# 12. 整理分组变量
data <- cleaned_data
data$group <- data$Cluster

# 13. 构建生存对象
surv_object <- survival::Surv(
  time = data$OS.time,
  event = data$OS
)

# 14. 拟合 KM 生存曲线
fit <- survival::survfit(
  surv_object ~ group,
  data = data
)

# 15. 进行生存曲线差异检验
diff_test <- survival::survdiff(
  surv_object ~ group,
  data = data
)

p_value <- 1 - pchisq(
  diff_test$chisq,
  length(diff_test$n) - 1
)

print(paste0("KM log-rank P value: ", signif(p_value, 4)))

# 16. 绘制 KM 生存曲线
km_plot <- survminer::ggsurvplot(
  fit,
  data = data,
  pval = TRUE,
  pval.method = TRUE,
  risk.table = TRUE,
  conf.int = TRUE,
  conf.int.style = "ribbon"
)

# 17. 保存 KM 生存曲线
pdf_path <- file.path(
  output_dir,
  "KM_survival_Cluster.pdf"
)

pdf(
  pdf_path,
  width = cluster_pdf_width,
  height = cluster_pdf_height
)

gridExtra::grid.arrange(
  km_plot$plot,
  km_plot$table,
  nrow = 2,
  heights = c(3, 1)
)

dev.off()

# 18. 输出完成提示
message("一致性聚类 KM 生存曲线分析完成。")
message("合并后的有效样本数：", nrow(cleaned_data))
message("结果文件已保存：", pdf_path)