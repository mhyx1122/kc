# 1. 加载必要 R 包
library(Mime1)
library(ggplot2)
library(survival)
library(aplot)

# 2. 设置输入数据文件路径
dataset_files <- c(
  "预后的示例数据(训练集).csv",
  "预后的示例数据(验证集1).csv",
  "预后的示例数据(验证集2) .csv"
)

# 3. 设置每个数据集名称
dataset_names <- c(
  "Dataset1",
  "Dataset2",
  "Dataset3"
)

# 4. 设置训练集名称
train_data_name <- "Dataset1"

# 5. 设置候选基因文件路径
gene_file_path <- "genelist.csv"

# 6. 设置既往模型信息文件路径
pre_sig_file_path <- "过往模型收集的示例数据.csv"

# 7. 检查数据集文件数量和名称数量是否一致
if (length(dataset_files) != length(dataset_names)) {
  stop("dataset_files 和 dataset_names 的数量不一致。")
}

# 8. 读取多个数据集，生成 list_train_vali_Data
list_train_vali_Data <- list()

for (i in seq_along(dataset_files)) {
  dataset <- read.csv(dataset_files[i])
  list_train_vali_Data[[dataset_names[i]]] <- dataset
}

# 9. 检查训练集名称是否存在
if (!train_data_name %in% names(list_train_vali_Data)) {
  stop("train_data_name 不在 list_train_vali_Data 的名称中。")
}

# 10. 读取候选基因列表
gene_data <- read.csv(gene_file_path, header = FALSE)

genelist_kk <- as.vector(gene_data[[1]])

cat("纳入的基因个数:\n")
print(length(genelist_kk))
print(genelist_kk)

# 11. 设置预后机器学习参数
unicox_filter <- TRUE
unicox_p_cutoff <- 0.05
nodesizemh <- 5
seedmh <- 5201314

# 12. 提取训练集数据
train_F_data <- list_train_vali_Data[[train_data_name]]

# 13. 运行预后组合机器学习
res <- ML.Dev.Prog.Sig(
  train_data = train_F_data,
  list_train_vali_Data = list_train_vali_Data,
  unicox.filter.for.candi = unicox_filter,
  unicox_p_cutoff = unicox_p_cutoff,
  candidate_genes = genelist_kk,
  mode = "all",
  nodesize = nodesizemh,
  seed = seedmh
)

# 14. 输出 C-index 可选模型名称
cat("C-index 可选模型名称:\n")
print(res$Cindex.res$Model)

# 15. 输出 riskscore 中的模型或结果名称
cat("res$riskscore 中的名称:\n")
print(names(res$riskscore))

# 16. 设置整体 C-index 图参数
cindex_all_colors <- c(
  "#A1D99BFF",
  "#F5FACD",
  "#FD8D3CFF",
  "#79AF97",
  "#8491B4"
)

cindex_bar_width <- 0.35
cindex_all_width <- 8
cindex_all_height <- 6

# 17. 绘制并保存整体 C-index 图
pdf(
  file = paste0("Cindex_plot_All_", Sys.Date(), ".pdf"),
  width = cindex_all_width,
  height = cindex_all_height
)

p_cindex_all <- cindex_dis_all(
  res,
  validate_set = names(list_train_vali_Data)[-1],
  order = names(list_train_vali_Data),
  width = cindex_bar_width,
  color = cindex_all_colors
)

if (!is.null(p_cindex_all)) {
  print(p_cindex_all)
}

dev.off()

# 18. 设置指定模型 C-index 图参数
model_MLnameYH <- "StepCox[forward] + Enet[α=0.1]"

dataset_colors <- c(
  "#9ECAE1FF",
  "#FDAE6BFF",
  "#A1D99BFF"
)

single_cindex_width <- 8
single_cindex_height <- 6

# 19. 检查指定模型 C-index 图颜色数量是否和数据集数量一致
if (length(dataset_colors) != length(list_train_vali_Data)) {
  stop("dataset_colors 的颜色数量必须和 list_train_vali_Data 的数据集数量一致。")
}

# 20. 绘制并保存指定模型 C-index 图
pdf(
  file = paste0("Cindex_plot_single_", Sys.Date(), ".pdf"),
  width = single_cindex_width,
  height = single_cindex_height
)

p_cindex_single <- cindex_dis_select(
  res,
  model = model_MLnameYH,
  dataset_col = dataset_colors,
  order = names(list_train_vali_Data)
)

if (!is.null(p_cindex_single)) {
  print(p_cindex_single)
}

dev.off()

# 21. 设置 KM 曲线参数
km_model_name <- "GBM"

km_risk_colors <- c(
  "#868686",
  "#B24745"
)

km_numb_data <- length(list_train_vali_Data)
km_width <- 8
km_height <- 6

# 22. 检查 KM 曲线数量是否超过数据集数量
if (km_numb_data > length(list_train_vali_Data)) {
  stop("km_numb_data 不能大于 list_train_vali_Data 中的数据集数量。")
}

# 23. 绘制多个数据集的 KM 曲线
survplot <- vector("list", km_numb_data)

for (i in seq_len(km_numb_data)) {
  survplot[[i]] <- rs_sur(
    res,
    model_name = km_model_name,
    dataset = names(list_train_vali_Data)[i],
    median.line = "hv",
    color = km_risk_colors,
    cutoff = 0.5,
    conf.int = TRUE,
    xlab = "Day",
    pval.coord = c(1000, 0.9)
  )
}

# 24. 保存 KM 曲线组合图
pdf(
  file = paste0("survival_plot_", Sys.Date(), ".pdf"),
  width = km_width,
  height = km_height
)

p_km <- aplot::plot_list(
  gglist = survplot,
  ncol = 2
)

if (!is.null(p_km)) {
  print(p_km)
}

dev.off()

# 25. 设置整体 AUC 图参数
auc_all_train_dataset_name <- "Dataset1"

custom_colors_AUC <- c(
  "#A1D99BFF",
  "#F5FACD",
  "#FD8D3CFF"
)

mean_colors_AUC <- c(
  "#79AF97",
  "#8491B4",
  "grey",
  "black"
)

AUC_time_YH <- 3
auc_cal_method <- "NNE"
auc_bar_width <- 0.35
auc_all_width <- 8
auc_all_height <- 6

# 26. 检查整体 AUC 使用的数据集是否存在
if (!auc_all_train_dataset_name %in% names(list_train_vali_Data)) {
  stop("auc_all_train_dataset_name 不在 list_train_vali_Data 的名称中。")
}

# 27. 计算整体 AUC
all.auc.main <- cal_AUC_ml_res(
  res.by.ML.Dev.Prog.Sig = res,
  train_data = list_train_vali_Data[[auc_all_train_dataset_name]],
  inputmatrix.list = list_train_vali_Data,
  mode = "all",
  AUC_time = AUC_time_YH,
  auc_cal_method = auc_cal_method
)

# 28. 绘制并保存整体 AUC 图
pdf(
  file = paste0("AUC_plot_All_", Sys.Date(), ".pdf"),
  width = auc_all_width,
  height = auc_all_height
)

p_auc_all <- auc_dis_all(
  all.auc.main,
  dataset = names(list_train_vali_Data),
  validate_set = names(list_train_vali_Data)[-1],
  order = names(list_train_vali_Data),
  width = auc_bar_width,
  color = c(custom_colors_AUC, mean_colors_AUC),
  year = AUC_time_YH
)

if (!is.null(p_auc_all)) {
  print(p_auc_all)
}

dev.off()

# 29. 设置特定模型 AUC 曲线参数
model_nameAUC <- "GBM"
year_YH <- 3
AUcal_method <- "NNE"
auc_single_width <- 8
auc_single_height <- 6

# 30. 计算特定模型 AUC 曲线数据
all.auc.single <- cal_AUC_ml_res(
  res.by.ML.Dev.Prog.Sig = res,
  train_data = list_train_vali_Data[[train_data_name]],
  inputmatrix.list = list_train_vali_Data,
  mode = "all",
  AUC_time = year_YH,
  auc_cal_method = AUcal_method
)

# 31. 绘制并保存特定模型 AUC 曲线图
pdf(
  file = paste0("AUC_plot_single_", Sys.Date(), ".pdf"),
  width = auc_single_width,
  height = auc_single_height
)

p_auc_single <- roc_vis(
  all.auc.single,
  model_name = model_nameAUC,
  dataset = names(list_train_vali_Data),
  order = names(list_train_vali_Data),
  anno_position = c(0.65, 0.55),
  year = year_YH
)

if (!is.null(p_auc_single)) {
  print(p_auc_single)
}

dev.off()

# 32. 设置特定模型 1、3、5 年 AUC 图参数
model_nameAUC135 <- "StepCox[forward] + plsRcox"
AUcal_method135 <- "NNE"
auc135_width <- 8
auc135_height <- 6

# 33. 分别计算 1、3、5 年 AUC
all.auc.1y <- cal_AUC_ml_res(
  res.by.ML.Dev.Prog.Sig = res,
  train_data = list_train_vali_Data[[train_data_name]],
  inputmatrix.list = list_train_vali_Data,
  mode = "all",
  AUC_time = 1,
  auc_cal_method = AUcal_method135
)

all.auc.3y <- cal_AUC_ml_res(
  res.by.ML.Dev.Prog.Sig = res,
  train_data = list_train_vali_Data[[train_data_name]],
  inputmatrix.list = list_train_vali_Data,
  mode = "all",
  AUC_time = 3,
  auc_cal_method = AUcal_method135
)

all.auc.5y <- cal_AUC_ml_res(
  res.by.ML.Dev.Prog.Sig = res,
  train_data = list_train_vali_Data[[train_data_name]],
  inputmatrix.list = list_train_vali_Data,
  mode = "all",
  AUC_time = 5,
  auc_cal_method = AUcal_method135
)

# 34. 绘制并保存特定模型 1、3、5 年 AUC 图
pdf(
  file = paste0("1_3_5_AUC_plot_", Sys.Date(), ".pdf"),
  width = auc135_width,
  height = auc135_height
)

p_auc135 <- auc_dis_select(
  list(all.auc.1y, all.auc.3y, all.auc.5y),
  model_name = model_nameAUC135,
  dataset = names(list_train_vali_Data),
  order = names(list_train_vali_Data),
  year = c(1, 3, 5)
)

if (!is.null(p_auc135)) {
  print(p_auc135)
}

dev.off()

# 35. 设置特定模型单因素 Cox 回归 Meta 分析参数
model_namemeta <- "StepCox[forward] + plsRcox"
type_modemeta <- "categorical"
meta_width <- 8
meta_height <- 6

# 36. 计算单因素 Cox 回归结果
unicox.rs.res <- cal_unicox_ml_res(
  res.by.ML.Dev.Prog.Sig = res,
  optimal.model = model_namemeta,
  type = type_modemeta
)

# 37. 计算 Meta 分析结果
metamodel <- cal_unicox_meta_ml_res(
  input = unicox.rs.res
)

# 38. 绘制并保存 Meta 分析图
pdf(
  file = paste0("meta_plot_", Sys.Date(), ".pdf"),
  width = meta_width,
  height = meta_height
)

p_meta <- meta_unicox_vis(
  metamodel,
  dataset = names(list_train_vali_Data)
)

if (!is.null(p_meta)) {
  print(p_meta)
}

dev.off()

# 39. 设置既往模型 C-index 比较参数
model_prebi <- "StepCox[forward] + plsRcox"
pre_width <- 10
pre_height <- 14

# 40. 读取既往模型信息表
pre_sig <- read.csv(pre_sig_file_path)

# 41. 计算既往模型 C-index
cc.glioma.lgg.gbm <- cal_cindex_pre.prog.sig(
  use_your_own_collected_sig = TRUE,
  collected_sig_table = pre_sig,
  list_input_data = list_train_vali_Data
)

# 42. 绘制并保存当前最佳模型与既往模型的 C-index 比较图
pdf(
  file = paste0("pre_plot_", Sys.Date(), ".pdf"),
  width = pre_width,
  height = pre_height
)

p_pre <- cindex_comp(
  cc.glioma.lgg.gbm,
  res,
  model_name = model_prebi,
  dataset = names(list_train_vali_Data)
)

if (!is.null(p_pre)) {
  print(p_pre)
}

dev.off()

# 43. 输出完成提示
message("全部分析完成。")
message("已生成预后机器学习结果 res。")
message("已保存整体 C-index 图。")
message("已保存指定模型 C-index 图。")
message("已保存 KM 曲线图。")
message("已保存整体 AUC 图。")
message("已保存指定模型 AUC 曲线图。")
message("已保存 1、3、5 年 AUC 图。")
message("已保存单因素 Cox Meta 分析图。")
message("已保存当前模型与既往模型 C-index 比较图。")