# 1. 加载必要 R 包

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(future)
})


# 2. 检查全局环境中的 seurat 对象

if (!exists("seurat", envir = .GlobalEnv)) {
  stop("全局环境中没有 seurat 对象，请先加载。")
}

seurat_obj <- get("seurat", envir = .GlobalEnv)


# 3. 设置 CellChat 运行参数

# 3.1 输出目录
out_dir <- "细胞通讯分析"

# 3.2 物种参数
# 可选：
# species <- "human"
# species <- "mouse"
species <- "human"

# 3.3 并行 workers 数量
workers <- 1

# 3.4 filterCommunication 的 min.cells 参数
min_cells <- 10


# 4. 检查 meta.data 中必要列

# 4.1 检查 cellType 列
if (!"cellType" %in% colnames(seurat_obj@meta.data)) {
  stop("seurat@meta.data 中不存在 cellType 列")
}

# 4.2 检查 group 列
if (!"group" %in% colnames(seurat_obj@meta.data)) {
  stop("seurat@meta.data 中不存在 group 列")
}

# 4.3 检查 CellChat 运行所需的元数据列
need_cols <- c(
  "orig.ident",
  "nCount_RNA",
  "nFeature_RNA",
  "group",
  "cellType"
)

miss_cols <- setdiff(
  need_cols,
  colnames(seurat_obj@meta.data)
)

if (length(miss_cols) > 0) {
  stop(paste0(
    "seurat@meta.data 缺少以下列：",
    paste(miss_cols, collapse = ", ")
  ))
}


# 5. 整理细胞类型信息

# 5.1 将 cellType 转换为 factor，并去除未使用水平
seurat_obj$cellType <- droplevels(
  factor(seurat_obj$cellType)
)


# 6. 提取表达矩阵

# 6.1 从 Seurat V5 的 data 层提取表达矩阵
data_mat <- GetAssayData(
  object = seurat_obj,
  layer = "data"
)

# 6.2 检查表达矩阵是否成功提取
if (is.null(data_mat)) {
  stop("无法从 seurat 对象中提取 data 层表达矩阵。")
}


# 7. 整理 CellChat 所需 meta.data

# 7.1 提取需要的 meta.data 列
meta_data <- seurat_obj@meta.data[, need_cols, drop = FALSE]

# 7.2 将 group 列重命名为 samples
colnames(meta_data)[colnames(meta_data) == "group"] <- "samples"

# 7.3 将 samples 转换为 factor
meta_data$samples <- as.factor(meta_data$samples)


# 8. 创建 CellChat 对象

cellchat_obj <- createCellChat(
  object = data_mat,
  meta = meta_data,
  group.by = "cellType"
)


# 9. 加载 CellChat 数据库

# 9.1 根据物种选择数据库
if (species == "human") {
  cellchat_db <- CellChatDB.human
} else if (species == "mouse") {
  cellchat_db <- CellChatDB.mouse
} else {
  stop("species 仅支持 human 或 mouse")
}

# 9.2 设置 CellChat 数据库
cellchat_obj@DB <- cellchat_db


# 10. 预处理 CellChat 数据

# 10.1 筛选用于通讯分析的表达数据
cellchat_obj <- subsetData(cellchat_obj)


# 11. 设置并行计算

# 11.1 设置 future 并行模式
future::plan(
  strategy = "multisession",
  workers = workers
)


# 12. 识别高表达基因和高表达相互作用

# 12.1 识别高表达基因
cellchat_obj <- identifyOverExpressedGenes(cellchat_obj)

# 12.2 识别高表达配体-受体相互作用
cellchat_obj <- identifyOverExpressedInteractions(cellchat_obj)


# 13. 计算细胞通讯概率

# 13.1 计算通讯概率
cellchat_obj <- computeCommunProb(cellchat_obj)

# 13.2 按最小细胞数过滤通讯结果
cellchat_obj <- filterCommunication(
  cellchat_obj,
  min.cells = min_cells
)

# 13.3 计算信号通路水平的通讯概率
cellchat_obj <- computeCommunProbPathway(cellchat_obj)

# 13.4 汇总细胞通讯网络
cellchat_obj <- aggregateNet(cellchat_obj)


# 14. 保存 CellChat 对象到全局环境

cellchat <- cellchat_obj

assign(
  "cellchat",
  cellchat_obj,
  envir = .GlobalEnv
)


# 15. 整理运行结果信息

# 15.1 整理可展示的信号通路
pathways_txt <- if (length(cellchat_obj@netP$pathways) > 0) {
  paste(cellchat_obj@netP$pathways, collapse = ", ")
} else {
  "无"
}

# 15.2 生成运行结果文本
run_text <- paste0(
  "运行完成。
",
  "cellchat 已保存到全局环境变量：cellchat
",
  "细胞群数量：", length(levels(cellchat_obj@idents)), "
",
"可展示的信号通路：", pathways_txt
)

# 15.3 打印运行结果
cat(run_text, "\n")


# 16. 保存参数记录

# 16.1 参数文件名
name_params <- "cellchat_run_parameters"

# 16.2 创建输出目录
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# 16.3 整理参数文本
param_text <- paste0(
  "本次 CellChat 运行参数总结：
",
  "- species：", species, "
",
"- workers：", workers, "
",
"- min.cells：", min_cells, "
",
"- cellchat 对象已保存到全局环境：cellchat"
)

# 16.4 保存参数文件
param_file <- file.path(
  out_dir,
  paste0(name_params, ".txt")
)

writeLines(
  param_text,
  con = param_file
)

cat("参数保存成功：", param_file, "\n")


# 17. 恢复 future 为顺序运行

future::plan(future::sequential)