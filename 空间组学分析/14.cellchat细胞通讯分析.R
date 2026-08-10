suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(future)
  library(jsonlite)
})

# 1. 读取Spatial_Data并准备表达矩阵

# 1.1 参数设置

out_dir <- "细胞通讯分析"

# 使用Spatial_Data当前DefaultAssay中的data层
assay_use <- DefaultAssay(Spatial_Data)
layer_use <- "data"

# 1.2 创建输出文件夹

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# 1.3 检查Spatial_Data对象

if (!exists("Spatial_Data", envir = .GlobalEnv)) {
  stop("全局环境中没有 Spatial_Data 对象。")
}

Spatial_Data <- get("Spatial_Data", envir = .GlobalEnv)

if (!inherits(Spatial_Data, "Seurat")) {
  stop("Spatial_Data 不是有效的Seurat对象。")
}

if (length(Spatial_Data@images) == 0) {
  stop("Spatial_Data中不存在空间图像信息，无法进行空间CellChat分析。")
}

# 1.4 提取表达矩阵

expression_data <- GetAssayData(
  object = Spatial_Data,
  assay = assay_use,
  layer = layer_use
)

if (nrow(expression_data) == 0 || ncol(expression_data) == 0) {
  stop("从Spatial_Data中提取的表达矩阵为空。")
}

# 2. 准备细胞分组信息

# 2.1 使用当前Idents作为cellType

cell_idents <- Idents(Spatial_Data)

meta <- data.frame(
  cellType = as.character(cell_idents),
  row.names = names(cell_idents),
  stringsAsFactors = FALSE
)

if (any(is.na(meta$cellType)) || any(meta$cellType == "")) {
  stop("Spatial_Data的Idents中存在NA或空值，请先完成细胞类型注释。")
}

meta$cellType <- factor(meta$cellType)

cat(
  "输入的细胞类型如下：",
  paste(levels(meta$cellType), collapse = " "),
  "\n"
)

# 3. 读取空间坐标和空间尺度因子

# 3.1 参数设置

# 当Spatial_Data包含多个空间切片时，可手动指定具体切片名称
# 设置为NULL时默认使用第一个切片
image_name <- NULL

scalefactor_path <- "10X_Spatial/spatial/scalefactors_json.json"
spot_size <- 65

# 3.2 确定使用的空间切片

image_names <- names(Spatial_Data@images)

if (is.null(image_name) || length(image_name) == 0 || trimws(image_name) == "") {
  image_name <- image_names[1]
}

if (!image_name %in% image_names) {
  stop(
    paste0(
      "指定的image_name不存在。当前可选空间切片为：",
      paste(image_names, collapse = ", ")
    )
  )
}

# 3.3 提取空间坐标

spatial_locs <- Seurat::GetTissueCoordinates(
  Spatial_Data,
  image = image_name
)

if (is.null(rownames(spatial_locs))) {
  stop("提取的空间坐标没有spot名称，无法与表达矩阵匹配。")
}

if (ncol(spatial_locs) < 2) {
  stop("空间坐标数据少于两列，无法构建二维空间位置矩阵。")
}

spatial_location <- as.matrix(
  spatial_locs[, 1:2, drop = FALSE]
)

# 3.4 读取空间尺度文件

if (!file.exists(scalefactor_path)) {
  stop(
    paste0(
      "未找到空间尺度文件：",
      scalefactor_path
    )
  )
}

scalefactors <- jsonlite::fromJSON(
  txt = scalefactor_path
)

if (
  is.null(scalefactors$spot_diameter_fullres) ||
  length(scalefactors$spot_diameter_fullres) != 1 ||
  is.na(scalefactors$spot_diameter_fullres) ||
  scalefactors$spot_diameter_fullres <= 0
) {
  stop("scalefactors_json.json中不存在有效的spot_diameter_fullres。")
}

conversion_factor <- spot_size / scalefactors$spot_diameter_fullres

spatial_factors <- data.frame(
  ratio = conversion_factor,
  tol = spot_size / 2
)

# 4. 统一表达矩阵、分组信息和空间坐标的spot顺序

# 4.1 获取三类数据共有的spot

common_cells <- colnames(expression_data)
common_cells <- common_cells[common_cells %in% rownames(meta)]
common_cells <- common_cells[common_cells %in% rownames(spatial_location)]

if (length(common_cells) == 0) {
  stop("表达矩阵、meta和空间坐标之间没有共同的spot名称。")
}

# 4.2 显示未匹配spot数量

removed_expression_cells <- setdiff(
  colnames(expression_data),
  common_cells
)

removed_meta_cells <- setdiff(
  rownames(meta),
  common_cells
)

removed_coordinate_cells <- setdiff(
  rownames(spatial_location),
  common_cells
)

if (
  length(removed_expression_cells) > 0 ||
  length(removed_meta_cells) > 0 ||
  length(removed_coordinate_cells) > 0
) {
  warning(
    paste0(
      "部分spot未能在表达矩阵、meta和空间坐标中同时匹配。",
      "表达矩阵移除：", length(removed_expression_cells), "个；",
      "meta移除：", length(removed_meta_cells), "个；",
      "空间坐标移除：", length(removed_coordinate_cells), "个。"
    )
  )
}

# 4.3 按同一顺序整理三类数据

expression_data <- expression_data[
  ,
  common_cells,
  drop = FALSE
]

meta <- meta[
  common_cells,
  ,
  drop = FALSE
]

spatial_location <- spatial_location[
  common_cells,
  ,
  drop = FALSE
]

if (
  !identical(colnames(expression_data), rownames(meta)) ||
  !identical(colnames(expression_data), rownames(spatial_location))
) {
  stop("表达矩阵、meta和空间坐标的spot顺序未能统一。")
}

# 5. 创建空间CellChat对象并设置数据库

# 5.1 参数设置

# 可选：human、mouse
species <- "human"

# 5.2 创建空间CellChat对象

cellchat <- createCellChat(
  object = expression_data,
  meta = meta,
  group.by = "cellType",
  datatype = "spatial",
  coordinates = spatial_location,
  spatial.factors = spatial_factors
)

# 5.3 根据物种选择CellChat数据库

if (species == "human") {
  cellchat_db <- CellChatDB.human
} else if (species == "mouse") {
  cellchat_db <- CellChatDB.mouse
} else {
  stop("species仅支持human或mouse。")
}

cellchat@DB <- cellchat_db

# 5.4 提取CellChat数据库相关表达数据

cellchat <- subsetData(cellchat)

# 6. 识别高表达基因和配体受体互作

# 6.1 参数设置

workers <- 1

if (workers < 1) {
  stop("workers必须大于或等于1。")
}

# 6.2 设置并行并运行高表达分析

previous_future_plan <- future::plan()

tryCatch(
  {
    future::plan(
      future::multisession,
      workers = workers
    )
    
    cellchat <- identifyOverExpressedGenes(
      cellchat
    )
    
    cellchat <- identifyOverExpressedInteractions(
      cellchat
    )
  },
  finally = {
    future::plan(previous_future_plan)
  }
)

# 7. 计算空间细胞通讯概率

# 7.1 参数设置

# 可选：triMean、truncatedMean、median
compute_type <- "truncatedMean"

trim <- 0.1
distance_use <- TRUE
interaction_range <- 250
scale_distance <- 0.01
contact_dependent <- TRUE
contact_range <- 100
min_cells <- 10

if (!compute_type %in% c("triMean", "truncatedMean", "median")) {
  stop("compute_type只能设置为triMean、truncatedMean或median。")
}

if (trim < 0 || trim > 0.5) {
  stop("trim必须设置在0至0.5之间。")
}

# 7.2 计算细胞通讯概率

cellchat <- computeCommunProb(
  cellchat,
  type = compute_type,
  trim = trim,
  distance.use = distance_use,
  interaction.range = interaction_range,
  scale.distance = scale_distance,
  contact.dependent = contact_dependent,
  contact.range = contact_range
)

# 7.3 过滤细胞数量过少的通讯

cellchat <- filterCommunication(
  cellchat,
  min.cells = min_cells
)

# 7.4 计算信号通路水平的通讯概率

cellchat <- computeCommunProbPathway(
  cellchat
)

# 7.5 汇总细胞通讯网络

cellchat <- aggregateNet(
  cellchat
)

# 8. 保存CellChat对象到全局环境

assign(
  "cellchat",
  cellchat,
  envir = .GlobalEnv
)

pathways_text <- if (length(cellchat@netP$pathways) > 0) {
  paste(cellchat@netP$pathways, collapse = ", ")
} else {
  "无"
}

cat(
  "运行完成。\n",
  "空间CellChat对象已保存到全局环境变量：cellchat\n",
  "使用的空间切片：", image_name, "\n",
  "参与分析的spot数量：", ncol(expression_data), "\n",
  "细胞群数量：", length(levels(cellchat@idents)), "\n",
  "可展示的信号通路：", pathways_text, "\n",
  sep = ""
)

# 9. 保存运行参数

# 9.1 参数设置

name_params <- "cellchat_run_parameters"

# 9.2 生成参数记录

param_text <- paste0(
  "本次Spatial CellChat运行参数总结：\n",
  "- species：", species, "\n",
  "- assay：", assay_use, "\n",
  "- layer：", layer_use, "\n",
  "- workers：", workers, "\n",
  "- min.cells：", min_cells, "\n",
  "- datatype：spatial\n",
  "- image_name：", image_name, "\n",
  "- scalefactor_path：", scalefactor_path, "\n",
  "- spot.size：", spot_size, "\n",
  "- spot_diameter_fullres：", scalefactors$spot_diameter_fullres, "\n",
  "- conversion.factor：", conversion_factor, "\n",
  "- 参与分析的spot数量：", ncol(expression_data), "\n",
  "- type：", compute_type, "\n",
  "- trim：", trim, "\n",
  "- distance.use：", distance_use, "\n",
  "- interaction.range：", interaction_range, "\n",
  "- scale.distance：", scale_distance, "\n",
  "- contact.dependent：", contact_dependent, "\n",
  "- contact.range：", contact_range, "\n",
  "- cellchat对象已保存到全局环境：cellchat\n",
  "- 可展示的信号通路：", pathways_text
)

# 9.3 保存参数文件

param_file <- file.path(
  out_dir,
  paste0(name_params, ".txt")
)

writeLines(
  param_text,
  con = param_file
)