suppressPackageStartupMessages({
  library(Seurat)
})

# 1. 读取Spatial_Data并设置比较方式

# 1.1 参数设置

compare_type <- "区域间比较"
group_by_value <- "cellType"
out_dir <- "5.1细胞间差异分析"

# 1.2 检查并读取Spatial_Data对象

if (!exists("Spatial_Data", envir = .GlobalEnv)) {
  stop("全局环境中没有 Spatial_Data 对象。")
}

srt <- get("Spatial_Data", envir = .GlobalEnv)

if (!inherits(srt, "Seurat")) {
  stop("Spatial_Data 不是有效的 Seurat 对象。")
}

# 1.3 检查分组字段

if (!group_by_value %in% colnames(srt@meta.data)) {
  stop(
    paste0(
      "Spatial_Data@meta.data 中不存在列：",
      group_by_value
    )
  )
}

# 2. 运行区域间差异表达分析

# 2.1 参数设置

only_pos <- FALSE
min_pct <- 0.01
logfc_threshold <- 0.1

# 2.2 运行FindAllMarkers

all_markers <- FindAllMarkers(
  srt,
  group.by = group_by_value,
  only.pos = only_pos,
  min.pct = min_pct,
  logfc.threshold = logfc_threshold
)

# 3. 保存差异基因结果

# 3.1 参数设置

out_file_name <- "所有的差异基因.csv"

# 3.2 创建输出文件夹

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# 3.3 保存全部差异基因

out_file <- file.path(out_dir, out_file_name)

write.csv(
  all_markers,
  file = out_file,
  row.names = TRUE
)

# 4. 生成分析结果说明

result_text <- paste0(
  "运行完成。\n",
  "比较方式：", compare_type, "\n",
  "group.by：", group_by_value, "\n",
  "差异基因结果已保存至：", out_file, "\n",
  "结果行数：", nrow(all_markers), "\n",
  "结果列数：", ncol(all_markers)
)

# 5. 保存分析参数记录

# 5.1 参数设置

name_params <- "DEG_parameters"

# 5.2 生成参数记录

param_text <- paste0(
  "本次分析参数总结：\n",
  "- 比较方式：", compare_type, "\n",
  "- group.by：", group_by_value, "\n",
  "- only.pos：", only_pos, "\n",
  "- min.pct：", min_pct, "\n",
  "- logfc.threshold：", logfc_threshold, "\n",
  "- 输出目录：", out_dir, "\n",
  "- 输出文件：", out_file_name, "\n",
  "写作提示词（自行组装语言，或借助AI组装）：\n",
  "1. 基于指定分组变量进行差异表达分析（FindAllMarkers）。\n",
  "2. 筛选并导出所有比较组之间的差异基因结果。\n",
  "3. 后续可结合显著差异基因开展功能富集分析与生物学解释。"
)

# 5.3 保存参数记录

param_file <- file.path(
  out_dir,
  paste0(name_params, ".txt")
)

writeLines(
  param_text,
  con = param_file
)