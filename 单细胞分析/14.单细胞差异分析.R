suppressPackageStartupMessages({
  library(Seurat)
})

# 1. 检查全局环境中是否存在 seurat 对象

if (!exists("seurat", envir = .GlobalEnv)) {
  stop("全局环境中没有 seurat 对象")
}

srt <- get("seurat", envir = .GlobalEnv)


# 2. 差异分析参数模块

# 2.1 比较方式参数
# 可选："细胞间比较" 或 "分组间比较"
compare_type <- "细胞间比较"

# 2.2 FindAllMarkers 参数
min_pct <- 0.01
logfc_threshold <- 0.1
only_pos <- FALSE

# 2.3 根据比较方式设置分组字段和输出目录
if (compare_type == "细胞间比较") {
  out_dir <- "5.1细胞间差异分析"
  group_by_value <- "cellType"
} else if (compare_type == "分组间比较") {
  out_dir <- "5.2分组间差异分析"
  group_by_value <- "group"
} else {
  stop("compare_type 只能设置为：细胞间比较 或 分组间比较")
}

# 2.4 创建输出目录
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}


# 3. 参数检查模块

# 3.1 检查分组字段是否存在
if (!group_by_value %in% colnames(srt@meta.data)) {
  stop(paste0("seurat@meta.data 中不存在列：", group_by_value))
}


# 4. 差异分析模块

# 4.1 运行 FindAllMarkers
all_markers <- FindAllMarkers(
  srt,
  group.by = group_by_value,
  only.pos = only_pos,
  min.pct = min_pct,
  logfc.threshold = logfc_threshold
)

# 4.2 保存差异基因结果
out_file <- file.path(out_dir, "所有的差异基因.csv")

write.csv(
  all_markers,
  file = out_file,
  row.names = TRUE
)


# 5. 参数记录模块

# 5.1 参数文件保存参数
name_params <- "DEG_parameters"

# 5.2 生成参数记录文本
param_text <- paste0(
  "本次分析参数总结：\n",
  "- 比较方式：", compare_type, "\n",
  "- group.by：", group_by_value, "\n",
  "- only.pos：", only_pos, "\n",
  "- min.pct：", min_pct, "\n",
  "- logfc.threshold：", logfc_threshold, "\n",
  "- 输出目录：", out_dir, "\n",
  "- 输出文件：所有的差异基因.csv\n",
  "- 结果行数：", nrow(all_markers), "\n",
  "- 结果列数：", ncol(all_markers), "\n",
  "写作提示词（自行组装语言，或借助AI组装）：\n",
  "1.基于指定分组变量进行差异表达分析（FindAllMarkers）。\n",
  "2.筛选并导出所有比较组之间的差异基因结果。\n",
  "3.后续可结合显著差异基因开展功能富集分析与生物学解释。"
)

# 5.3 保存参数记录
writeLines(
  param_text,
  con = file.path(out_dir, paste0(name_params, ".txt"))
)