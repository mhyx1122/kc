# 1. 加载必要 R 包

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggpubr)
  library(rlang)
  library(Seurat)
})


# 2. 检查全局环境中的 Seurat 对象

if (!exists("seurat", envir = .GlobalEnv)) {
  stop("全局环境中没有 seurat 对象，请先加载。")
}

srt <- get("seurat", envir = .GlobalEnv)


# 3. 读取 meta.data 并设置基础参数

# 3.1 基础参数
x_var <- "cellType"
y_var <- "singscore_Pyroptosis"

# 可选：
# facet_by <- "None"
# facet_by <- "group"
# facet_by <- "cluster"
facet_by <- "None"

# 3.2 读取 meta.data
md <- srt@meta.data

# 3.3 检查 X 轴变量是否存在
if (!x_var %in% colnames(md)) {
  stop(paste0("meta.data 中不存在列：", x_var))
}

# 3.4 检查分面变量是否存在
if (facet_by != "None" && !facet_by %in% colnames(md)) {
  stop(paste0("meta.data 中不存在分面列：", facet_by))
}


# 4. 检查 Y 轴变量

# 4.1 如果 Y 轴变量不在 meta.data 中，则尝试从 Seurat 的 data 层提取
if (!y_var %in% colnames(md)) {
  expr_mat <- GetAssayData(
    srt,
    layer = "data"
  )
  
  if (y_var %in% rownames(expr_mat)) {
    srt@meta.data[[y_var]] <- FetchData(
      srt,
      vars = y_var,
      layer = "data"
    )[, 1]
    
    md <- srt@meta.data
    
  } else {
    stop(paste0(
      "Y轴变量既不在 meta.data 中，也不在 data 层基因列表中：",
      y_var
    ))
  }
}

# 4.2 更新当前 seurat 对象
seurat <- srt


# 5. 设置小提琴图参数

# 5.1 是否添加散点和箱线图
add_jitter <- TRUE
add_boxplot <- TRUE

# 5.2 小提琴图参数
violin_alpha <- 0.4
violin_alpha <- 0.8
violin_width <- 0.6

# 5.3 散点图参数
jitter_width <- 0.2
jitter_size <- 0.2

# 5.4 箱线图参数
boxplot_width <- 0.2
boxplot_alpha <- 0.1

# 5.5 X 轴标签旋转角度
angle <- 90

# 5.6 颜色集合
palette_text <- "#FBB4AE,#B3CDE3,#CCEBC5,#DECBE4,#FED9A6,#FFFFCC,#E5D8BD,#FDDAEC,#F2F2F2"

color_palette <- unlist(strsplit(palette_text, ","))
color_palette <- trimws(color_palette)
color_palette <- color_palette[color_palette != ""]

if (length(color_palette) == 0) {
  stop("颜色集合不能为空。")
}

# 5.7 分面参数
split_by <- if (facet_by == "None") {
  NULL
} else {
  facet_by
}


# 6. 设置显著性比较参数

# 6.1 是否添加显著性比较
# 如果不需要显著性比较，保持 list()
comparisons_list <- list()

# 6.2 如果需要显著性比较，按下面格式手动填写
# 例如：
# comparisons_list <- list(
#   c("Fibroblasts", "Monocytes"),
#   c("Fibroblasts", "Macrophages")
# )

# 6.3 检查比较组是否合理
group_values <- unique(as.character(md[[x_var]]))
group_values <- group_values[!is.na(group_values) & group_values != ""]

if (length(comparisons_list) > 0) {
  comparisons_list <- comparisons_list[
    vapply(
      comparisons_list,
      function(x) {
        length(x) == 2 &&
          all(!is.null(x)) &&
          all(!is.na(x)) &&
          all(x != "") &&
          all(x %in% group_values) &&
          x[1] != x[2]
      },
      logical(1)
    )
  ]
}


# 7. 绘制小提琴图

# 7.1 准备动态变量
x_sym <- rlang::sym(x_var)
y_sym <- rlang::sym(y_var)

# 7.2 绘制基础小提琴图
p <- ggviolin(
  md,
  alpha = violin_alpha,
  width = violin_width,
  facet.by = split_by,
  x = x_var,
  y = y_var,
  color = x_var,
  add = "none",
  fill = x_var,
  add.params = list(color = "black")
) +
  scale_color_manual(values = color_palette) +
  scale_fill_manual(values = color_palette) +
  theme(
    axis.text.x = element_text(
      angle = angle,
      vjust = 0.5,
      hjust = 1
    ),
    legend.position = "none"
  ) +
  labs(x = "")

# 7.3 添加散点图
if (add_jitter) {
  p <- p +
    geom_jitter(
      aes(
        color = !!x_sym,
        y = !!y_sym,
        x = !!x_sym
      ),
      width = jitter_width,
      size = jitter_size
    )
}

# 7.4 添加箱线图
if (add_boxplot) {
  p <- p +
    geom_boxplot(
      aes(
        fill = !!x_sym,
        y = !!y_sym,
        x = !!x_sym
      ),
      color = "black",
      width = boxplot_width,
      alpha = boxplot_alpha,
      outlier.shape = NA
    )
}

# 7.5 添加显著性比较
if (length(comparisons_list) > 0) {
  p <- p +
    stat_compare_means(
      comparisons = comparisons_list,
      label = "p.signif",
      method = "t.test"
    )
}

# 7.6 显示图片
print(p)


# 8. 保存小提琴图

# 8.1 保存参数
out_dir <- "13.MHviolin可视化"
plot_name <- "violin_plot"
img_width <- 12
img_height <- 9

# 8.2 创建输出文件夹
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# 8.3 保存 PDF
plot_file <- file.path(out_dir, paste0(plot_name, ".pdf"))

ggsave(
  filename = plot_file,
  plot = p,
  device = "pdf",
  width = img_width,
  height = img_height
)

cat("图片保存成功：", plot_file, "\n")


# 9. 保存参数记录

# 9.1 参数文件名
name_params <- "MHviolin_parameters"

# 9.2 整理显著性比较信息
comparisons_text <- if (length(comparisons_list) == 0) {
  "无"
} else {
  paste(
    vapply(
      comparisons_list,
      function(x) paste0(x[1], " vs ", x[2]),
      character(1)
    ),
    collapse = " ; "
  )
}

# 9.3 生成参数文本
param_text <- paste0(
  "本次分析参数总结：
",
  "- 数据来源：seurat@meta.data
",
  "- X轴变量：", x_var, "
",
"- Y轴变量：", y_var, "
",
"- 若Y轴变量原本不在meta.data中，则已从data层提取
",
"- 分面变量：", facet_by, "
",
"- 添加散点图：", add_jitter, "
",
"- 添加箱线图：", add_boxplot, "
",
"- 小提琴透明度：", violin_alpha, "
",
"- 散点宽度：", jitter_width, "
",
"- 散点大小：", jitter_size, "
",
"- 箱线图宽度：", boxplot_width, "
",
"- 箱线图透明度：", boxplot_alpha, "
",
"- X轴标签旋转角度：", angle, "
",
"- 显著性比较数量：", length(comparisons_list), "
",
"- 比较组：", comparisons_text, "
",
"- 颜色集合：", palette_text, "
"
)

# 9.4 保存参数文件
param_file <- file.path(out_dir, paste0(name_params, ".txt"))

writeLines(
  param_text,
  con = param_file
)

cat("参数保存成功：", param_file, "\n")