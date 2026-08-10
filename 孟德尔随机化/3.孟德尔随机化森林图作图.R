# 1. 加载R包

suppressPackageStartupMessages({
  library(forestploter)
})


# 2. 选择分析类型并读取汇总结果

# 可设置为"exposure"或"outcomes"
analysis_type <- "exposure"

# exposure对应final_results(exposure).csv
# outcomes对应final_results(outcome).csv
file_path <- "final_results(exposure).csv"

final_results <- read.csv(file_path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

if (!analysis_type %in% c("exposure", "outcomes")) {
  stop("analysis_type只能设置为'exposure'或'outcomes'。")
}

if (analysis_type == "exposure") {
  id_col <- "id.exposure"
  id_title <- "exposure"
} else {
  id_col <- "id.outcome"
  id_title <- "outcome"
}


# 3. 检查输入文件并整理数值列

required_cols <- c(
  id_col, "nSNP", "P3", "or", "or_lci95", "or_uci95",
  "heterogeneity1", "heterogeneity2", "pleiotropy"
)

missing_cols <- setdiff(required_cols, colnames(final_results))

if (length(missing_cols) > 0) {
  stop("输入文件缺少以下必要列：", paste(missing_cols, collapse = "、"))
}

numeric_cols <- c(
  "nSNP", "P3", "or", "or_lci95", "or_uci95",
  "heterogeneity1", "heterogeneity2", "pleiotropy"
)

for (col_name in numeric_cols) {
  final_results[[col_name]] <- suppressWarnings(as.numeric(final_results[[col_name]]))
}


# 4. 筛选可以用于森林图的数据

# 只根据森林图必须使用的列删除缺失行
plot_required_cols <- c(id_col, "nSNP", "P3", "or", "or_lci95", "or_uci95")
keep_rows <- complete.cases(final_results[, plot_required_cols])
keep_rows <- keep_rows & !is.na(final_results[[id_col]]) & trimws(final_results[[id_col]]) != ""

final_results_clean <- final_results[keep_rows, , drop = FALSE]

if (nrow(final_results_clean) == 0) {
  stop("绘图必要列去除缺失值后，没有剩余数据可以绘图。")
}

# 按照数值型IVW方法P值从小到大排序
final_results_clean <- final_results_clean[order(final_results_clean$P3), , drop = FALSE]


# 5. 格式化P值显示文字

# P3小于0.001时显示为<0.001，否则保留3位小数
p3_text <- ifelse(
  final_results_clean$P3 < 0.001,
  "<0.001",
  sprintf("%.3f", final_results_clean$P3)
)

# MR-Egger异质性检验P值
heterogeneity1_text <- ifelse(
  is.na(final_results_clean$heterogeneity1),
  "",
  ifelse(
    final_results_clean$heterogeneity1 < 0.001,
    "<0.001",
    sprintf("%.3f", final_results_clean$heterogeneity1)
  )
)

# IVW异质性检验P值
heterogeneity2_text <- ifelse(
  is.na(final_results_clean$heterogeneity2),
  "",
  ifelse(
    final_results_clean$heterogeneity2 < 0.001,
    "<0.001",
    sprintf("%.3f", final_results_clean$heterogeneity2)
  )
)

# 水平多效性检验P值
pleiotropy_text <- ifelse(
  is.na(final_results_clean$pleiotropy),
  "",
  ifelse(
    final_results_clean$pleiotropy < 0.001,
    "<0.001",
    sprintf("%.3f", final_results_clean$pleiotropy)
  )
)

heterogeneity_text <- paste0(
  heterogeneity1_text, " (MR Egger)\n",
  heterogeneity2_text, " (IVW)"
)

heterogeneity_text[
  heterogeneity1_text == "" & heterogeneity2_text == ""
] <- ""


# 6. 构建森林图数据

empty_df <- data.frame(
  analysis_id = as.character(final_results_clean[[id_col]]),
  nSNP = as.integer(final_results_clean$nSNP),
  P.value = p3_text,
  or = final_results_clean$or,
  or_lci95 = final_results_clean$or_lci95,
  or_uci95 = final_results_clean$or_uci95,
  Heterogeneity.Q_pval = heterogeneity_text,
  Pleiotropy.Pval = pleiotropy_text,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

colnames(empty_df)[1] <- id_title


# 7. 构建森林图展示表格

dt <- empty_df

# 空白列用于放置森林图中的点估计和置信区间
dt$` ` <- paste(rep(" ", 25), collapse = " ")

# OR及95%置信区间只保留3位小数，不转换成<0.001
dt$`OR(95%CI)` <- sprintf(
  "%.3f(%.3f to %.3f)",
  dt$or,
  dt$or_lci95,
  dt$or_uci95
)

# 展示列顺序：ID、SNP数量、P值、森林图、OR及置信区间、异质性、多效性
forest_table <- dt[, c(1:3, 9, 10, 7:8), drop = FALSE]
forest_table[is.na(forest_table)] <- " "


# 8. 设置置信区间样式

base_size <- 10
ci_pch <- 20
ci_col <- "#4575b4"
ci_lty <- 1
ci_lwd <- 2.3
ci_Theight <- 0.2


# 9. 设置参考线样式

refline_lwd <- 1.5
refline_lty <- 2
refline_col <- "red"


# 10. 设置汇总菱形样式

summary_fill <- "#4575b4"
summary_col <- "#4575b4"


# 11. 设置脚注样式

footnote_cex <- 1
footnote_fontface <- "italic"
footnote_col <- "blue"


# 12. 创建森林图主题

forest_theme_setting <- forest_theme(
  base_size = base_size,
  ci_pch = ci_pch,
  ci_col = ci_col,
  ci_lty = ci_lty,
  ci_lwd = ci_lwd,
  ci_Theight = ci_Theight,
  refline_lwd = refline_lwd,
  refline_lty = refline_lty,
  refline_col = refline_col,
  summary_fill = summary_fill,
  summary_col = summary_col,
  footnote_cex = footnote_cex,
  footnote_fontface = footnote_fontface,
  footnote_col = footnote_col
)


# 13. 设置森林图坐标轴参数

# 点估计大小
sizes <- 0.6

# X轴范围
xlim <- c(0, 2)

# X轴刻度
ticks_at <- c(0, 0.5, 1, 1.5, 2)

# X轴左右方向标签
arrow_lab <- c("protective factor", "risk factor")

# 脚注文字
footnote <- "P<0.05 was considered statistically significant"

if (length(xlim) != 2 || xlim[1] >= xlim[2]) {
  stop("xlim必须包含两个数值，并且最小值必须小于最大值。")
}

if (length(arrow_lab) != 2) {
  stop("arrow_lab必须包含左侧和右侧两个标签。")
}


# 14. 绘制森林图

forest_plot <- forest(
  forest_table,
  est = dt$or,
  lower = dt$or_lci95,
  upper = dt$or_uci95,
  sizes = sizes,
  ci_column = 4,
  ref_line = 1,
  xlim = xlim,
  ticks_at = ticks_at,
  arrow_lab = arrow_lab,
  footnote = footnote,
  theme = forest_theme_setting
)


# 15. 保存森林图

output_file <- "forest森林图.pdf"
width <- 15
height <- 12

pdf(output_file, width = width, height = height)
print(forest_plot)
dev.off()

message("森林图绘制完成，结果已保存至：", output_file)