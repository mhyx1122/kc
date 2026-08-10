# 1. 加载R包

suppressPackageStartupMessages({
  library(forestploter)
})


# 2. 选择作图类型并读取数据

# 可设置为"exposure"或"outcomes"
analysis_type <- "exposure"

# exposure一般读取MRPRESSO_result_exposure.csv
# outcomes一般读取MRPRESSO_result_outcome.csv
file_path <- "MRPRESSO_result_exposure.csv"

FINPRESSO_results <- read.csv(
  file_path,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

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


# 3. 检查输入数据

required_cols <- c(
  id_col, "nSNP", "MR.PRESSO", "P3", "or", "or_lci95", "or_uci95",
  "heterogeneity1", "heterogeneity2", "pleiotropy"
)

missing_cols <- setdiff(required_cols, colnames(FINPRESSO_results))

if (length(missing_cols) > 0) {
  stop("输入文件缺少以下必要列：", paste(missing_cols, collapse = "、"))
}


# 4. 将作图相关列转换为数值型

numeric_cols <- c(
  "nSNP", "MR.PRESSO", "P3", "or", "or_lci95", "or_uci95",
  "heterogeneity1", "heterogeneity2", "pleiotropy"
)

for (col_name in numeric_cols) {
  FINPRESSO_results[[col_name]] <- suppressWarnings(as.numeric(FINPRESSO_results[[col_name]]))
}


# 5. 筛选能够用于森林图的数据

# 异质性和多效性P值允许缺失，缺失时在森林图中显示为空白
plot_required_cols <- c(id_col, "nSNP", "MR.PRESSO", "P3", "or", "or_lci95", "or_uci95")

keep_rows <- complete.cases(FINPRESSO_results[, plot_required_cols, drop = FALSE])
keep_rows <- keep_rows & !is.na(FINPRESSO_results[[id_col]])
keep_rows <- keep_rows & trimws(as.character(FINPRESSO_results[[id_col]])) != ""

final_results_clean <- FINPRESSO_results[keep_rows, , drop = FALSE]

if (nrow(final_results_clean) == 0) {
  stop("作图必需列去除缺失值后，没有剩余数据可以绘图。")
}

# 按照数值型IVW方法P值从小到大排序
final_results_clean <- final_results_clean[order(final_results_clean$P3), , drop = FALSE]


# 6. 格式化MR-PRESSO Global Test P值

mrpresso_text <- ifelse(
  final_results_clean$MR.PRESSO < 0.001,
  "<0.001",
  sprintf("%.3f", final_results_clean$MR.PRESSO)
)


# 7. 格式化IVW方法P值

p3_text <- ifelse(
  final_results_clean$P3 < 0.001,
  "<0.001",
  sprintf("%.3f", final_results_clean$P3)
)


# 8. 格式化异质性检验P值

heterogeneity1_text <- ifelse(
  is.na(final_results_clean$heterogeneity1),
  "",
  ifelse(
    final_results_clean$heterogeneity1 < 0.001,
    "<0.001",
    sprintf("%.3f", final_results_clean$heterogeneity1)
  )
)

heterogeneity2_text <- ifelse(
  is.na(final_results_clean$heterogeneity2),
  "",
  ifelse(
    final_results_clean$heterogeneity2 < 0.001,
    "<0.001",
    sprintf("%.3f", final_results_clean$heterogeneity2)
  )
)

heterogeneity_text <- paste0(
  ifelse(heterogeneity1_text == "", "", paste0(heterogeneity1_text, " (MR Egger)")),
  ifelse(heterogeneity1_text != "" & heterogeneity2_text != "", "\n", ""),
  ifelse(heterogeneity2_text == "", "", paste0(heterogeneity2_text, " (IVW)"))
)


# 9. 格式化水平多效性检验P值

pleiotropy_text <- ifelse(
  is.na(final_results_clean$pleiotropy),
  "",
  ifelse(
    final_results_clean$pleiotropy < 0.001,
    "<0.001",
    sprintf("%.3f", final_results_clean$pleiotropy)
  )
)


# 10. 构建森林图数据

empty_df <- data.frame(
  analysis_id = as.character(final_results_clean[[id_col]]),
  nSNP = as.integer(final_results_clean$nSNP),
  `MR.PRESSO\nGlobal Test P-value` = mrpresso_text,
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


# 11. 构建森林图展示表格

dt <- empty_df

# 空白列用于放置森林图中的点估计和置信区间
dt$` ` <- paste(rep(" ", 25), collapse = " ")

# OR和95%置信区间仅保留3位小数，不转换成<0.001
dt$`OR(95%CI)` <- sprintf("%.3f(%.3f to %.3f)", dt$or, dt$or_lci95, dt$or_uci95)

# 展示列顺序：ID、nSNP、MR-PRESSO、IVW P值、森林图、OR、异质性、多效性
forest_table <- dt[, c(1:4, 10, 11, 8:9), drop = FALSE]

for (col_name in colnames(forest_table)) {
  forest_table[[col_name]] <- as.character(forest_table[[col_name]])
  forest_table[[col_name]][is.na(forest_table[[col_name]])] <- " "
}


# 12. 设置置信区间样式

base_size <- 10
ci_pch <- 20
ci_col <- "#4575b4"
ci_lty <- 1
ci_lwd <- 2.3
ci_Theight <- 0.2


# 13. 设置参考线样式

refline_lwd <- 1.5
refline_lty <- 2
refline_col <- "red"


# 14. 设置汇总菱形样式

summary_fill <- "#4575b4"
summary_col <- "#4575b4"


# 15. 设置脚注样式

footnote_cex <- 1
footnote_fontface <- "italic"
footnote_col <- "blue"


# 16. 创建森林图主题

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


# 17. 设置森林图坐标轴参数

# 点估计大小
sizes <- 0.6

# X轴范围
xlim <- c(0, 2)

# X轴刻度
ticks_at <- c(0, 0.5, 1, 1.5, 2)

# X轴左右标签
arrow_lab <- c("protective factor", "risk factor")

# 脚注文字
footnote <- "P<0.05 was considered statistically significant"

if (length(xlim) != 2 || xlim[1] >= xlim[2]) {
  stop("xlim必须包含两个数值，并且最小值必须小于最大值。")
}

if (length(arrow_lab) != 2) {
  stop("arrow_lab必须包含左侧和右侧两个标签。")
}


# 18. 绘制森林图

forest_plot <- forest(
  forest_table,
  est = dt$or,
  lower = dt$or_lci95,
  upper = dt$or_uci95,
  sizes = sizes,
  ci_column = 5,
  ref_line = 1,
  xlim = xlim,
  ticks_at = ticks_at,
  arrow_lab = arrow_lab,
  footnote = footnote,
  theme = forest_theme_setting
)


# 19. 保存森林图

output_file <- "forest森林图带MRPRESSO.pdf"
width <- 15
height <- 12

pdf(output_file, width = width, height = height)
print(forest_plot)
dev.off()

message("森林图绘制完成，结果已保存至：", output_file)