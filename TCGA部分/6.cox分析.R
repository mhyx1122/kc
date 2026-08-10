# 1. 加载必要 R 包

library(data.table)
library(dplyr)
library(stringr)
library(survival)
library(grid)
library(forestploter)


# 2. 设置输入文件参数

input_file <- "合并后的标准肿瘤分析数据.csv"


# 3. 设置 Cox 分析参数

gene_columns <- c(
  "stage",
  "FUT4",
  "PIWIL4",
  "CST1"
)

output_folder <- "Cox分析结果"


# 4. 设置森林图颜色和样式参数

summary_fill <- "#4575b4"

summary_col <- "#4575b4"

refline_lwd <- 1.5

refline_col <- "red"

sizes <- 0.6

xlimss <- c(0, 2)

ticks_at <- c(0, 0.5, 1, 1.5, 2)

arrow_lab <- c(
  "protective factor",
  "risk factor"
)

footnote <- "P<0.05 was considered statistically significant"

widthkkm <- 15

heightkkm <- 12


# 5. 创建结果保存文件夹

if (!dir.exists(output_folder)) {
  dir.create(output_folder)
}


# 6. 读取数据

datacoxs <- fread(
  input_file,
  data.table = FALSE
)

rownames(datacoxs) <- datacoxs[, 1]

datacoxs <- datacoxs[, -1]


# 7. 准备 Cox 分析数据

data <- datacoxs

merged_data <- data


# 8. 处理 stage 分期变量

if ("stage" %in% colnames(merged_data)) {
  
  merged_data <- merged_data %>%
    mutate(
      stage = case_when(
        stage %in% c("StageI", "StageIA", "StageIB", "StageIC") ~ "1",
        stage %in% c("StageII", "StageIIA", "StageIIB", "StageIIC") ~ "2",
        stage %in% c("StageIII", "StageIIIA", "StageIIIB", "StageIIIC") ~ "3",
        stage %in% c("StageIV", "StageIVA", "StageIVB", "StageIVC") ~ "4",
        TRUE ~ stage
      )
    )
  
} else {
  
  print("Column 'stage' is missing, skipping stage modification.")
}


# 9. 处理 T、N、M 分期变量

if (all(c("T", "N", "M") %in% colnames(merged_data))) {
  
  merged_data <- merged_data %>%
    mutate(
      T = str_extract(T, "[0-9]+"),
      N = str_extract(N, "[0-9]+"),
      M = str_extract(M, "[0-9]+")
    )
  
} else {
  
  print("One or more of the columns 'T', 'N', 'M' are missing, skipping their modification.")
}


# 10. 处理 gender 变量

if ("gender" %in% colnames(merged_data)) {
  
  merged_data <- merged_data %>%
    dplyr::mutate(
      gender = case_when(
        gender == "FEMALE" ~ 0,
        gender == "MALE" ~ 1
      )
    )
  
} else {
  
  print("Column 'gender' is missing, skipping gender modification.")
}


# 11. 将数据转换为数值型

data <- merged_data

data <- data.frame(
  lapply(
    data,
    function(x) {
      if (is.factor(x) || is.character(x)) {
        as.numeric(as.character(x))
      } else {
        as.numeric(x)
      }
    }
  )
)


# 12. 检查转换后是否均为数值型

check_numeric <- sapply(
  data,
  function(x) {
    is.numeric(na.omit(x))
  }
)

if (all(check_numeric)) {
  
  print("All columns are numeric (ignoring NA).")
  
} else {
  
  non_numeric_cols <- names(check_numeric[!check_numeric])
  
  print(
    paste(
      "The following columns are not numeric (ignoring NA):",
      paste(non_numeric_cols, collapse = ", ")
    )
  )
}


# 13. 单因素 Cox 回归分析

uniCoxResults <- data.frame()

for (i in gene_columns) {
  
  single_data <- data[, c("futime", "fustat", i), drop = FALSE]
  
  single_data <- na.omit(single_data)
  
  cat("单因素 Cox 分析变量:", i, "；纳入样本数:", nrow(single_data), "\n")
  
  if (nrow(single_data) == 0) {
    message("变量 ", i, " 去除 NA 后没有剩余样本，已跳过。")
    next
  }
  
  single_formula <- as.formula(
    paste(
      "Surv(futime, fustat) ~",
      i
    )
  )
  
  cox <- coxph(
    single_formula,
    data = single_data
  )
  
  coxSummary <- summary(cox)
  
  uniCoxResults <- rbind(
    uniCoxResults,
    data.frame(
      Variable = i,
      HR = coxSummary$conf.int[, "exp(coef)"],
      HR.95L = coxSummary$conf.int[, "lower .95"],
      HR.95H = coxSummary$conf.int[, "upper .95"],
      pvalue = coxSummary$coefficients[, "Pr(>|z|)"]
    )
  )
}

print(uniCoxResults)
uniCoxResults_path <- file.path(output_folder, "uniCoxResults.csv")
write.csv(uniCoxResults, file = uniCoxResults_path, row.names = FALSE)

# 14. 多因素 Cox 回归分析

all_vars <- colnames(data)

predictors <- intersect(
  gene_columns,
  all_vars
)

multi_data <- data[, c("futime", "fustat", predictors), drop = FALSE]

multi_data <- na.omit(multi_data)

cat("多因素 Cox 分析纳入样本数:", nrow(multi_data), "\n")

if (nrow(multi_data) == 0) {
  stop("多因素 Cox 分析去除 NA 后没有剩余样本，请检查 futime、fustat 或模型变量。")
}

formula <- as.formula(
  paste(
    "Surv(futime, fustat) ~",
    paste(predictors, collapse = " + ")
  )
)

multiCox <- coxph(
  formula,
  data = multi_data
)

multiCoxSummary <- summary(multiCox)

multiCoxResults <- data.frame(
  Variable = rownames(multiCoxSummary$coefficients),
  HR = multiCoxSummary$conf.int[, "exp(coef)"],
  HR.95L = multiCoxSummary$conf.int[, "lower .95"],
  HR.95H = multiCoxSummary$conf.int[, "upper .95"],
  pvalue = multiCoxSummary$coefficients[, "Pr(>|z|)"]
)

print(multiCoxResults)
multiCoxResults_path <- file.path(output_folder, "multiCoxResults.csv")
write.csv(multiCoxResults, file = multiCoxResults_path, row.names = FALSE)

# 15. 设置森林图通用参数

base_size <- 10

ci_pch <- 20

ci_col <- "#4575b4"

ci_lty <- 1

ci_lwd <- 2.3

ci_Theight <- 0.2

refline_lty <- "dashed"

footnote_cex <- 1

footnote_fontface <- "italic"

footnote_col <- "blue"

tm <- forest_theme(
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


# 16. 定义数值格式化函数

format_numbers <- function(x) {
  
  formatted <- format(
    round(x, 3),
    nsmall = 3
  )
  
  formatted[as.numeric(formatted) < 0.001] <- "<0.001"
  
  return(formatted)
}


# 17. 绘制多因素 Cox 森林图

final_results <- multiCoxResults

final_results_clean <- na.omit(final_results)

numeric_cols <- sapply(
  final_results_clean,
  is.numeric
)

final_results_clean[numeric_cols] <- lapply(
  final_results_clean[numeric_cols],
  format_numbers
)

dt <- final_results_clean

dt <- dt[, 1:5]

dt$` ` <- paste(
  rep(" ", 25),
  collapse = " "
)

dt$`HR(95%CI)` <- ifelse(
  is.na(dt$HR),
  "",
  sprintf(
    "%.3f(%.3f to %.3f)",
    as.numeric(dt$HR),
    as.numeric(dt$HR.95L),
    as.numeric(dt$HR.95H)
  )
)

dt[is.na(dt)] <- " "

dt$HR <- as.numeric(dt$HR)

dt$HR.95L <- as.numeric(dt$HR.95L)

dt$HR.95H <- as.numeric(dt$HR.95H)

dt <- dt[order(dt$pvalue), ]

p <- forest(
  dt[, c(1, 5, 6, 7)],
  est = dt$HR,
  lower = dt$HR.95L,
  upper = dt$HR.95H,
  sizes = sizes,
  ci_column = 3,
  ref_line = 1,
  xlim = xlimss,
  ticks_at = ticks_at,
  arrow_lab = arrow_lab,
  theme = tm
)

multi_forest_path <- file.path(
  output_folder,
  "forest森林图(多因素).pdf"
)

pdf(
  multi_forest_path,
  width = widthkkm,
  height = heightkkm
)

print(p)

dev.off()


# 18. 绘制单因素 Cox 森林图

final_results <- uniCoxResults

final_results_clean <- na.omit(final_results)

numeric_cols <- sapply(
  final_results_clean,
  is.numeric
)

final_results_clean[numeric_cols] <- lapply(
  final_results_clean[numeric_cols],
  format_numbers
)

dt <- final_results_clean

dt <- dt[, 1:5]

dt$` ` <- paste(
  rep(" ", 25),
  collapse = " "
)

dt$`HR(95%CI)` <- ifelse(
  is.na(dt$HR),
  "",
  sprintf(
    "%.3f(%.3f to %.3f)",
    as.numeric(dt$HR),
    as.numeric(dt$HR.95L),
    as.numeric(dt$HR.95H)
  )
)

dt[is.na(dt)] <- " "

dt$HR <- as.numeric(dt$HR)

dt$HR.95L <- as.numeric(dt$HR.95L)

dt$HR.95H <- as.numeric(dt$HR.95H)

dt <- dt[order(dt$pvalue), ]

p <- forest(
  dt[, c(1, 5, 6, 7)],
  est = dt$HR,
  lower = dt$HR.95L,
  upper = dt$HR.95H,
  sizes = sizes,
  ci_column = 3,
  ref_line = 1,
  xlim = xlimss,
  ticks_at = ticks_at,
  arrow_lab = arrow_lab,
  theme = tm
)

uni_forest_path <- file.path(
  output_folder,
  "forest森林图(单因素).pdf"
)

pdf(
  uni_forest_path,
  width = widthkkm,
  height = heightkkm
)

print(p)

dev.off()


# 19. 输出结果路径

cat("所有结果已保存到文件夹：", output_folder, "\n")

cat("单因素分析结果：", uniCoxResults_path, "\n")

cat("多因素分析结果：", multiCoxResults_path, "\n")

cat("多因素森林图：", multi_forest_path, "\n")

cat("单因素森林图：", uni_forest_path, "\n")