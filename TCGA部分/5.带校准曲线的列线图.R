# 1. 加载必要 R 包

library(data.table)
library(dplyr)
library(survival)
library(rms)
library(regplot)


# 2. 设置通用参数

input_file <- "合并后的标准肿瘤分析数据.csv"

model_vars <- c(
  "age",
  "gender",
  "stage",
  "T",
  "N",
  "FUT4"
)

plot_width <- 8
plot_height <- 6

sample_index <- 2

save_dir <- "列线图结果"

regression_type <- "coxph"
# 可选：
# regression_type <- "coxph"


# 3. 读取数据

data11221 <- fread(
  input_file,
  data.table = FALSE
)

rownames(data11221) <- data11221[, 1]

data11221 <- data11221[, -1, drop = FALSE]


# 4. 数据预处理

data_clean <- data11221

data_clean$futime <- data_clean$futime / 365

data_clean <- na.omit(data_clean)

data_clean[] <- lapply(
  data_clean,
  function(x) {
    if (is.character(x) | is.logical(x)) {
      as.factor(x)
    } else {
      x
    }
  }
)

dd <- suppressWarnings(datadist(data_clean))

options(datadist = "dd")


# 5. 构建 Cox 模型公式

model_formula <- as.formula(
  paste(
    "Surv(futime, fustat) ~",
    paste(model_vars, collapse = " + ")
  )
)


# 6. 创建结果保存文件夹

if (!dir.exists(save_dir)) {
  dir.create(save_dir)
}


# 7. 使用 rms::cph 进行 Cox 回归并绘制列线图

if (regression_type == "cox") {
  
  cox_model_clean <- cph(
    model_formula,
    data = data_clean,
    x = TRUE,
    y = TRUE,
    surv = TRUE
  )
  
  regplot_output_path <- file.path(
    save_dir,
    "regplot_ALL(cox).pdf"
  )
  
  regplot_sam_path <- file.path(
    save_dir,
    "regplot_指定(cox).pdf"
  )
  
  risk_score_file_path <- file.path(
    save_dir,
    "risk_score(cox).csv"
  )
  
  risk_score_file <- file.path(
    save_dir,
    "总表risk_score(cox).csv"
  )
  
  nom <- regplot(
    cox_model_clean,
    plots = c("density", "boxes"),
    clickable = FALSE,
    title = "",
    points = TRUE,
    droplines = TRUE,
    observation = NULL,
    rank = "sd",
    failtime = c(1, 3, 5),
    prfail = FALSE
  )
  
  dev.copy2pdf(
    file = regplot_output_path,
    width = plot_width,
    height = plot_height
  )
  
  nom1 <- regplot(
    cox_model_clean,
    plots = c("density", "boxes"),
    clickable = FALSE,
    title = "",
    points = TRUE,
    droplines = TRUE,
    observation = data_clean[sample_index, ],
    rank = "sd",
    failtime = c(1, 3, 5),
    prfail = FALSE
  )
  
  dev.copy2pdf(
    file = regplot_sam_path,
    width = plot_width,
    height = plot_height
  )
  
  nomoRisk <- predict(
    cox_model_clean,
    type = "lp"
  )
  
  relative_risk <- exp(nomoRisk)
  
  data_with_risk <- data_clean %>%
    dplyr::mutate(
      risk_score = relative_risk
    )
  
  dataRSK2 <- data_with_risk %>%
    dplyr::mutate(
      row_name = rownames(data_with_risk)
    ) %>%
    dplyr::select(
      risk_score
    )
  
  write.csv(
    dataRSK2,
    risk_score_file_path,
    row.names = TRUE
  )
  
  data_with_risk$futime <- data_with_risk$futime * 365
  
  write.csv(
    data_with_risk,
    risk_score_file,
    row.names = TRUE
  )
}


# 8. 使用 survival::coxph 进行 Cox 比例风险回归并绘制列线图

if (regression_type == "coxph") {
  
  coxph_model_clean <- coxph(
    model_formula,
    data = data_clean
  )
  
  regplot_output_path <- file.path(
    save_dir,
    "regplot_ALL(coxph).pdf"
  )
  
  regplot_sam_path <- file.path(
    save_dir,
    "regplot_指定(coxph).pdf"
  )
  
  risk_score_file_path <- file.path(
    save_dir,
    "risk_score(coxph).csv"
  )
  
  risk_score_file <- file.path(
    save_dir,
    "总表risk_score(coxph).csv"
  )
  
  nom <- regplot(
    coxph_model_clean,
    plots = c("density", "boxes"),
    clickable = FALSE,
    title = "",
    points = TRUE,
    droplines = TRUE,
    observation = NULL,
    rank = "sd",
    failtime = c(1, 3, 5),
    prfail = FALSE
  )
  
  dev.copy2pdf(
    file = regplot_output_path,
    width = plot_width,
    height = plot_height
  )
  
  nom1 <- regplot(
    coxph_model_clean,
    plots = c("density", "boxes"),
    clickable = FALSE,
    title = "",
    points = TRUE,
    droplines = TRUE,
    observation = data_clean[sample_index, ],
    rank = "sd",
    failtime = c(1, 3, 5),
    prfail = FALSE
  )
  
  dev.copy2pdf(
    file = regplot_sam_path,
    width = plot_width,
    height = plot_height
  )
  
  nomoRisk <- predict(
    coxph_model_clean,
    type = "risk"
  )
  
  data_with_risk <- data_clean %>%
    dplyr::mutate(
      risk_score = nomoRisk
    )
  
  dataRSK1 <- data_with_risk %>%
    dplyr::mutate(
      row_name = rownames(data_with_risk)
    ) %>%
    dplyr::select(
      risk_score
    )
  
  write.csv(
    dataRSK1,
    risk_score_file_path,
    row.names = TRUE
  )
  data_with_risk$futime <- data_with_risk$futime * 365
  write.csv(
    data_with_risk,
    risk_score_file,
    row.names = TRUE
  )
}

# 9. 绘制列线图校准曲线

if (regression_type == "cox") {
  
  cal_time_points <- c(1, 3, 5)
  
  cal_B <- 1000
  
  cal_m <- floor(nrow(data_clean) / 5)
  
  if (cal_m < 10) {
    cal_m <- 10
  }
  
  if (cal_m >= nrow(data_clean)) {
    cal_m <- max(2, floor(nrow(data_clean) / 2))
  }
  
  for (cal_time in cal_time_points) {
    
    cal_result <- calibrate(
      cox_model_clean,
      cmethod = "KM",
      method = "boot",
      u = cal_time,
      m = cal_m,
      B = cal_B
    )
    
    cal_output_path <- file.path(
      save_dir,
      paste0("Calibration_", cal_time, "year(cox).pdf")
    )
    
    pdf(
      cal_output_path,
      width = plot_width,
      height = plot_height
    )
    
    plot(
      cal_result,
      xlab = paste0("Predicted ", cal_time, "-year survival probability"),
      ylab = paste0("Observed ", cal_time, "-year survival probability"),
      subtitles = FALSE
    )
    
    abline(
      0,
      1,
      lty = 2,
      col = "gray"
    )
    
    dev.off()
    
    message(
      paste0(
        cal_time,
        " 年校准曲线已保存: ",
        cal_output_path
      )
    )
  }
  
} else {
  
  message(
    "当前 regression_type 不是 'cox'，未运行 rms::cph() 模型，因此没有生成 calibrate() 校准曲线。"
  )
}