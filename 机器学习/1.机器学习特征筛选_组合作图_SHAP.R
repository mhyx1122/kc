# 1. 加载必要R包
suppressPackageStartupMessages({
  library(caret)
  library(doParallel)
  library(ggplot2)
  library(glmnet)
  library(randomForest)
  library(pROC)
  library(DALEX)
  library(rmda)
  library(ggsci)
  library(lattice)
  library(rpart)
  library(C50)
  library(gbm)
  library(xgboost)
  library(nnet)
  library(kernlab)
})

# 2. 设置输入参数
file_path <- "示例数据1.csv"
case_num <- 541
control_num <- 59
do_log2 <- TRUE
train_prop <- 0.7
random_seed <- 123

number_gene <- 10
cv_folds <- 10
cv_repeats <- 1
cores <- 6

data_option <- "train_data"

positive_class <- "Case"
negative_class <- "Control"

# 2.1 设置需要运行的模型
# 可选模型名称：
# lasso, rf, svm, gbm, xgb, dt, knn, nnet, glm, c50
# 写进去的模型会运行，不写进去的模型会自动跳过

run_models <- c(
  "lasso", "rf", "svm", "gbm", "xgb",
  "dt", "knn", "nnet", "glm", "c50"
)

# 示例1：只运行 RF、SVM、XGB
# run_models <- c("rf", "svm", "xgb")

# 示例2：只运行 LASSO 和 GLM
# run_models <- c("lasso", "glm")

# 示例3：只运行 RF
# run_models <- c("rf")

all_model_names <- c(
  "lasso", "rf", "svm", "gbm", "xgb",
  "dt", "knn", "nnet", "glm", "c50"
)

run_models <- unique(tolower(run_models))

invalid_models <- setdiff(run_models, all_model_names)

if (length(invalid_models) > 0) {
  stop(
    paste0(
      "run_models 中存在不支持的模型名称：",
      paste(invalid_models, collapse = ", ")
    )
  )
}

if (length(run_models) == 0) {
  stop("run_models 为空，没有需要运行的模型。")
}

should_run_model <- function(model_name) {
  tolower(model_name) %in% run_models
}

# 3. 定义预测函数
p_fun <- function(object, newdata) {
  pred <- predict(object, newdata = newdata, type = "prob")
  
  if (!positive_class %in% colnames(pred)) {
    stop("预测概率中没有 positive_class 对应的列，请检查分组因子水平。")
  }
  
  pred[, positive_class]
}

# 4. 自动判断并执行log2转换
auto_log2_transform <- function(expr_data, do_log2 = TRUE) {
  if (!do_log2) {
    return(expr_data)
  }
  
  qx <- as.numeric(
    quantile(
      as.matrix(expr_data),
      c(0.00, 0.25, 0.5, 0.75, 0.99, 1.0),
      na.rm = TRUE
    )
  )
  
  need_log2 <- (qx[5] > 100) ||
    (qx[6] - qx[1] > 50 && qx[2] > 0)
  
  if (need_log2) {
    message("log2 transform finished")
    return(log2(expr_data + 1))
  }
  
  message("log2 transform not needed")
  expr_data
}

# 5. 保存特征重要性条形图
save_importance_barplot <- function(importance_data,
                                    plot_title,
                                    out_file) {
  p <- ggplot(
    importance_data,
    aes(x = reorder(var, impor), y = impor)
  ) +
    geom_col(fill = "#69b3a2", color = "black") +
    labs(
      title = plot_title,
      x = "",
      y = "Importance"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 12, hjust = 0.5),
      axis.text.x = element_text(size = 8, angle = 0, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 10, color = "black", hjust = 1),
      axis.title = element_text(size = 12)
    ) +
    coord_flip() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
  
  ggsave(out_file, plot = p, width = 8, height = 6, units = "in")
  
  p
}

# 6. 模型评估函数
evaluate_model <- function(model,
                           test_data,
                           model_label,
                           out_dir,
                           roc_file,
                           calibration_file,
                           dca_file,
                           metrics_file,
                           set_colors,
                           positive_class = "Case") {
  
  x_test <- test_data[, setdiff(colnames(test_data), "Type"), drop = FALSE]
  y_test <- ifelse(test_data$Type == positive_class, 1, 0)
  
  if (length(unique(test_data$Type)) < 2) {
    stop(paste0(model_label, "测试集中只有一个类别，无法计算ROC、校准曲线、DCA和混淆矩阵。"))
  }
  
  pred_prob <- predict(model, newdata = x_test, type = "prob")
  
  if (!positive_class %in% colnames(pred_prob)) {
    stop(paste0(model_label, "模型预测概率中没有", positive_class, "列。"))
  }
  
  prob_case <- pred_prob[, positive_class]
  
  roc_obj <- pROC::roc(y_test, prob_case, quiet = TRUE)
  
  pdf(file = file.path(out_dir, roc_file), width = 5, height = 5)
  plot(roc_obj, print.auc = FALSE, legacy.axes = TRUE, main = "", col = "red")
  legend(
    x = 0.6,
    y = 0.1,
    paste0(model_label, " (AUC = ", sprintf("%.03f", roc_obj$auc), ")"),
    col = "red",
    lwd = 2,
    bty = "n"
  )
  dev.off()
  
  brier_score <- mean((prob_case - y_test)^2)
  
  test_cal <- test_data
  test_cal$prob_case <- prob_case
  
  cal <- caret::calibration(
    Type ~ prob_case,
    data = test_cal,
    class = positive_class
  )
  
  pdf(file = file.path(out_dir, calibration_file), width = 8, height = 6)
  print(lattice::xyplot(
    Percent ~ midpoint,
    data = cal$data,
    auto.key = list(columns = 1),
    main = paste("Calibration Curve (Brier Score =", round(brier_score, 3), ")"),
    panel = function(x, y, ...) {
      lattice::panel.abline(a = 0, b = 1, col = "darkgray", lwd = 2, lty = 2)
      lattice::panel.xyplot(x, y, col = "red", lwd = 3, type = "l", ...)
      lattice::panel.points(x, y, col = "red", pch = 19, cex = 1.2)
    },
    scales = list(
      x = list(draw = TRUE, tck = c(0.8, 0)),
      y = list(draw = TRUE, tck = c(0.8, 0))
    ),
    par.settings = list(
      axis.line = list(col = "black", lty = 1, lwd = 1)
    ),
    xlab = "Mean Predicted Probability",
    ylab = "Observed Event Probability"
  ))
  dev.off()
  
  decision_data <- data.frame(
    outcome = y_test,
    model_prob = prob_case
  )
  
  dca_obj <- rmda::decision_curve(
    outcome ~ model_prob,
    data = decision_data,
    family = binomial(link = "logit"),
    thresholds = seq(0, 1, by = 0.01),
    fitted.risk = TRUE
  )
  
  pdf(file = file.path(out_dir, dca_file), width = 8, height = 6)
  rmda::plot_decision_curve(
    dca_obj,
    curve.names = model_label,
    xlab = "Threshold Probability",
    ylab = "Net Benefit",
    col = set_colors,
    lwd = 2,
    legend.position = "left"
  )
  dev.off()
  
  pred_label <- predict(model, newdata = x_test)
  pred_label <- factor(pred_label, levels = levels(test_data$Type))
  true_label <- factor(test_data$Type, levels = levels(test_data$Type))
  
  cm <- caret::confusionMatrix(
    pred_label,
    true_label,
    positive = positive_class
  )
  
  metrics <- data.frame(
    Metric = c("AUROC", "Sensitivity", "Specificity", "Accuracy", "Kappa", "F1"),
    Value = round(c(
      as.numeric(pROC::auc(roc_obj)),
      cm$byClass["Sensitivity"],
      cm$byClass["Specificity"],
      cm$overall["Accuracy"],
      cm$overall["Kappa"],
      cm$byClass["F1"]
    ), 3)
  )
  
  rownames(metrics) <- NULL
  
  metrics_wide <- as.data.frame(t(metrics))
  colnames(metrics_wide) <- metrics_wide[1, ]
  metrics_wide <- metrics_wide[-1, , drop = FALSE]
  
  write.csv(
    metrics_wide,
    file.path(out_dir, metrics_file),
    row.names = TRUE
  )
  
  explainer <- DALEX::explain(
    model,
    label = model_label,
    data = x_test,
    y = y_test,
    predict_function = p_fun,
    verbose = FALSE
  )
  
  model_perf <- DALEX::model_performance(explainer)
  
  list(
    roc = roc_obj,
    dca = dca_obj,
    explainer = explainer,
    model_performance = model_perf,
    metrics = metrics_wide
  )
}

# 7. 读取表达矩阵
raw_expr <- read.csv(
  file_path,
  header = TRUE,
  check.names = FALSE,
  row.names = 1
)

raw_expr <- as.data.frame(t(raw_expr))
raw_expr[] <- lapply(raw_expr, as.numeric)

if (anyNA(raw_expr)) {
  stop("表达矩阵中存在NA。请检查是否有非数字字符、空值或异常符号。")
}

# 8. 设置分组
sample_group <- factor(
  c(
    rep(positive_class, case_num),
    rep(negative_class, control_num)
  ),
  levels = c(positive_class, negative_class)
)

if (length(sample_group) != nrow(raw_expr)) {
  stop("分组数量和样本数量不一致：case_num + control_num 必须等于样本数。")
}

# 9. log2转换并添加分组
data_all <- auto_log2_transform(raw_expr, do_log2 = do_log2)
data_all$Type <- sample_group

# 10. 划分训练集和测试集
set.seed(random_seed)

train_index <- caret::createDataPartition(
  y = data_all$Type,
  p = train_prop,
  list = FALSE
)

train_index <- as.vector(train_index)

train_data <- data_all[train_index, , drop = FALSE]
test_data <- data_all[-train_index, , drop = FALSE]

# 11. 选择特征筛选使用的数据
if (!data_option %in% c("data_type", "train_data")) {
  stop("data_option 只能是 'data_type' 或 'train_data'。")
}

if (data_option == "data_type") {
  feature_data <- data_all
} else {
  feature_data <- train_data
}

# 12. 检查交叉验证折数
min_class_n <- min(table(feature_data$Type))

if (cv_folds > min_class_n) {
  stop("cv_folds 大于训练数据中最小类别样本数，请调小 cv_folds。")
}

# 13. 注册并行计算
if (cores > parallel::detectCores()) {
  stop("错误：核心数超出系统限制，请重新设置。")
}

doParallel::registerDoParallel(cores = cores)

# 14. 设置决策曲线颜色
set_colors <- c(
  ggsci::pal_npg()(9),
  ggsci::pal_jco()(9),
  ggsci::pal_jama()(7),
  ggsci::pal_nejm()(8)
)

# 15. 初始化结果列表
results <- list()

# 16. LASSO模型
if (should_run_model("lasso")) {
  
  message("LASSO运行中···")
  
  results$lasso_gene <- tryCatch({
    
    out_dir <- "lasso"
    if (!dir.exists(out_dir)) {
      dir.create(out_dir)
    }
    
    group1 <- feature_data$Type
    data0 <- feature_data[, setdiff(colnames(feature_data), "Type"), drop = FALSE]
    
    x <- as.matrix(data0)
    y <- group1
    
    fit <- glmnet::glmnet(x, y, family = "binomial", alpha = 1)
    
    pdf(file.path(out_dir, "lambda.pdf"), width = 6, height = 5.5)
    plot(fit, xvar = "lambda", label = FALSE)
    dev.off()
    
    cvfit <- glmnet::cv.glmnet(
      x,
      y,
      family = "binomial",
      alpha = 1,
      type.measure = "deviance",
      nfolds = cv_folds
    )
    
    pdf(file.path(out_dir, "cvfit.pdf"), width = 6, height = 5.5)
    plot(cvfit)
    dev.off()
    
    coef_mat <- as.matrix(coef(fit, s = cvfit$lambda.min))
    gene_index <- which(coef_mat != 0)
    
    lasso_genes <- rownames(coef_mat)[gene_index]
    lasso_genes <- setdiff(lasso_genes, "(Intercept)")
    
    if (length(lasso_genes) == 0) {
      stop("LASSO没有筛选出特征基因。")
    }
    
    write.table(
      lasso_genes,
      file = file.path(out_dir, "LASSO特征基因.csv"),
      sep = ",",
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE
    )
    
    gene_expr <- t(raw_expr)
    gene_expr <- auto_log2_transform(gene_expr, do_log2 = do_log2)
    gene_expr <- gene_expr[lasso_genes, , drop = FALSE]
    gene_expr <- as.data.frame(gene_expr)
    
    write.table(
      gene_expr,
      file = file.path(out_dir, "LASSO特征基因表达量.csv"),
      sep = ",",
      quote = FALSE,
      row.names = TRUE,
      col.names = TRUE
    )
    
    lasso_data <- as.data.frame(t(gene_expr))
    lasso_data$Type <- sample_group
    
    train <- lasso_data[train_index, , drop = FALSE]
    test <- lasso_data[-train_index, , drop = FALSE]
    
    control <- caret::trainControl(
      method = "repeatedcv",
      number = cv_folds,
      repeats = cv_repeats,
      savePredictions = TRUE
    )
    
    mod_lasso <- caret::train(
      Type ~ .,
      data = train,
      method = "glmnet",
      trControl = control
    )
    
    saveRDS(mod_lasso, file.path(out_dir, "mod_lasso.rds"))
    
    eval_result <- evaluate_model(
      model = mod_lasso,
      test_data = test,
      model_label = "LASSO",
      out_dir = out_dir,
      roc_file = "ROC_LASSO.pdf",
      calibration_file = "Calibration_Curve.pdf",
      dca_file = "Decision_Curve_LASSO.pdf",
      metrics_file = "LASSO_Performance_Results.csv",
      set_colors = set_colors,
      positive_class = positive_class
    )
    
    list(
      exp = gene_expr,
      model = mod_lasso,
      roc = eval_result$roc,
      dca = eval_result$dca,
      explainer = eval_result$explainer,
      model_performance = eval_result$model_performance,
      metrics = eval_result$metrics
    )
    
  }, error = function(e) {
    message("Error in LASSO gene selection: ", e$message)
    NULL
  })
  
} else {
  message("LASSO未选择运行。")
  results$lasso_gene <- NULL
}

# 17. RF模型
if (should_run_model("rf")) {
  
  message("RF运行中···")
  
  results$rf_gene <- tryCatch({
    
    out_dir <- "rf"
    if (!dir.exists(out_dir)) {
      dir.create(out_dir)
    }
    
    group1 <- feature_data$Type
    data0 <- feature_data[, setdiff(colnames(feature_data), "Type"), drop = FALSE]
    
    rf_data <- data0
    rf_data$Type <- group1
    
    rf <- randomForest::randomForest(Type ~ ., data = rf_data, ntree = 1000)
    
    pdf(file.path(out_dir, "随机森林误差图.pdf"), width = 6, height = 6)
    plot(rf, main = "Random forest", lwd = 2)
    dev.off()
    
    option_trees <- which.min(rf$err.rate[, 1])
    rf2 <- randomForest::randomForest(Type ~ ., data = rf_data, ntree = option_trees)
    
    importance_mat <- randomForest::importance(rf2)
    
    if (!"MeanDecreaseGini" %in% colnames(importance_mat)) {
      stop("随机森林重要性结果中没有MeanDecreaseGini列。")
    }
    
    varimpdf_top <- data.frame(
      var = rownames(importance_mat),
      impor = importance_mat[, "MeanDecreaseGini"]
    )
    
    varimpdf_top <- varimpdf_top[order(-varimpdf_top$impor), , drop = FALSE]
    varimpdf_top <- head(varimpdf_top, number_gene)
    
    write.csv(
      varimpdf_top,
      file.path(out_dir, "基因重要度得分表.csv"),
      row.names = FALSE,
      quote = FALSE
    )
    
    save_importance_barplot(
      importance_data = varimpdf_top,
      plot_title = "Gene Importance (Random Forest)",
      out_file = file.path(out_dir, "rf重要度排序.pdf")
    )
    
    pdf(file.path(out_dir, "基因重要度排序.pdf"), width = 6.2, height = 5.8)
    randomForest::varImpPlot(rf2, main = "", n.var = 14)
    dev.off()
    
    rf_genes <- varimpdf_top$var
    
    write.table(
      rf_genes,
      file = file.path(out_dir, "随机森林特征基因.csv"),
      sep = ",",
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE
    )
    
    gene_expr <- t(raw_expr)
    gene_expr <- auto_log2_transform(gene_expr, do_log2 = do_log2)
    gene_expr <- gene_expr[rf_genes, , drop = FALSE]
    gene_expr <- as.data.frame(gene_expr)
    
    write.table(
      gene_expr,
      file = file.path(out_dir, "随机森林特征基因表达量.csv"),
      sep = ",",
      quote = FALSE,
      row.names = TRUE,
      col.names = TRUE
    )
    
    rf_model_data <- as.data.frame(t(gene_expr))
    rf_model_data$Type <- sample_group
    
    train <- rf_model_data[train_index, , drop = FALSE]
    test <- rf_model_data[-train_index, , drop = FALSE]
    
    control <- caret::trainControl(
      method = "repeatedcv",
      number = cv_folds,
      repeats = cv_repeats,
      savePredictions = TRUE
    )
    
    mod_rf <- caret::train(
      Type ~ .,
      data = train,
      method = "rf",
      trControl = control
    )
    
    saveRDS(mod_rf, file.path(out_dir, "mod_rf.rds"))
    
    eval_result <- evaluate_model(
      model = mod_rf,
      test_data = test,
      model_label = "RF",
      out_dir = out_dir,
      roc_file = "ROC_RF.pdf",
      calibration_file = "Calibration_Curve.pdf",
      dca_file = "Decision_Curve_RF.pdf",
      metrics_file = "RF_Performance_Results.csv",
      set_colors = set_colors,
      positive_class = positive_class
    )
    
    list(
      exp = gene_expr,
      model = mod_rf,
      roc = eval_result$roc,
      dca = eval_result$dca,
      explainer = eval_result$explainer,
      model_performance = eval_result$model_performance,
      metrics = eval_result$metrics,
      gene_importance = varimpdf_top
    )
    
  }, error = function(e) {
    message("Error in RF gene selection: ", e$message)
    NULL
  })
  
} else {
  message("RF未选择运行。")
  results$rf_gene <- NULL
}

# 18. SVM-RFE模型
if (should_run_model("svm")) {
  
  message("SVM运行中···")
  
  results$svm_gene <- tryCatch({
    
    out_dir <- "svm"
    if (!dir.exists(out_dir)) {
      dir.create(out_dir)
    }
    
    group1 <- feature_data$Type
    data0 <- feature_data[, setdiff(colnames(feature_data), "Type"), drop = FALSE]
    
    ctrl <- caret::rfeControl(
      functions = caretFuncs,
      method = "cv",
      number = cv_folds
    )
    
    sizes_use <- seq(1, min(30, ncol(data0)), by = 3)
    
    profile <- caret::rfe(
      x = data0,
      y = group1,
      sizes = sizes_use,
      rfeControl = ctrl,
      method = "svmRadial"
    )
    
    pdf(file.path(out_dir, "SVM-RFE递归消除验证.pdf"), width = 6, height = 5.5)
    par(las = 1)
    x <- profile$results$Variables
    y <- profile$results$Accuracy
    plot(x, y, xlab = "Variables", ylab = "", col = "darkgreen")
    lines(x, y, col = "darkgreen")
    wmax <- which.max(y)
    points(x[wmax], y[wmax], col = "blue", pch = 16)
    text(x[wmax], y[wmax], paste0("N=", x[wmax]), pos = 4, col = 2)
    mtext(
      "Accuracy (Cross-Validation)",
      side = 3,
      line = 1.5,
      at = mean(range(x)) - 1,
      col = "black",
      cex = 1.2
    )
    dev.off()
    
    svm_genes <- profile$optVariables
    
    write.table(
      svm_genes,
      file = file.path(out_dir, "SVM-RFE.gene.csv"),
      sep = ",",
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE
    )
    
    gene_expr <- t(raw_expr)
    gene_expr <- auto_log2_transform(gene_expr, do_log2 = do_log2)
    gene_expr <- gene_expr[svm_genes, , drop = FALSE]
    gene_expr <- as.data.frame(gene_expr)
    
    write.table(
      gene_expr,
      file = file.path(out_dir, "SVM.geneExp.csv"),
      sep = ",",
      quote = FALSE,
      row.names = TRUE,
      col.names = TRUE
    )
    
    svm_model_data <- as.data.frame(t(gene_expr))
    svm_model_data$Type <- sample_group
    
    train <- svm_model_data[train_index, , drop = FALSE]
    test <- svm_model_data[-train_index, , drop = FALSE]
    
    control <- caret::trainControl(
      method = "repeatedcv",
      number = cv_folds,
      repeats = cv_repeats,
      savePredictions = TRUE
    )
    
    mod_svm <- caret::train(
      Type ~ .,
      data = train,
      method = "svmRadial",
      prob.model = TRUE,
      trControl = control
    )
    
    saveRDS(mod_svm, file.path(out_dir, "mod_svm.rds"))
    
    eval_result <- evaluate_model(
      model = mod_svm,
      test_data = test,
      model_label = "SVM",
      out_dir = out_dir,
      roc_file = "ROC_svm.pdf",
      calibration_file = "Calibration_Curve.pdf",
      dca_file = "Decision_Curve_SVM.pdf",
      metrics_file = "SVM_Performance_Results.csv",
      set_colors = set_colors,
      positive_class = positive_class
    )
    
    list(
      exp = gene_expr,
      model = mod_svm,
      roc = eval_result$roc,
      dca = eval_result$dca,
      explainer = eval_result$explainer,
      model_performance = eval_result$model_performance,
      metrics = eval_result$metrics
    )
    
  }, error = function(e) {
    message("Error in SVM gene selection: ", e$message)
    NULL
  })
  
} else {
  message("SVM未选择运行。")
  results$svm_gene <- NULL
}

# 19. GBM模型
if (should_run_model("gbm")) {
  
  message("GBM运行中···")
  
  results$gbm_gene <- tryCatch({
    
    out_dir <- "gbm"
    if (!dir.exists(out_dir)) {
      dir.create(out_dir)
    }
    
    group1 <- feature_data$Type
    data0 <- feature_data[, setdiff(colnames(feature_data), "Type"), drop = FALSE]
    
    fit_control <- caret::trainControl(
      method = "repeatedcv",
      number = cv_folds,
      repeats = cv_repeats
    )
    
    gbm_feature_model <- caret::train(
      x = data0,
      y = group1,
      method = "gbm",
      trControl = fit_control,
      verbose = FALSE
    )
    
    importance <- caret::varImp(gbm_feature_model, scale = TRUE)
    importance_df <- as.data.frame(importance$importance)
    
    varimpdf_top <- data.frame(
      var = rownames(importance_df),
      impor = importance_df[, 1]
    )
    
    varimpdf_top <- varimpdf_top[order(-varimpdf_top$impor), , drop = FALSE]
    varimpdf_top <- head(varimpdf_top, number_gene)
    
    write.csv(
      varimpdf_top,
      file.path(out_dir, "基因重要度数值表.csv"),
      row.names = FALSE
    )
    
    save_importance_barplot(
      importance_data = varimpdf_top,
      plot_title = "Gene Importance (Gradient Boosting Machine)",
      out_file = file.path(out_dir, "GBM重要度排序.pdf")
    )
    
    gbm_genes <- varimpdf_top$var
    
    write.table(
      gbm_genes,
      file = file.path(out_dir, "GBM特征基因.gene.csv"),
      sep = ",",
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE
    )
    
    gene_expr <- t(raw_expr)
    gene_expr <- auto_log2_transform(gene_expr, do_log2 = do_log2)
    gene_expr <- gene_expr[gbm_genes, , drop = FALSE]
    gene_expr <- as.data.frame(gene_expr)
    
    write.table(
      gene_expr,
      file = file.path(out_dir, "GBM特征基因表达量.csv"),
      sep = ",",
      quote = FALSE,
      row.names = TRUE,
      col.names = TRUE
    )
    
    gbm_model_data <- as.data.frame(t(gene_expr))
    gbm_model_data$Type <- sample_group
    
    train <- gbm_model_data[train_index, , drop = FALSE]
    test <- gbm_model_data[-train_index, , drop = FALSE]
    
    final_control <- caret::trainControl(
      method = "repeatedcv",
      number = cv_folds,
      repeats = cv_repeats,
      savePredictions = TRUE
    )
    
    mod_gbm <- caret::train(
      Type ~ .,
      data = train,
      method = "gbm",
      trControl = final_control,
      verbose = FALSE
    )
    
    saveRDS(mod_gbm, file.path(out_dir, "mod_gbm.rds"))
    
    eval_result <- evaluate_model(
      model = mod_gbm,
      test_data = test,
      model_label = "GBM",
      out_dir = out_dir,
      roc_file = "ROC_GBM.pdf",
      calibration_file = "Calibration_Curve.pdf",
      dca_file = "Decision_Curve_GBM.pdf",
      metrics_file = "GBM_Performance_Results.csv",
      set_colors = set_colors,
      positive_class = positive_class
    )
    
    list(
      exp = gene_expr,
      model = mod_gbm,
      roc = eval_result$roc,
      dca = eval_result$dca,
      explainer = eval_result$explainer,
      model_performance = eval_result$model_performance,
      metrics = eval_result$metrics,
      gene_importance = varimpdf_top
    )
    
  }, error = function(e) {
    message("Error in GBM gene selection: ", e$message)
    NULL
  })
  
} else {
  message("GBM未选择运行。")
  results$gbm_gene <- NULL
}

# 20. XGBoost模型
if (should_run_model("xgb")) {
  
  message("XGBoost运行中···")
  
  results$xgb_gene <- tryCatch({
    
    out_dir <- "xgboost"
    if (!dir.exists(out_dir)) {
      dir.create(out_dir)
    }
    
    group1 <- feature_data$Type
    data0 <- feature_data[, setdiff(colnames(feature_data), "Type"), drop = FALSE]
    
    fit_control <- caret::trainControl(
      method = "repeatedcv",
      number = cv_folds,
      repeats = cv_repeats
    )
    
    xgb_feature_model <- caret::train(
      x = data0,
      y = group1,
      method = "xgbTree",
      trControl = fit_control
    )
    
    importance <- caret::varImp(xgb_feature_model, scale = TRUE)
    importance_df <- as.data.frame(importance$importance)
    
    varimpdf_top <- data.frame(
      var = rownames(importance_df),
      impor = importance_df[, 1]
    )
    
    varimpdf_top <- varimpdf_top[order(-varimpdf_top$impor), , drop = FALSE]
    varimpdf_top <- head(varimpdf_top, number_gene)
    
    write.csv(
      varimpdf_top,
      file.path(out_dir, "基因重要度得分表.csv"),
      row.names = FALSE,
      quote = FALSE
    )
    
    save_importance_barplot(
      importance_data = varimpdf_top,
      plot_title = "Gene Importance (eXtreme Gradient Boosting)",
      out_file = file.path(out_dir, "xgboost重要度排序.pdf")
    )
    
    xgb_genes <- varimpdf_top$var
    
    write.table(
      xgb_genes,
      file = file.path(out_dir, "xgboost特征基因.gene.csv"),
      sep = ",",
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE
    )
    
    gene_expr <- t(raw_expr)
    gene_expr <- auto_log2_transform(gene_expr, do_log2 = do_log2)
    gene_expr <- gene_expr[xgb_genes, , drop = FALSE]
    gene_expr <- as.data.frame(gene_expr)
    
    write.table(
      gene_expr,
      file = file.path(out_dir, "xgboost特征基因表达量.csv"),
      sep = ",",
      quote = FALSE,
      row.names = TRUE,
      col.names = TRUE
    )
    
    xgb_model_data <- as.data.frame(t(gene_expr))
    xgb_model_data$Type <- sample_group
    
    train <- xgb_model_data[train_index, , drop = FALSE]
    test <- xgb_model_data[-train_index, , drop = FALSE]
    
    final_control <- caret::trainControl(
      method = "repeatedcv",
      number = cv_folds,
      repeats = cv_repeats,
      savePredictions = TRUE
    )
    
    mod_xgb <- caret::train(
      Type ~ .,
      data = train,
      method = "xgbTree",
      trControl = final_control
    )
    
    saveRDS(mod_xgb, file.path(out_dir, "mod_xgb.rds"))
    
    eval_result <- evaluate_model(
      model = mod_xgb,
      test_data = test,
      model_label = "XGB",
      out_dir = out_dir,
      roc_file = "ROC_xgboost.pdf",
      calibration_file = "Calibration_Curve.pdf",
      dca_file = "Decision_Curve_XGB.pdf",
      metrics_file = "XGB_Performance_Results.csv",
      set_colors = set_colors,
      positive_class = positive_class
    )
    
    list(
      exp = gene_expr,
      model = mod_xgb,
      roc = eval_result$roc,
      dca = eval_result$dca,
      explainer = eval_result$explainer,
      model_performance = eval_result$model_performance,
      metrics = eval_result$metrics,
      gene_importance = varimpdf_top
    )
    
  }, error = function(e) {
    message("Error in XGBoost gene selection: ", e$message)
    NULL
  })
  
} else {
  message("XGBoost未选择运行。")
  results$xgb_gene <- NULL
}

# 21. 决策树模型
if (should_run_model("dt")) {
  
  message("DecisionTree运行中···")
  
  results$dt_gene <- tryCatch({
    
    out_dir <- "DecisionTree"
    if (!dir.exists(out_dir)) {
      dir.create(out_dir)
    }
    
    group1 <- feature_data$Type
    data0 <- feature_data[, setdiff(colnames(feature_data), "Type"), drop = FALSE]
    
    dt_data <- as.data.frame(data0)
    dt_data$Type <- group1
    
    mod_dt_feature <- rpart::rpart(
      Type ~ .,
      data = dt_data,
      method = "class"
    )
    
    importance <- caret::varImp(mod_dt_feature, scale = TRUE)
    
    varimpdf_top <- data.frame(
      var = rownames(importance),
      impor = importance$Overall
    )
    
    varimpdf_top <- varimpdf_top[order(-varimpdf_top$impor), , drop = FALSE]
    varimpdf_top <- head(varimpdf_top, number_gene)
    
    write.csv(
      varimpdf_top,
      file.path(out_dir, "基因重要度数值表.csv"),
      row.names = FALSE
    )
    
    save_importance_barplot(
      importance_data = varimpdf_top,
      plot_title = "Gene Importance (Decision Tree)",
      out_file = file.path(out_dir, "DecisionTree重要度排序.pdf")
    )
    
    dt_genes <- varimpdf_top$var
    
    write.table(
      dt_genes,
      file = file.path(out_dir, "DecisionTree特征基因.gene.csv"),
      sep = ",",
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE
    )
    
    gene_expr <- t(raw_expr)
    gene_expr <- auto_log2_transform(gene_expr, do_log2 = do_log2)
    gene_expr <- gene_expr[dt_genes, , drop = FALSE]
    gene_expr <- as.data.frame(gene_expr)
    
    write.table(
      gene_expr,
      file = file.path(out_dir, "DecisionTree特征基因表达量.csv"),
      sep = ",",
      quote = FALSE,
      row.names = TRUE,
      col.names = TRUE
    )
    
    dt_model_data <- as.data.frame(t(gene_expr))
    dt_model_data$Type <- sample_group
    
    train <- dt_model_data[train_index, , drop = FALSE]
    test <- dt_model_data[-train_index, , drop = FALSE]
    
    final_control <- caret::trainControl(
      method = "repeatedcv",
      number = cv_folds,
      repeats = cv_repeats,
      savePredictions = TRUE
    )
    
    mod_dt <- caret::train(
      Type ~ .,
      data = train,
      method = "rpart",
      trControl = final_control
    )
    
    saveRDS(mod_dt, file.path(out_dir, "mod_dt.rds"))
    
    eval_result <- evaluate_model(
      model = mod_dt,
      test_data = test,
      model_label = "DT",
      out_dir = out_dir,
      roc_file = "ROC_DT.pdf",
      calibration_file = "Calibration_Curve.pdf",
      dca_file = "Decision_Curve_DT.pdf",
      metrics_file = "DT_Performance_Results.csv",
      set_colors = set_colors,
      positive_class = positive_class
    )
    
    list(
      exp = gene_expr,
      model = mod_dt,
      roc = eval_result$roc,
      dca = eval_result$dca,
      explainer = eval_result$explainer,
      model_performance = eval_result$model_performance,
      metrics = eval_result$metrics,
      gene_importance = varimpdf_top
    )
    
  }, error = function(e) {
    message("Error in Decision Tree gene selection: ", e$message)
    NULL
  })
  
} else {
  message("DecisionTree未选择运行。")
  results$dt_gene <- NULL
}

# 22. KNN模型
if (should_run_model("knn")) {
  
  message("KNN运行中···")
  
  results$knn_gene <- tryCatch({
    
    out_dir <- "knn"
    if (!dir.exists(out_dir)) {
      dir.create(out_dir)
    }
    
    group1 <- feature_data$Type
    data0 <- feature_data[, setdiff(colnames(feature_data), "Type"), drop = FALSE]
    
    fit_control <- caret::trainControl(
      method = "repeatedcv",
      number = cv_folds,
      repeats = cv_repeats
    )
    
    knn_feature_model <- caret::train(
      x = data0,
      y = group1,
      method = "knn",
      trControl = fit_control
    )
    
    importance <- caret::varImp(knn_feature_model, scale = TRUE)
    importance_df <- as.data.frame(importance$importance)
    
    varimpdf_top <- data.frame(
      var = rownames(importance_df),
      impor = importance_df[, 1]
    )
    
    varimpdf_top <- varimpdf_top[order(-varimpdf_top$impor), , drop = FALSE]
    varimpdf_top <- head(varimpdf_top, number_gene)
    
    write.csv(
      varimpdf_top,
      file.path(out_dir, "基因重要度数值表.csv"),
      row.names = FALSE
    )
    
    save_importance_barplot(
      importance_data = varimpdf_top,
      plot_title = "Gene Importance (K-Nearest Neighbors)",
      out_file = file.path(out_dir, "knn重要度排序.pdf")
    )
    
    knn_genes <- varimpdf_top$var
    
    write.table(
      knn_genes,
      file = file.path(out_dir, "knn特征基因.gene.csv"),
      sep = ",",
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE
    )
    
    gene_expr <- t(raw_expr)
    gene_expr <- auto_log2_transform(gene_expr, do_log2 = do_log2)
    gene_expr <- gene_expr[knn_genes, , drop = FALSE]
    gene_expr <- as.data.frame(gene_expr)
    
    write.table(
      gene_expr,
      file = file.path(out_dir, "knn特征基因表达量.csv"),
      sep = ",",
      quote = FALSE,
      row.names = TRUE,
      col.names = TRUE
    )
    
    knn_model_data <- as.data.frame(t(gene_expr))
    knn_model_data$Type <- sample_group
    
    train <- knn_model_data[train_index, , drop = FALSE]
    test <- knn_model_data[-train_index, , drop = FALSE]
    
    final_control <- caret::trainControl(
      method = "repeatedcv",
      number = cv_folds,
      repeats = cv_repeats,
      savePredictions = TRUE
    )
    
    mod_knn <- caret::train(
      Type ~ .,
      data = train,
      method = "knn",
      trControl = final_control
    )
    
    saveRDS(mod_knn, file.path(out_dir, "mod_knn.rds"))
    
    eval_result <- evaluate_model(
      model = mod_knn,
      test_data = test,
      model_label = "KNN",
      out_dir = out_dir,
      roc_file = "ROC_KNN.pdf",
      calibration_file = "Calibration_Curve.pdf",
      dca_file = "Decision_Curve_KNN.pdf",
      metrics_file = "KNN_Performance_Results.csv",
      set_colors = set_colors,
      positive_class = positive_class
    )
    
    list(
      exp = gene_expr,
      model = mod_knn,
      roc = eval_result$roc,
      dca = eval_result$dca,
      explainer = eval_result$explainer,
      model_performance = eval_result$model_performance,
      metrics = eval_result$metrics,
      gene_importance = varimpdf_top
    )
    
  }, error = function(e) {
    message("Error in KNN gene selection: ", e$message)
    NULL
  })
  
} else {
  message("KNN未选择运行。")
  results$knn_gene <- NULL
}

# 23. NNET模型
if (should_run_model("nnet")) {
  
  message("NNET运行中···")
  
  results$nnet_gene <- tryCatch({
    
    out_dir <- "NNET"
    if (!dir.exists(out_dir)) {
      dir.create(out_dir)
    }
    
    group1 <- feature_data$Type
    data0 <- feature_data[, setdiff(colnames(feature_data), "Type"), drop = FALSE]
    
    fit_control <- caret::trainControl(
      method = "repeatedcv",
      number = cv_folds,
      repeats = cv_repeats
    )
    
    nnet_feature_model <- caret::train(
      x = data0,
      y = group1,
      method = "nnet",
      trControl = fit_control,
      trace = FALSE
    )
    
    importance <- caret::varImp(nnet_feature_model, scale = TRUE)
    importance_df <- as.data.frame(importance$importance)
    
    varimpdf_top <- data.frame(
      var = rownames(importance_df),
      impor = importance_df[, 1]
    )
    
    varimpdf_top <- varimpdf_top[order(-varimpdf_top$impor), , drop = FALSE]
    varimpdf_top <- head(varimpdf_top, number_gene)
    
    write.csv(
      varimpdf_top,
      file.path(out_dir, "基因重要度数值表.csv"),
      row.names = FALSE
    )
    
    save_importance_barplot(
      importance_data = varimpdf_top,
      plot_title = "Gene Importance (Neural Network)",
      out_file = file.path(out_dir, "NNET重要度排序.pdf")
    )
    
    nnet_genes <- varimpdf_top$var
    
    write.table(
      nnet_genes,
      file = file.path(out_dir, "NNET特征基因.gene.csv"),
      sep = ",",
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE
    )
    
    gene_expr <- t(raw_expr)
    gene_expr <- auto_log2_transform(gene_expr, do_log2 = do_log2)
    gene_expr <- gene_expr[nnet_genes, , drop = FALSE]
    gene_expr <- as.data.frame(gene_expr)
    
    write.table(
      gene_expr,
      file = file.path(out_dir, "NNET特征基因表达量.csv"),
      sep = ",",
      quote = FALSE,
      row.names = TRUE,
      col.names = TRUE
    )
    
    nnet_model_data <- as.data.frame(t(gene_expr))
    nnet_model_data$Type <- sample_group
    
    train <- nnet_model_data[train_index, , drop = FALSE]
    test <- nnet_model_data[-train_index, , drop = FALSE]
    
    final_control <- caret::trainControl(
      method = "repeatedcv",
      number = cv_folds,
      repeats = cv_repeats,
      savePredictions = TRUE
    )
    
    mod_nnet <- caret::train(
      Type ~ .,
      data = train,
      method = "nnet",
      trControl = final_control,
      trace = FALSE
    )
    
    saveRDS(mod_nnet, file.path(out_dir, "mod_nnet.rds"))
    
    eval_result <- evaluate_model(
      model = mod_nnet,
      test_data = test,
      model_label = "NNET",
      out_dir = out_dir,
      roc_file = "ROC_NNET.pdf",
      calibration_file = "Calibration_Curve.pdf",
      dca_file = "Decision_Curve_NNET.pdf",
      metrics_file = "NNET_Performance_Results.csv",
      set_colors = set_colors,
      positive_class = positive_class
    )
    
    list(
      exp = gene_expr,
      model = mod_nnet,
      roc = eval_result$roc,
      dca = eval_result$dca,
      explainer = eval_result$explainer,
      model_performance = eval_result$model_performance,
      metrics = eval_result$metrics,
      gene_importance = varimpdf_top
    )
    
  }, error = function(e) {
    message("Error in NNET gene selection: ", e$message)
    NULL
  })
  
} else {
  message("NNET未选择运行。")
  results$nnet_gene <- NULL
}

# 24. GLM模型
if (should_run_model("glm")) {
  
  message("GLM运行中···")
  
  results$glm_gene <- tryCatch({
    
    out_dir <- "GLM"
    if (!dir.exists(out_dir)) {
      dir.create(out_dir)
    }
    
    group1 <- feature_data$Type
    data0 <- feature_data[, setdiff(colnames(feature_data), "Type"), drop = FALSE]
    
    fit_control <- caret::trainControl(
      method = "repeatedcv",
      number = cv_folds,
      repeats = cv_repeats
    )
    
    glm_feature_model <- caret::train(
      x = data0,
      y = group1,
      method = "glm",
      family = "binomial",
      trControl = fit_control
    )
    
    importance <- caret::varImp(glm_feature_model, scale = TRUE)
    importance_df <- as.data.frame(importance$importance)
    
    varimpdf_top <- data.frame(
      var = rownames(importance_df),
      impor = importance_df[, 1]
    )
    
    varimpdf_top <- varimpdf_top[order(-varimpdf_top$impor), , drop = FALSE]
    varimpdf_top <- head(varimpdf_top, number_gene)
    
    write.csv(
      varimpdf_top,
      file.path(out_dir, "基因重要度数值表.csv"),
      row.names = FALSE
    )
    
    save_importance_barplot(
      importance_data = varimpdf_top,
      plot_title = "Gene Importance(GLM)",
      out_file = file.path(out_dir, "GLM重要度排序.pdf")
    )
    
    glm_genes <- varimpdf_top$var
    
    write.table(
      glm_genes,
      file = file.path(out_dir, "GLM特征基因.gene.csv"),
      sep = ",",
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE
    )
    
    gene_expr <- t(raw_expr)
    gene_expr <- auto_log2_transform(gene_expr, do_log2 = do_log2)
    gene_expr <- gene_expr[glm_genes, , drop = FALSE]
    gene_expr <- as.data.frame(gene_expr)
    
    write.table(
      gene_expr,
      file = file.path(out_dir, "GLM特征基因表达量.csv"),
      sep = ",",
      quote = FALSE,
      row.names = TRUE,
      col.names = TRUE
    )
    
    glm_model_data <- as.data.frame(t(gene_expr))
    glm_model_data$Type <- sample_group
    
    train <- glm_model_data[train_index, , drop = FALSE]
    test <- glm_model_data[-train_index, , drop = FALSE]
    
    final_control <- caret::trainControl(
      method = "repeatedcv",
      number = cv_folds,
      repeats = cv_repeats,
      savePredictions = TRUE
    )
    
    mod_glm <- caret::train(
      Type ~ .,
      data = train,
      method = "glm",
      family = "binomial",
      trControl = final_control
    )
    
    saveRDS(mod_glm, file.path(out_dir, "mod_glm.rds"))
    
    eval_result <- evaluate_model(
      model = mod_glm,
      test_data = test,
      model_label = "GLM",
      out_dir = out_dir,
      roc_file = "ROC_GLM.pdf",
      calibration_file = "Calibration_Curve.pdf",
      dca_file = "Decision_Curve_GLM.pdf",
      metrics_file = "GLM_Performance_Results.csv",
      set_colors = set_colors,
      positive_class = positive_class
    )
    
    list(
      exp = gene_expr,
      model = mod_glm,
      roc = eval_result$roc,
      dca = eval_result$dca,
      explainer = eval_result$explainer,
      model_performance = eval_result$model_performance,
      metrics = eval_result$metrics,
      gene_importance = varimpdf_top
    )
    
  }, error = function(e) {
    message("Error in GLM gene selection: ", e$message)
    NULL
  })
  
} else {
  message("GLM未选择运行。")
  results$glm_gene <- NULL
}

# 25. C5.0模型
if (should_run_model("c50")) {
  
  message("C5.0运行中···")
  
  results$c50_gene <- tryCatch({
    
    out_dir <- "C5.0"
    if (!dir.exists(out_dir)) {
      dir.create(out_dir)
    }
    
    group1 <- feature_data$Type
    data0 <- feature_data[, setdiff(colnames(feature_data), "Type"), drop = FALSE]
    
    c50_tune_grid <- expand.grid(
      .trials = c(1, 2, 5, 10),
      .model = c("tree", "rules"),
      .winnow = c(TRUE, FALSE)
    )
    
    fit_control <- caret::trainControl(
      method = "repeatedcv",
      number = cv_folds,
      repeats = cv_repeats
    )
    
    c50_feature_model <- caret::train(
      x = data0,
      y = group1,
      method = "C5.0",
      trControl = fit_control,
      tuneGrid = c50_tune_grid
    )
    
    importance <- caret::varImp(c50_feature_model, scale = TRUE)
    importance_df <- as.data.frame(importance$importance)
    
    varimpdf_top <- data.frame(
      var = rownames(importance_df),
      impor = importance_df[, 1]
    )
    
    varimpdf_top <- varimpdf_top[order(-varimpdf_top$impor), , drop = FALSE]
    varimpdf_top <- head(varimpdf_top, number_gene)
    
    write.csv(
      varimpdf_top,
      file.path(out_dir, "基因重要度数值表.csv"),
      row.names = FALSE
    )
    
    save_importance_barplot(
      importance_data = varimpdf_top,
      plot_title = "Gene Importance (C5.0)",
      out_file = file.path(out_dir, "C5.0重要度排序.pdf")
    )
    
    c50_genes <- varimpdf_top$var
    
    write.table(
      c50_genes,
      file = file.path(out_dir, "C5.0特征基因.gene.csv"),
      sep = ",",
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE
    )
    
    gene_expr <- t(raw_expr)
    gene_expr <- auto_log2_transform(gene_expr, do_log2 = do_log2)
    gene_expr <- gene_expr[c50_genes, , drop = FALSE]
    gene_expr <- as.data.frame(gene_expr)
    
    write.table(
      gene_expr,
      file = file.path(out_dir, "C5.0特征基因表达量.csv"),
      sep = ",",
      quote = FALSE,
      row.names = TRUE,
      col.names = TRUE
    )
    
    c50_model_data <- as.data.frame(t(gene_expr))
    c50_model_data$Type <- sample_group
    
    train <- c50_model_data[train_index, , drop = FALSE]
    test <- c50_model_data[-train_index, , drop = FALSE]
    
    final_control <- caret::trainControl(
      method = "repeatedcv",
      number = cv_folds,
      repeats = cv_repeats,
      savePredictions = TRUE
    )
    
    mod_c50 <- caret::train(
      Type ~ .,
      data = train,
      method = "C5.0",
      trControl = final_control,
      tuneGrid = c50_tune_grid
    )
    
    saveRDS(mod_c50, file.path(out_dir, "mod_C50.rds"))
    
    eval_result <- evaluate_model(
      model = mod_c50,
      test_data = test,
      model_label = "C5.0",
      out_dir = out_dir,
      roc_file = "ROC_C50.pdf",
      calibration_file = "Calibration_Curve.pdf",
      dca_file = "Decision_Curve_C50.pdf",
      metrics_file = "C50_Performance_Results.csv",
      set_colors = set_colors,
      positive_class = positive_class
    )
    
    list(
      exp = gene_expr,
      model = mod_c50,
      roc = eval_result$roc,
      dca = eval_result$dca,
      explainer = eval_result$explainer,
      model_performance = eval_result$model_performance,
      metrics = eval_result$metrics,
      gene_importance = varimpdf_top
    )
    
  }, error = function(e) {
    message("Error in C5.0 gene selection: ", e$message)
    NULL
  })
  
} else {
  message("C5.0未选择运行。")
  results$c50_gene <- NULL
}

# 26. 释放并行资源
stopImplicitCluster()

# 27. 输出结果
message("全部已选择模型运行结束。")
results

# 28. 组合作图：根据 run_models 自动选择需要展示的模型
model_name_map <- c(
  lasso = "LASSO",
  rf = "RF",
  svm = "SVM",
  gbm = "GBM",
  xgb = "XGB",
  dt = "DT",
  knn = "KNN",
  nnet = "NNET",
  glm = "GLM",
  c50 = "C50"
)

selected_model_names <- unname(model_name_map[run_models])
selected_model_names <- selected_model_names[!is.na(selected_model_names)]

# 29. 组合作图：提取已经成功运行的模型结果
model_performance_list <- list()
roc_list <- list()
dca_list <- list()
final_model_names <- character(0)

if ("RF" %in% selected_model_names && !is.null(results$rf_gene)) {
  model_performance_list$RF <- results$rf_gene$model_performance
  roc_list$RF <- results$rf_gene$roc
  dca_list$RF <- results$rf_gene$dca
  final_model_names <- c(final_model_names, "RF")
}

if ("SVM" %in% selected_model_names && !is.null(results$svm_gene)) {
  model_performance_list$SVM <- results$svm_gene$model_performance
  roc_list$SVM <- results$svm_gene$roc
  dca_list$SVM <- results$svm_gene$dca
  final_model_names <- c(final_model_names, "SVM")
}

if ("XGB" %in% selected_model_names && !is.null(results$xgb_gene)) {
  model_performance_list$XGB <- results$xgb_gene$model_performance
  roc_list$XGB <- results$xgb_gene$roc
  dca_list$XGB <- results$xgb_gene$dca
  final_model_names <- c(final_model_names, "XGB")
}

if ("GLM" %in% selected_model_names && !is.null(results$glm_gene)) {
  model_performance_list$GLM <- results$glm_gene$model_performance
  roc_list$GLM <- results$glm_gene$roc
  dca_list$GLM <- results$glm_gene$dca
  final_model_names <- c(final_model_names, "GLM")
}

if ("GBM" %in% selected_model_names && !is.null(results$gbm_gene)) {
  model_performance_list$GBM <- results$gbm_gene$model_performance
  roc_list$GBM <- results$gbm_gene$roc
  dca_list$GBM <- results$gbm_gene$dca
  final_model_names <- c(final_model_names, "GBM")
}

if ("KNN" %in% selected_model_names && !is.null(results$knn_gene)) {
  model_performance_list$KNN <- results$knn_gene$model_performance
  roc_list$KNN <- results$knn_gene$roc
  dca_list$KNN <- results$knn_gene$dca
  final_model_names <- c(final_model_names, "KNN")
}

if ("NNET" %in% selected_model_names && !is.null(results$nnet_gene)) {
  model_performance_list$NNET <- results$nnet_gene$model_performance
  roc_list$NNET <- results$nnet_gene$roc
  dca_list$NNET <- results$nnet_gene$dca
  final_model_names <- c(final_model_names, "NNET")
}

if ("LASSO" %in% selected_model_names && !is.null(results$lasso_gene)) {
  model_performance_list$LASSO <- results$lasso_gene$model_performance
  roc_list$LASSO <- results$lasso_gene$roc
  dca_list$LASSO <- results$lasso_gene$dca
  final_model_names <- c(final_model_names, "LASSO")
}

if ("DT" %in% selected_model_names && !is.null(results$dt_gene)) {
  model_performance_list$DT <- results$dt_gene$model_performance
  roc_list$DT <- results$dt_gene$roc
  dca_list$DT <- results$dt_gene$dca
  final_model_names <- c(final_model_names, "DT")
}

if ("C50" %in% selected_model_names && !is.null(results$c50_gene)) {
  model_performance_list$C50 <- results$c50_gene$model_performance
  roc_list$C50 <- results$c50_gene$roc
  dca_list$C50 <- results$c50_gene$dca
  final_model_names <- c(final_model_names, "C50")
}

if (length(final_model_names) == 0) {
  stop("错误：没有可用于组合作图的模型结果。")
}

# 30. 组合作图：残差反向累计分布图
pdf(file = "residual.pdf", width = 6, height = 6)
p1 <- plot(model_performance_list)
print(p1)
dev.off()

# 31. 组合作图：残差箱线图
pdf(file = "boxplot.pdf", width = 6, height = 6)
p2 <- plot(model_performance_list, geom = "boxplot")
print(p2)
dev.off()

# 32. 组合作图：组合ROC曲线
pdf(file = "ROC.pdf", width = 6, height = 6)

set.seed(123)
roc_colors <- rainbow(length(roc_list))

for (i in seq_along(roc_list)) {
  plot(
    roc_list[[i]],
    print.auc = FALSE,
    legacy.axes = TRUE,
    main = "",
    col = roc_colors[i],
    add = i != 1
  )
}

legend_labels <- character(length(roc_list))

for (i in seq_along(roc_list)) {
  legend_labels[i] <- paste0(
    names(roc_list)[i],
    ": ",
    sprintf("%.03f", roc_list[[i]]$auc)
  )
}

legend(
  "bottomright",
  legend_labels,
  col = roc_colors,
  lwd = 2,
  bty = "n"
)

dev.off()

# 33. 组合作图：组合决策曲线
pdf(file = "Decision_Curve.pdf", width = 8, height = 6)

rmda::plot_decision_curve(
  dca_list,
  curve.names = final_model_names,
  xlab = "Threshold Probability",
  confidence.intervals = FALSE,
  ylab = "Benefit",
  col = set_colors,
  lwd = 2,
  legend.position = "left"
)

dev.off()

message("组合作图完成，已生成 residual.pdf、boxplot.pdf、ROC.pdf 和 Decision_Curve.pdf。")

# 34. 加载SHAP分析需要的R包
suppressPackageStartupMessages({
  library(ggplot2)
  library(DALEX)
  library(kernelshap)
  library(shapviz)
})

# 35. 设置SHAP分析参数
Model_name <- "nnet"
xbgsm <- 50

barplot_width <- 6
barplot_height <- 6

bee_width <- 7
bee_height <- 6

dep_width <- 9
dep_height <- 6

water_width <- 7
water_height <- 5

force_width <- 9
force_height <- 5

# 36. 设置可选模型名称
available_models <- c(
  "c50", "dt", "gbm", "glm", "knn",
  "lasso", "nnet", "rf", "svm", "xgb"
)

Model_name <- tolower(Model_name)

if (!Model_name %in% available_models) {
  stop("错误：Model_name 不在可选模型列表中。")
}

if (!should_run_model(Model_name)) {
  stop(
    paste0(
      "错误：当前 SHAP 指定模型 Model_name = '",
      Model_name,
      "'，但该模型没有在 run_models 中运行。请修改 Model_name 或 run_models。"
    )
  )
}

# 37. 检查必要对象
if (!exists("results")) {
  stop("错误：没有找到 results 对象，请先运行前面的模型训练代码。")
}

if (!exists("sample_group")) {
  stop("错误：没有找到 sample_group 对象。")
}

if (!exists("train_index")) {
  stop("错误：没有找到 train_index 对象。")
}

# 38. 生成模型结果名称
Model_exp <- paste0(Model_name, "_gene")

if (!Model_exp %in% names(results)) {
  stop(paste0("错误：results 中没有找到 ", Model_exp, "。"))
}

if (is.null(results[[Model_exp]])) {
  stop(paste0("错误：", Model_exp, " 的结果为空，可能该模型前面运行失败。"))
}

if (!"exp" %in% names(results[[Model_exp]])) {
  stop(paste0("错误：", Model_exp, " 中没有 exp 特征表达矩阵。"))
}

if (!"explainer" %in% names(results[[Model_exp]])) {
  stop(paste0("错误：", Model_exp, " 中没有 explainer 对象。"))
}

# 39. 提取特征基因表达矩阵和解释器对象
exp <- results[[Model_exp]]$exp
explainer_obj <- results[[Model_exp]]$explainer

if (is.null(exp)) {
  stop("错误：特征表达矩阵 exp 为空。")
}

if (is.null(explainer_obj)) {
  stop("错误：explainer_obj 为空。")
}

# 40. 整理SHAP输入数据
data1 <- t(exp)
data1 <- as.data.frame(data1)
data1$Type <- sample_group

if (length(sample_group) != nrow(data1)) {
  stop("错误：sample_group 长度和SHAP输入表达矩阵的样本数不一致。")
}

# 41. 按前面建模时的训练集索引划分数据
train <- data1[train_index, , drop = FALSE]
test <- data1[-train_index, , drop = FALSE]

# 42. 创建SHAP结果目录
if (!dir.exists("shap_results")) {
  dir.create("shap_results")
}

if (!dir.exists("shap_results/permshap")) {
  dir.create("shap_results/permshap", recursive = TRUE)
}

if (!dir.exists("shap_results/kernelshap")) {
  dir.create("shap_results/kernelshap", recursive = TRUE)
}

# 43. 运行permshap
message("计算 permshap...")

fit_permshap <- permshap(
  explainer_obj,
  train[, -ncol(train), drop = FALSE],
  bg_n = xbgsm
)

shp_permshap <- shapviz(
  fit_permshap,
  X = train[, -ncol(train), drop = FALSE]
)

# 44. 输出permshap barplot
pdf(
  file = paste0("shap_results/permshap/", Model_name, "_barplot.pdf"),
  width = barplot_width,
  height = barplot_height
)

print(
  sv_importance(
    shp_permshap,
    kind = "bar",
    show_numbers = TRUE
  ) +
    theme_bw()
)

dev.off()

# 45. 输出permshap beeplot
pdf(
  file = paste0("shap_results/permshap/", Model_name, "_bee.pdf"),
  width = bee_width,
  height = bee_height
)

print(
  sv_importance(
    shp_permshap,
    kind = "bee",
    show_numbers = TRUE
  ) +
    theme_bw()
)

dev.off()

# 46. 输出permshap dependence plot
pdf(
  file = paste0("shap_results/permshap/", Model_name, "_dependence.pdf"),
  width = dep_width,
  height = dep_height
)

print(
  sv_dependence(
    shp_permshap,
    v = names(sort(colMeans(abs(shp_permshap$S)), decreasing = TRUE))
  ) +
    theme_bw()
)

dev.off()

# 47. 输出permshap waterfall plot
pdf(
  file = paste0("shap_results/permshap/", Model_name, "_waterfall.pdf"),
  width = water_width,
  height = water_height
)

print(
  sv_waterfall(
    shp_permshap,
    row_id = 1
  )
)

dev.off()

# 48. 输出permshap force plot
pdf(
  file = paste0("shap_results/permshap/", Model_name, "_force.pdf"),
  width = force_width,
  height = force_height
)

print(
  sv_force(
    shp_permshap,
    row_id = 1
  )
)

dev.off()

message("Permshap 分析已完成。")

# 49. 运行kernelshap
message("计算 kernelshap...")

fit_kernelshap <- kernelshap(
  explainer_obj,
  train[, -ncol(train), drop = FALSE],
  bg_n = xbgsm
)

shp_kernelshap <- shapviz(
  fit_kernelshap,
  X = train[, -ncol(train), drop = FALSE]
)

# 50. 输出kernelshap barplot
pdf(
  file = paste0("shap_results/kernelshap/", Model_name, "_barplot.pdf"),
  width = barplot_width,
  height = barplot_height
)

print(
  sv_importance(
    shp_kernelshap,
    kind = "bar",
    show_numbers = TRUE
  ) +
    theme_bw()
)

dev.off()

# 51. 输出kernelshap beeplot
pdf(
  file = paste0("shap_results/kernelshap/", Model_name, "_bee.pdf"),
  width = bee_width,
  height = bee_height
)

print(
  sv_importance(
    shp_kernelshap,
    kind = "bee",
    show_numbers = TRUE
  ) +
    theme_bw()
)

dev.off()

# 52. 输出kernelshap dependence plot
pdf(
  file = paste0("shap_results/kernelshap/", Model_name, "_dependence.pdf"),
  width = dep_width,
  height = dep_height
)

print(
  sv_dependence(
    shp_kernelshap,
    v = names(sort(colMeans(abs(shp_kernelshap$S)), decreasing = TRUE))
  ) +
    theme_bw()
)

dev.off()

# 53. 输出kernelshap waterfall plot
pdf(
  file = paste0("shap_results/kernelshap/", Model_name, "_waterfall.pdf"),
  width = water_width,
  height = water_height
)

print(
  sv_waterfall(
    shp_kernelshap,
    row_id = 1
  )
)

dev.off()

# 54. 输出kernelshap force plot
pdf(
  file = paste0("shap_results/kernelshap/", Model_name, "_force.pdf"),
  width = force_width,
  height = force_height
)

print(
  sv_force(
    shp_kernelshap,
    row_id = 1
  )
)

dev.off()

# 55. 输出完成提示
message("Kernelshap 分析已完成。")
message("SHAP分析完成，结果已保存到 shap_results 文件夹。")