# 1. 加载必要 R 包

suppressPackageStartupMessages({
  library(IOBR)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(reshape2)
  library(ggpubr)
  library(pheatmap)
})

# 2. 设置表达矩阵文件

expr_file_path <- "tpms.csv"

# 3. 设置分组信息

group_names <- c("Stage1", "Stage2", "Stage3")
group_sizes <- c(5, 8, 7)

# 多组示例：
# group_names <- c("Stage1", "Stage2", "Stage3", "Stage4")
# group_sizes <- c(5, 5, 5, 5)

# 4. 设置表达矩阵预处理参数

iflog <- TRUE
ifmad <- FALSE
top_nim <- 10000

# 5. 设置免疫浸润分析方法

# 1 = mcpcounter
# 2 = epic
# 3 = xcell
# 4 = ips
# 5 = estimate
# 6 = quantiseq
# 7 = timer
# 8 = cibersort
# 9 = cibersort_abs

method <- c(1, 3, 5)

# 6. 设置免疫浸润方法参数

tumor_type <- "lusc"
tumor <- FALSE
scale_mrna <- TRUE
arrays <- FALSE
perm <- 1000
platform <- "affymetrix"
P_value_CIBERSORT <- 0.05

# 7. 设置需要做免疫相关性的基因

gene.duo <- c("TP53", "DPM1", "ABCEFD")

# 8. 设置图1：分面箱线图参数

values778 <- c("#4DBBD5FF", "#E64B35FF")

# compare_mode 可选：
# "overall"  = 只显示整体组间差异检验
# "pairwise" = 只显示自动两两比较
# "both"     = 同时显示整体检验和两两比较

compare_mode <- "both"

overall_test_method <- "auto"
pairwise_test_method <- "wilcox.test"

overall_p_label <- "p.format"
pairwise_p_label <- "p.signif"

hide_ns <- FALSE

overall_label_y_npc <- 0.95
pairwise_step_increase <- 0.08

size.fig1 <- 8
colour.fig1 <- "black"
angle.fig1x <- 45
width.fig1 <- 16
height.fig1 <- 12
facet_ncol <- 3

# 9. 设置图2：单样本比例堆叠图参数

width.figt2 <- 0.7
linewidth.fig2 <- 0.5
colour.fig2 <- "#222222"
color.fig2 <- "black"
linewidth.figw2 <- 0.5
width.fig2 <- 16
height.fig2 <- 6

# 10. 设置图3：多基因免疫相关性热图参数

method.fig3 <- "spearman"
cluster1 <- FALSE
color.fig3 <- c("steelblue", "white", "firebrick")
fontsize <- 9
cellheight <- 30
cellwidth <- 30
width.fig3 <- 20
height.fig3 <- 20

# 11. 设置图4：单基因免疫相关性气泡图参数

low.fig4 <- "blue"
mid.fig4 <- "white"
high.fig4 <- "red"
width.fig4 <- 8
height.fig4 <- 8

# 12. 创建文件夹函数

create_dir_if_needed <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
}

# 13. 生成分组颜色函数

get_group_colors <- function(group_names, base_colors) {
  group_number <- length(group_names)
  
  if (length(base_colors) >= group_number) {
    final_colors <- base_colors[seq_len(group_number)]
  } else {
    final_colors <- grDevices::hcl.colors(
      n = group_number,
      palette = "Dark 3"
    )
  }
  
  names(final_colors) <- group_names
  
  return(final_colors)
}

# 14. 自动选择整体统计检验方法

get_group_test_method <- function(group_vector, preferred_method = "auto") {
  group_vector <- droplevels(factor(group_vector))
  group_number <- nlevels(group_vector)
  
  if (group_number < 2) {
    stop("分组数量少于 2，无法进行组间比较。")
  }
  
  if (preferred_method != "auto") {
    return(preferred_method)
  }
  
  if (group_number == 2) {
    return("wilcox.test")
  } else {
    return("kruskal.test")
  }
}

# 15. 生成两两比较列表

get_pairwise_comparisons <- function(group_names) {
  if (length(group_names) < 2) {
    stop("分组数量少于 2，无法进行两两比较。")
  }
  
  comparisons_list <- combn(
    group_names,
    2,
    simplify = FALSE
  )
  
  return(comparisons_list)
}

# 16. 读取表达矩阵

rtim <- read.csv(
  expr_file_path,
  header = TRUE,
  row.names = 1,
  check.names = TRUE
)

# 17. 生成分组信息

if (length(group_names) != length(group_sizes)) {
  stop("group_names 和 group_sizes 的长度不一致。")
}

group_labels <- rep(
  group_names,
  times = group_sizes
)

if (length(group_labels) != ncol(rtim)) {
  stop("分组数量和表达矩阵样本数量不一致，请检查 group_names、group_sizes 和表达矩阵列数。")
}

new_db <- data.frame(
  Group = factor(group_labels, levels = group_names),
  stringsAsFactors = FALSE
)

rownames(new_db) <- colnames(rtim)

group_colors <- get_group_colors(
  group_names = group_names,
  base_colors = values778
)

write.csv(
  new_db,
  file = "分组匹配信息(小提琴图用).csv",
  row.names = TRUE
)

# 18. 表达矩阵 log2 转换判断

rt <- rtim
ex <- rt

if (iflog) {
  qx <- as.numeric(
    quantile(
      ex,
      c(0.00, 0.25, 0.5, 0.75, 0.99, 1.0),
      na.rm = TRUE
    )
  )
  
  LogC <- (qx[5] > 100) ||
    (qx[6] - qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
  
  if (LogC) {
    rt <- log2(ex + 1)
  } else {
    rt <- ex
  }
} else {
  rt <- ex
}

# 19. MAD 筛选

if (ifmad) {
  mad_values <- apply(
    rt,
    1,
    function(x) median(abs(x - median(x)))
  )
  
  top_n_use <- min(top_nim, nrow(rt))
  
  expr <- rt[
    order(mad_values, decreasing = TRUE)[seq_len(top_n_use)],
    ,
    drop = FALSE
  ]
} else {
  expr <- rt
}

exprim <- expr

# 20. 通用函数：免疫浸润结果作图

save_immune_result_plots <- function(score_data,
                                     out_dir,
                                     return_name,
                                     expr_for_correlation = exprim,
                                     cell_name_clean = FALSE,
                                     box_id_replace_hyphen = FALSE,
                                     stack_id_replace_hyphen = FALSE) {
  
  create_dir_if_needed(out_dir)
  
  cluster <- cluster1
  
  data <- as.data.frame(score_data, check.names = FALSE)
  id_col <- colnames(data)[1]
  
  if (box_id_replace_hyphen) {
    data[[id_col]] <- gsub("-", ".", data[[id_col]])
  }
  
  data$Group <- new_db[data[[id_col]], "Group"]
  data$Group <- factor(data$Group, levels = group_names)
  
  if (anyNA(data$Group)) {
    bad_id <- data[[id_col]][is.na(data$Group)]
    stop(
      paste0(
        out_dir,
        " 中有样本没有匹配到分组：",
        paste(head(bad_id, 20), collapse = ", ")
      )
    )
  }
  
  box_data <- data[, -1, drop = FALSE]
  
  df <- reshape2::melt(
    box_data,
    id.vars = "Group",
    variable.name = "Cell",
    value.name = "score"
  )
  
  if (cell_name_clean) {
    df$Cell <- sub("_[^_]*$", "", df$Cell)
    df$Cell <- gsub("_", " ", df$Cell)
  }
  
  stat_method <- get_group_test_method(
    group_vector = df$Group,
    preferred_method = overall_test_method
  )
  
  comparisons_list <- get_pairwise_comparisons(group_names)
  
  p <- ggplot(
    data = df,
    aes(x = Group, y = score, fill = Group)
  ) +
    geom_boxplot(
      width = 0.5,
      outlier.colour = NA
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0.05, 0.30))
    ) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    scale_fill_manual(values = group_colors)
  
  if (compare_mode %in% c("overall", "both")) {
    p <- p +
      ggpubr::stat_compare_means(
        method = stat_method,
        label = overall_p_label,
        hide.ns = hide_ns,
        label.y.npc = overall_label_y_npc,
        size = 3
      )
  }
  
  if (compare_mode %in% c("pairwise", "both")) {
    p <- p +
      ggpubr::stat_compare_means(
        comparisons = comparisons_list,
        method = pairwise_test_method,
        label = pairwise_p_label,
        hide.ns = hide_ns,
        step.increase = pairwise_step_increase,
        size = 3
      )
  }
  
  p <- p +
    facet_wrap(
      ~ Cell,
      scales = "free_y",
      ncol = facet_ncol
    ) +
    theme(
      axis.text.x = element_text(
        size = size.fig1,
        colour = colour.fig1,
        angle = angle.fig1x,
        hjust = 1
      ),
      axis.text.y = element_text(size = size.fig1, angle = 0),
      axis.title.x = element_text(size = size.fig1),
      axis.title.y = element_text(size = size.fig1),
      legend.text = element_text(size = size.fig1),
      strip.text = element_text(size = size.fig1)
    ) +
    xlab("") +
    ylab("Immune score")
  
  ggsave(
    filename = file.path(out_dir, "1.immune_boxplot_facet.pdf"),
    plot = p,
    width = width.fig1,
    height = height.fig1
  )
  
  im_timer1 <- as.data.frame(score_data, check.names = FALSE)
  id_col_stack <- colnames(im_timer1)[1]
  
  im_timer1 <- tibble::column_to_rownames(
    im_timer1,
    var = id_col_stack
  )
  
  if (stack_id_replace_hyphen) {
    rownames(im_timer1) <- gsub("-", ".", rownames(im_timer1))
  }
  
  total_scores <- rowSums(im_timer1)
  
  if (any(total_scores == 0, na.rm = TRUE)) {
    stop("存在总评分为 0 的样本，无法计算比例堆叠图。")
  }
  
  proportion_scores <- sweep(im_timer1, 1, total_scores, "/")
  
  timer1 <- as.data.frame(proportion_scores)
  
  timer1 <- timer1 %>%
    tibble::rownames_to_column(var = "Sample")
  
  Cellratio <- timer1 %>%
    tidyr::pivot_longer(
      cols = -Sample,
      names_to = "CellName",
      values_to = "Proportion"
    )
  
  P1 <- ggplot(Cellratio) +
    geom_bar(
      aes(x = Sample, y = Proportion, fill = CellName),
      stat = "identity",
      width = width.figt2,
      linewidth = linewidth.fig2,
      colour = colour.fig2
    ) +
    theme_classic() +
    labs(x = "Sample", y = "Ratio") +
    coord_flip() +
    theme(
      panel.border = element_rect(
        fill = NA,
        color = color.fig2,
        linewidth = linewidth.figw2,
        linetype = "solid"
      )
    )
  
  pdf(
    file = file.path(out_dir, "2.immune_stacke.pdf"),
    width = width.fig2,
    height = height.fig2
  )
  print(P1)
  dev.off()
  
  write.csv(
    proportion_scores,
    file = file.path(out_dir, "proportion_scores.csv"),
    row.names = TRUE
  )
  
  gene <- gene.duo
  gene_keep <- gene[gene %in% rownames(expr_for_correlation)]
  
  if (length(gene_keep) == 0) {
    stop("gene.duo 中没有任何基因存在于表达矩阵行名中。")
  }
  
  if (length(gene_keep) < length(gene)) {
    warning("gene.duo 中部分基因不存在于表达矩阵行名中，已自动跳过。")
  }
  
  gene_exp <- expr_for_correlation[gene_keep, , drop = FALSE]
  gene_exp <- gene_exp[complete.cases(gene_exp), , drop = FALSE]
  
  if (nrow(gene_exp) == 0) {
    stop("目标基因表达矩阵经过 NA 过滤后为空，无法做相关性分析。")
  }
  
  gene_exp <- as.data.frame(t(gene_exp))
  
  common_rows <- intersect(
    rownames(gene_exp),
    rownames(im_timer1)
  )
  
  if (length(common_rows) == 0) {
    stop("表达矩阵样本名和免疫浸润结果样本名没有交集，无法做相关性分析。")
  }
  
  gene_exp_subset <- gene_exp[common_rows, , drop = FALSE]
  im_timer1_subset <- im_timer1[common_rows, , drop = FALSE]
  
  correlation_matrix <- cor(
    gene_exp_subset,
    im_timer1_subset,
    method = method.fig3
  )
  
  write.csv(
    correlation_matrix,
    file = file.path(out_dir, "相关性数值.csv"),
    row.names = TRUE
  )
  
  na_cols <- apply(is.na(correlation_matrix), 2, any)
  correlation_matrix <- correlation_matrix[, !na_cols, drop = FALSE]
  
  if (ncol(correlation_matrix) == 0) {
    stop("相关性矩阵全部为 NA，无法绘制热图。")
  }
  
  datacor <- t(correlation_matrix)
  
  pdf(
    file = file.path(out_dir, "3.多个基因免疫相关性图.pdf"),
    width = width.fig3,
    height = height.fig3
  )
  
  pheatmap::pheatmap(
    datacor,
    scale = "none",
    cluster_cols = cluster,
    cluster_rows = cluster,
    color = colorRampPalette(color.fig3)(50),
    legend = TRUE,
    show_colnames = TRUE,
    show_rownames = TRUE,
    fontsize = fontsize,
    cellheight = cellheight,
    cellwidth = cellwidth
  )
  
  dev.off()
  
  for (siggene in rownames(correlation_matrix)) {
    siggene_cor <- as.data.frame(
      t(correlation_matrix[siggene, , drop = FALSE])
    )
    
    siggene_cor$immune <- rownames(siggene_cor)
    cor_values <- siggene_cor[, 1]
    title_name <- names(siggene_cor)[1]
    
    p3 <- ggplot(
      data = siggene_cor,
      aes(x = cor_values, y = immune)
    ) +
      labs(x = "Correlation coefficient", y = "Immune cell") +
      geom_point(aes(size = abs(cor_values), color = cor_values)) +
      scale_color_gradient2(
        low = low.fig4,
        mid = mid.fig4,
        high = high.fig4,
        midpoint = 0
      ) +
      theme(
        panel.background = element_rect(
          fill = "white",
          linewidth = 1,
          color = "black"
        ),
        panel.grid = element_line(
          color = "grey75",
          linewidth = 0.5
        ),
        axis.ticks = element_line(linewidth = 0.5),
        axis.text.y = element_text(
          colour = "black",
          size = 9
        ),
        plot.title = element_text(
          hjust = 0.5,
          size = 14,
          face = "bold"
        )
      ) +
      xlab("Correlation coefficient") +
      ylab("Immune cell") +
      ggtitle(title_name)
    
    pdf(
      file = file.path(
        out_dir,
        paste0(siggene, "_immune_correlation.pdf")
      ),
      width = width.fig4,
      height = height.fig4
    )
    
    print(p3)
    dev.off()
  }
  
  return_list <- list(
    data = data,
    proportion_scores = proportion_scores,
    correlation_matrix = correlation_matrix
  )
  
  return_list[[return_name]] <- df
  
  return(return_list)
}

# 21. 运行免疫浸润分析

results_im <- list()

if (1 %in% method) {
  create_dir_if_needed("im_mcpcounter")
  
  im_mcpcounter <- IOBR::deconvo_tme(
    eset = exprim,
    method = "mcpcounter"
  )
  
  write.csv(
    im_mcpcounter,
    file = "im_mcpcounter/im_mcpcounter.csv",
    row.names = FALSE
  )
  
  results_im$im_mcpcounter1 <- save_immune_result_plots(
    score_data = im_mcpcounter,
    out_dir = "im_mcpcounter",
    return_name = "mcpcounter",
    expr_for_correlation = exprim,
    cell_name_clean = FALSE,
    box_id_replace_hyphen = FALSE,
    stack_id_replace_hyphen = TRUE
  )
}

if (2 %in% method) {
  create_dir_if_needed("im_epic")
  
  im_epic <- IOBR::deconvo_tme(
    eset = exprim,
    method = "epic",
    tumor = tumor,
    scale_mrna = scale_mrna
  )
  
  write.csv(
    im_epic,
    file = "im_epic/im_epic.csv",
    row.names = FALSE
  )
  
  results_im$im_epic2 <- save_immune_result_plots(
    score_data = im_epic,
    out_dir = "im_epic",
    return_name = "epic",
    expr_for_correlation = exprim,
    cell_name_clean = FALSE,
    box_id_replace_hyphen = FALSE,
    stack_id_replace_hyphen = FALSE
  )
}

if (3 %in% method) {
  create_dir_if_needed("im_xcell")
  
  im_xcell <- suppressWarnings(
    IOBR::deconvo_tme(
      eset = exprim,
      method = "xcell",
      arrays = arrays
    )
  )
  
  write.csv(
    im_xcell,
    file = "im_xcell/im_xcell.csv",
    row.names = FALSE
  )
  
  results_im$im_xcell3 <- save_immune_result_plots(
    score_data = im_xcell,
    out_dir = "im_xcell",
    return_name = "xcell",
    expr_for_correlation = exprim,
    cell_name_clean = FALSE,
    box_id_replace_hyphen = FALSE,
    stack_id_replace_hyphen = TRUE
  )
}

if (4 %in% method) {
  create_dir_if_needed("im_ips")
  
  im_ips <- IOBR::deconvo_tme(
    eset = exprim,
    method = "ips"
  )
  
  write.csv(
    im_ips,
    file = "im_ips/im_ips.csv",
    row.names = FALSE
  )
  
  results_im$im_ips4 <- save_immune_result_plots(
    score_data = im_ips,
    out_dir = "im_ips",
    return_name = "ips",
    expr_for_correlation = exprim,
    cell_name_clean = FALSE,
    box_id_replace_hyphen = FALSE,
    stack_id_replace_hyphen = TRUE
  )
}

if (5 %in% method) {
  create_dir_if_needed("im_estimate")
  
  im_estimate <- IOBR::deconvo_tme(
    eset = exprim,
    method = "estimate",
    platform = platform
  )
  
  write.csv(
    im_estimate,
    file = "im_estimate/im_estimate.csv",
    row.names = FALSE
  )
  
  results_im$im_estimate5 <- save_immune_result_plots(
    score_data = im_estimate,
    out_dir = "im_estimate",
    return_name = "estimate",
    expr_for_correlation = exprim,
    cell_name_clean = FALSE,
    box_id_replace_hyphen = TRUE,
    stack_id_replace_hyphen = TRUE
  )
}

if (6 %in% method) {
  create_dir_if_needed("im_quantiseq")
  
  im_quantiseq <- IOBR::deconvo_tme(
    eset = exprim,
    method = "quantiseq",
    scale_mrna = scale_mrna
  )
  
  write.csv(
    im_quantiseq,
    file = "im_quantiseq/im_quantiseq.csv",
    row.names = FALSE
  )
  
  results_im$im_quantiseq6 <- save_immune_result_plots(
    score_data = im_quantiseq,
    out_dir = "im_quantiseq",
    return_name = "quantiseq",
    expr_for_correlation = exprim,
    cell_name_clean = FALSE,
    box_id_replace_hyphen = FALSE,
    stack_id_replace_hyphen = TRUE
  )
}

if (7 %in% method) {
  create_dir_if_needed("im_timer")
  
  im_timer <- IOBR::deconvo_timer(
    eset = exprim,
    indications = rep(tumor_type, ncol(exprim))
  )
  
  write.csv(
    im_timer,
    file = "im_timer/im_timer免疫评分.csv",
    row.names = FALSE
  )
  
  results_im$im_timer7 <- save_immune_result_plots(
    score_data = im_timer,
    out_dir = "im_timer",
    return_name = "timer",
    expr_for_correlation = exprim,
    cell_name_clean = FALSE,
    box_id_replace_hyphen = FALSE,
    stack_id_replace_hyphen = FALSE
  )
}

if (8 %in% method) {
  create_dir_if_needed("im_cibersort")
  
  expr_cibersort <- exprim
  
  if (max(expr_cibersort, na.rm = TRUE) < 50) {
    expr_cibersort <- 2^expr_cibersort - 1
  }
  
  im_cibersort <- IOBR::deconvo_tme(
    eset = expr_cibersort,
    method = "cibersort",
    arrays = arrays,
    perm = perm
  )
  
  write.csv(
    im_cibersort,
    file = "im_cibersort/im_cibersort免疫评分.csv",
    row.names = FALSE
  )
  
  data_cibersort <- as.data.frame(im_cibersort, check.names = FALSE)
  
  if (!"P-value_CIBERSORT" %in% colnames(data_cibersort)) {
    stop("CIBERSORT 结果中未找到 P-value_CIBERSORT 列。")
  }
  
  data_cibersort <- data_cibersort[
    data_cibersort$`P-value_CIBERSORT` <= P_value_CIBERSORT,
    ,
    drop = FALSE
  ]
  
  if (nrow(data_cibersort) == 0) {
    stop("CIBERSORT 筛选后没有剩余样本，请放宽 P_value_CIBERSORT 或检查输入表达矩阵。")
  }
  
  score_cibersort <- data_cibersort[
    ,
    1:(ncol(data_cibersort) - 3),
    drop = FALSE
  ]
  
  results_im$im_cibersort8 <- save_immune_result_plots(
    score_data = score_cibersort,
    out_dir = "im_cibersort",
    return_name = "cibersort",
    expr_for_correlation = expr_cibersort,
    cell_name_clean = TRUE,
    box_id_replace_hyphen = FALSE,
    stack_id_replace_hyphen = FALSE
  )
}

if (9 %in% method) {
  create_dir_if_needed("im_cibersort_abs")
  
  expr_cibersort_abs <- exprim
  
  if (max(expr_cibersort_abs, na.rm = TRUE) < 50) {
    expr_cibersort_abs <- 2^expr_cibersort_abs - 1
  }
  
  im_cibersort_abs <- IOBR::deconvo_tme(
    eset = expr_cibersort_abs,
    method = "cibersort_abs",
    arrays = arrays,
    perm = perm
  )
  
  write.csv(
    im_cibersort_abs,
    file = "im_cibersort_abs/im_cibersort_abs.csv",
    row.names = FALSE
  )
  
  data_cibersort_abs <- as.data.frame(
    im_cibersort_abs,
    check.names = FALSE
  )
  
  if (!"P-value_CIBERSORT" %in% colnames(data_cibersort_abs)) {
    stop("CIBERSORT-ABS 结果中未找到 P-value_CIBERSORT 列。")
  }
  
  data_cibersort_abs <- data_cibersort_abs[
    data_cibersort_abs$`P-value_CIBERSORT` <= P_value_CIBERSORT,
    ,
    drop = FALSE
  ]
  
  if (nrow(data_cibersort_abs) == 0) {
    stop("CIBERSORT-ABS 筛选后没有剩余样本，请放宽 P_value_CIBERSORT 或检查输入表达矩阵。")
  }
  
  score_cibersort_abs <- data_cibersort_abs[
    ,
    1:(ncol(data_cibersort_abs) - 4),
    drop = FALSE
  ]
  
  results_im$im_cibersort_abs9 <- save_immune_result_plots(
    score_data = score_cibersort_abs,
    out_dir = "im_cibersort_abs",
    return_name = "cibersort_abs",
    expr_for_correlation = expr_cibersort_abs,
    cell_name_clean = TRUE,
    box_id_replace_hyphen = FALSE,
    stack_id_replace_hyphen = FALSE
  )
}

# 22. 输出完成提示

cat("IOBR 免疫浸润分析完成。\n")
cat("预处理后的表达矩阵对象：exprim\n")
cat("分组信息对象：new_db\n")
cat("免疫浸润结果对象：results_im\n")