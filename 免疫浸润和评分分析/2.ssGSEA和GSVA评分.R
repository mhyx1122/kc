# 1. 加载必要 R 包

suppressPackageStartupMessages({
  library(IOBR)
  library(GSVA)
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

# 5. 设置 ssGSEA / GSVA 分析参数

signature_type <- "signature_tme"
mini_gene_count <- 5
use_custom_data <- T
custom_gene_file_path <- "19PCDgenes.csv"

# 6. 设置需要做相关性的基因

gene.duo <- c("TP53", "DPM1", "ABCEFD")

# 7. 设置分面箱线图参数

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

# 8. 设置单样本评分组成堆叠图参数

width.figt2 <- 0.7
linewidth.fig2 <- 0.5
colour.fig2 <- "#222222"
color.fig2 <- "black"
linewidth.figw2 <- 0.5
width.fig2 <- 16
height.fig2 <- 6

# 9. 设置多基因相关性热图参数

method.fig3 <- "spearman"
cluster1 <- FALSE
color.fig3 <- c("steelblue", "white", "firebrick")
fontsize <- 9
cellheight <- 30
cellwidth <- 30
width.fig3 <- 20
height.fig3 <- 20

# 10. 设置单基因相关性气泡图参数

low.fig4 <- "blue"
mid.fig4 <- "white"
high.fig4 <- "red"
width.fig4 <- 8
height.fig4 <- 8

# 11. 设置小提琴图参数

violin_score_file_path <- "ssgsea/ssgsea_order.csv"
violin_group_file_path <- "分组匹配信息(小提琴图用).csv"

violin_colors <- c("#e5451d", "#9084bd")
violin_pdf_width <- 12
violin_pdf_height <- 15
violin_facet_ncol <- 3

plot_points <- TRUE
violin_alpha <- 0.6
scatter_alpha <- 0.6

p_display <- "numeric"

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

# 16. 计算单个指标的整体组间 P 值

calculate_group_pvalue <- function(value_vector,
                                   group_vector,
                                   preferred_method = "auto") {
  group_vector <- droplevels(factor(group_vector))
  value_vector <- as.numeric(value_vector)
  
  valid_index <- !is.na(value_vector) & !is.na(group_vector)
  value_vector <- value_vector[valid_index]
  group_vector <- droplevels(group_vector[valid_index])
  
  if (nlevels(group_vector) < 2) {
    return(NA_real_)
  }
  
  if (length(unique(value_vector)) <= 1) {
    return(NA_real_)
  }
  
  test_method <- get_group_test_method(
    group_vector = group_vector,
    preferred_method = preferred_method
  )
  
  p_value <- tryCatch(
    {
      if (test_method == "wilcox.test") {
        wilcox.test(value_vector ~ group_vector)$p.value
      } else if (test_method == "t.test") {
        t.test(value_vector ~ group_vector)$p.value
      } else if (test_method == "kruskal.test") {
        kruskal.test(value_vector ~ group_vector)$p.value
      } else if (test_method == "anova") {
        summary(aov(value_vector ~ group_vector))[[1]]$"Pr(>F)"[1]
      } else {
        stop("不支持的统计检验方法。")
      }
    },
    error = function(e) {
      NA_real_
    }
  )
  
  return(p_value)
}

# 17. 给 ggplot 添加统计检验图层

add_compare_layers <- function(p,
                               group_vector,
                               compare_mode,
                               group_names) {
  
  if (!compare_mode %in% c("overall", "pairwise", "both")) {
    stop("compare_mode 只能是 overall、pairwise 或 both。")
  }
  
  stat_method <- get_group_test_method(
    group_vector = group_vector,
    preferred_method = overall_test_method
  )
  
  comparisons_list <- get_pairwise_comparisons(group_names)
  
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
  
  return(p)
}

# 18. 读取表达矩阵

rtim <- read.csv(
  expr_file_path,
  header = TRUE,
  row.names = 1,
  check.names = TRUE
)

# 19. 生成分组信息

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

violin_colors_use <- get_group_colors(
  group_names = group_names,
  base_colors = violin_colors
)

write.csv(
  new_db,
  file = "分组匹配信息(小提琴图用).csv",
  row.names = TRUE
)

# 20. 表达矩阵 log2 转换判断

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

# 21. MAD 筛选

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

# 22. 创建 ssGSEA / GSVA 结果目录

create_dir_if_needed("ssgsea")

# 23. 选择内置基因集或自定义基因集

if (!use_custom_data) {
  if (signature_type == "signature_tme") {
    signatureAPP <- signature_tme
  } else if (signature_type == "signature_metabolism") {
    signatureAPP <- signature_metabolism
  } else if (signature_type == "go_bp") {
    signatureAPP <- go_bp
  } else if (signature_type == "go_mf") {
    signatureAPP <- go_mf
  } else if (signature_type == "go_cc") {
    signatureAPP <- go_cc
  } else if (signature_type == "kegg") {
    signatureAPP <- kegg
  } else if (signature_type == "hallmark") {
    signatureAPP <- hallmark
  } else {
    stop("signature_type 设置错误。")
  }
}

if (use_custom_data) {
  gene_data <- read.csv(
    custom_gene_file_path,
    header = TRUE,
    row.names = NULL,
    check.names = FALSE
  )
  
  PD_list <- list()
  
  for (col_name in colnames(gene_data)) {
    genes <- trimws(gene_data[[col_name]])
    genes <- genes[genes != ""]
    genes <- genes[!is.na(genes)]
    PD_list[[col_name]] <- genes
  }
  
  signatureAPP <- PD_list
}

signature <- signatureAPP
expr <- exprim

# 24. 运行 ssGSEA

ssgsea_param <- GSVA::ssgseaParam(
  exprData = as.matrix(expr),
  geneSets = signature,
  minSize = mini_gene_count
)

ssgsea_result <- GSVA::gsva(
  ssgsea_param,
  verbose = FALSE
)

ssgsearesult <- t(as.data.frame(ssgsea_result))

write.csv(
  ssgsearesult,
  file = "ssgsea/ssgsea_score.csv",
  row.names = TRUE
)

# 25. 运行 GSVA

gsva_params <- GSVA::gsvaParam(
  exprData = as.matrix(expr),
  geneSets = signature,
  kcdf = "Gaussian"
)

gsva_result <- GSVA::gsva(
  gsva_params,
  verbose = FALSE
)

gsvaresult <- t(as.data.frame(gsva_result))

write.csv(
  gsvaresult,
  file = "ssgsea/gsva_score.csv",
  row.names = TRUE
)

# 26. 保存基因集参考文献

if (exists("signature_collection_citation")) {
  write.csv(
    signature_collection_citation,
    file = "ssgsea/ssgsea所用基因集的参考文献.csv",
    row.names = FALSE
  )
}

# 27. 通用函数：按整体组间差异 P 值排序

sort_score_by_pvalue <- function(input_file,
                                 output_file,
                                 group_object,
                                 test_method = overall_test_method,
                                 result_name = "score") {
  
  data_DEG <- read.csv(
    input_file,
    header = TRUE,
    row.names = 1,
    check.names = FALSE
  )
  
  data_DEG$Group <- group_object[rownames(data_DEG), "Group"]
  data_DEG$Group <- factor(data_DEG$Group, levels = group_names)
  
  if (anyNA(data_DEG$Group)) {
    bad_id <- rownames(data_DEG)[is.na(data_DEG$Group)]
    stop(
      paste0(
        result_name,
        " 结果中有样本没有匹配到分组：",
        paste(head(bad_id, 20), collapse = ", ")
      )
    )
  }
  
  data_cols <- setdiff(names(data_DEG), "Group")
  
  p_values <- data.frame(
    Column = data_cols,
    P_value = NA_real_,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(data_cols)) {
    col <- data_cols[i]
    
    p_values$P_value[i] <- calculate_group_pvalue(
      value_vector = data_DEG[[col]],
      group_vector = data_DEG$Group,
      preferred_method = test_method
    )
  }
  
  data_DEG$Group <- NULL
  data_DEG <- rbind(data_DEG, p_values$P_value)
  rownames(data_DEG)[nrow(data_DEG)] <- "P_value"
  
  data_DEG <- t(data_DEG)
  data_DEG <- as.data.frame(data_DEG)
  data_DEG <- data_DEG[order(data_DEG$P_value, na.last = TRUE), ]
  
  write.csv(
    data_DEG,
    file = output_file,
    row.names = TRUE
  )
  
  return(data_DEG)
}

# 28. 通用函数：绘制分面箱线图

plot_score_facet_boxplot <- function(order_file,
                                     output_pdf,
                                     group_object,
                                     y_label = "Score") {
  
  score_data <- read.csv(
    order_file,
    header = TRUE,
    row.names = 1,
    check.names = FALSE
  )
  
  if ("P_value" %in% colnames(score_data)) {
    score_data$P_value <- NULL
  }
  
  if ("P_value" %in% rownames(score_data)) {
    score_data <- score_data[
      rownames(score_data) != "P_value",
      ,
      drop = FALSE
    ]
  }
  
  data_plot <- as.data.frame(t(score_data))
  
  data_plot$Group <- group_object[rownames(data_plot), "Group"]
  data_plot$Group <- factor(data_plot$Group, levels = group_names)
  
  if (anyNA(data_plot$Group)) {
    bad_id <- rownames(data_plot)[is.na(data_plot$Group)]
    stop(
      paste0(
        output_pdf,
        " 作图数据中有样本没有匹配到分组：",
        paste(head(bad_id, 20), collapse = ", ")
      )
    )
  }
  
  data_long <- reshape2::melt(
    data_plot,
    id.vars = "Group",
    variable.name = "Metric",
    value.name = "score"
  )
  
  p <- ggplot(
    data = data_long,
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
  
  p <- add_compare_layers(
    p = p,
    group_vector = data_long$Group,
    compare_mode = compare_mode,
    group_names = group_names
  )
  
  p <- p +
    facet_wrap(
      ~ Metric,
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
    ylab(y_label)
  
  ggsave(
    filename = output_pdf,
    plot = p,
    width = width.fig1,
    height = height.fig1
  )
  
  return(data_plot)
}

# 29. 对 ssGSEA / GSVA 结果按整体 P 值排序

ssgsea_order <- sort_score_by_pvalue(
  input_file = "ssgsea/ssgsea_score.csv",
  output_file = "ssgsea/ssgsea_order.csv",
  group_object = new_db,
  test_method = overall_test_method,
  result_name = "ssGSEA"
)

gsva_order <- sort_score_by_pvalue(
  input_file = "ssgsea/gsva_score.csv",
  output_file = "ssgsea/gsva_order.csv",
  group_object = new_db,
  test_method = overall_test_method,
  result_name = "GSVA"
)

# 30. 绘制 ssGSEA / GSVA 分面箱线图

ssgsea_plot_data <- plot_score_facet_boxplot(
  order_file = "ssgsea/ssgsea_order.csv",
  output_pdf = "ssgsea/1.ssgsea_boxplot_facet.pdf",
  group_object = new_db,
  y_label = "ssGSEA score"
)

gsva_plot_data <- plot_score_facet_boxplot(
  order_file = "ssgsea/gsva_order.csv",
  output_pdf = "ssgsea/1.gsva_boxplot_facet.pdf",
  group_object = new_db,
  y_label = "GSVA score"
)

# 31. 绘制 ssGSEA 单样本评分组成堆叠图

im_timer1 <- ssgsea_plot_data[
  ,
  setdiff(colnames(ssgsea_plot_data), "Group"),
  drop = FALSE
]

total_scores <- rowSums(im_timer1)

if (any(total_scores == 0, na.rm = TRUE)) {
  stop("存在总评分为 0 的样本，无法计算 ssGSEA 评分组成堆叠图。")
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
  file = "ssgsea/2.ssgsea_stacke.pdf",
  width = width.fig2,
  height = height.fig2
)

print(P1)
dev.off()

write.csv(
  proportion_scores,
  file = "ssgsea/proportion_scores.csv",
  row.names = TRUE
)

# 32. 目标基因和 ssGSEA 分数相关性分析

gene <- gene.duo
gene_keep <- gene[gene %in% rownames(expr)]

if (length(gene_keep) == 0) {
  stop("gene.duo 中没有任何基因存在于表达矩阵行名中。")
}

if (length(gene_keep) < length(gene)) {
  warning("gene.duo 中部分基因不存在于表达矩阵行名中，已自动跳过。")
}

gene_exp <- expr[gene_keep, , drop = FALSE]
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
  stop("表达矩阵样本名和 ssGSEA 分数样本名没有交集，无法做相关性分析。")
}

gene_exp_subset <- gene_exp[common_rows, , drop = FALSE]
im_timer1_subset <- im_timer1[common_rows, , drop = FALSE]

correlation_matrix <- cor(
  gene_exp_subset,
  im_timer1_subset,
  method = method.fig3
)

na_cols <- apply(is.na(correlation_matrix), 2, any)
correlation_matrix <- correlation_matrix[, !na_cols, drop = FALSE]

if (ncol(correlation_matrix) == 0) {
  stop("相关性矩阵全部为 NA，无法绘制热图。")
}

write.csv(
  correlation_matrix,
  file = "ssgsea/相关性数值.csv",
  row.names = TRUE
)

datacor <- t(correlation_matrix)

pdf(
  file = "ssgsea/3.多个基因免疫相关性图.pdf",
  width = width.fig3,
  height = height.fig3
)

pheatmap::pheatmap(
  datacor,
  scale = "none",
  cluster_cols = cluster1,
  cluster_rows = cluster1,
  color = colorRampPalette(color.fig3)(50),
  legend = TRUE,
  show_colnames = TRUE,
  show_rownames = TRUE,
  fontsize = fontsize,
  cellheight = cellheight,
  cellwidth = cellwidth
)

dev.off()

# 33. 绘制每个目标基因的相关性气泡图

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
    file = paste0("ssgsea/", siggene, "_immune_correlation.pdf"),
    width = width.fig4,
    height = height.fig4
  )
  
  print(p3)
  dev.off()
}

# 34. 读取小提琴图输入数据

violin_data <- read.csv(
  violin_score_file_path,
  header = TRUE,
  row.names = 1,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

violin_data <- t(violin_data)

if ("P_value" %in% rownames(violin_data)) {
  violin_data <- violin_data[
    rownames(violin_data) != "P_value",
    ,
    drop = FALSE
  ]
}

violin_data <- as.data.frame(violin_data)

# 35. 读取小提琴图分组文件

new_dbssgsea <- read.csv(
  violin_group_file_path,
  header = TRUE,
  row.names = 1,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

violin_data$group <- new_dbssgsea[rownames(violin_data), "Group"]
violin_data$group <- factor(violin_data$group, levels = group_names)

if (anyNA(violin_data$group)) {
  bad_id <- rownames(violin_data)[is.na(violin_data$group)]
  stop(
    paste0(
      "小提琴图数据中有样本没有匹配到分组：",
      paste(head(bad_id, 20), collapse = ", ")
    )
  )
}

# 36. 转换小提琴图长数据

cols_to_include <- setdiff(colnames(violin_data), "group")

data_long <- tidyr::pivot_longer(
  violin_data,
  cols = all_of(cols_to_include),
  names_to = "Metric",
  values_to = "Value"
)

# 37. 计算每个通路 / 基因集的整体 P 值

pvalues <- data_long %>%
  group_by(Metric) %>%
  summarise(
    p = calculate_group_pvalue(
      value_vector = Value,
      group_vector = group,
      preferred_method = overall_test_method
    ),
    .groups = "drop"
  )

# 38. 生成 facet 标签

if (p_display == "numeric") {
  pvalues <- pvalues %>%
    mutate(label = paste0("p = ", signif(p, digits = 3)))
} else if (p_display == "stars") {
  pvalues <- pvalues %>%
    mutate(
      label = case_when(
        is.na(p) ~ "NA",
        p < 0.001 ~ "***",
        p < 0.01 ~ "**",
        p < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    )
} else {
  stop("p_display 只能是 numeric 或 stars。")
}

new_labels <- setNames(
  paste0(pvalues$Metric, " (", pvalues$label, ")"),
  pvalues$Metric
)

# 39. 绘制小提琴图

p1 <- ggplot(
  data_long,
  aes(x = group, y = Value, fill = group)
) +
  geom_violin(
    trim = TRUE,
    alpha = violin_alpha,
    color = NA
  ) +
  {
    if (plot_points) {
      geom_jitter(
        shape = 16,
        position = position_jitter(0.2),
        aes(color = group),
        alpha = scatter_alpha,
        size = 2
      )
    } else {
      NULL
    }
  } +
  geom_boxplot(
    width = 0.5,
    color = "black",
    fill = NA,
    outlier.shape = NA
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.30))
  ) +
  scale_fill_manual(values = violin_colors_use) +
  scale_color_manual(values = violin_colors_use)

p1 <- add_compare_layers(
  p = p1,
  group_vector = data_long$group,
  compare_mode = compare_mode,
  group_names = group_names
)

p1 <- p1 +
  facet_wrap(
    ~ Metric,
    ncol = violin_facet_ncol,
    scales = "free_y",
    labeller = as_labeller(new_labels)
  ) +
  labs(title = " ", x = NULL, y = " ") +
  theme_minimal(base_size = 15) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.5
    ),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(
      fill = "grey",
      color = "black",
      linewidth = 0.5
    ),
    strip.text = element_text(size = 10),
    axis.ticks.y = element_line(color = "black")
  )

ggsave(
  filename = "score_violin_plot1.pdf",
  plot = p1,
  width = violin_pdf_width,
  height = violin_pdf_height,
  dpi = 300
)

# 40. 输出完成提示

cat("ssGSEA 和 GSVA 分析完成。\n")
cat("预处理后的表达矩阵对象：exprim\n")
cat("分组信息对象：new_db\n")
cat("ssGSEA 排序结果：ssgsea/ssgsea_order.csv\n")
cat("GSVA 排序结果：ssgsea/gsva_order.csv\n")
cat("ssGSEA 分面箱线图：ssgsea/1.ssgsea_boxplot_facet.pdf\n")
cat("GSVA 分面箱线图：ssgsea/1.gsva_boxplot_facet.pdf\n")
cat("小提琴图结果：score_violin_plot1.pdf\n")