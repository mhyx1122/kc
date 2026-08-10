# 1. 加载必要 R 包

library(data.table)
library(ggplot2)
library(ggpubr)


# 2. 数据读取参数

input_file <- "your_data.csv"

clinical_col <- "stage"

gene_names <- c("FUT4", "PIWIL4", "CST1")

diff_method <- "wilcox.test"
# 可选：
# diff_method <- "t.test"

stage_colors <- c(
  "#E69F00",
  "#56B4E9",
  "#009E73",
  "#F0E442"
)

plot_width <- 8
plot_height <- 6

output_dir_stage <- "临床相关的差异分析"


# 3. 读取数据

data11221 <- fread(
  input_file,
  header = TRUE,
  sep = ",",
  check.names = FALSE,
  data.table = FALSE
)

rownames(data11221) <- data11221[, 1]

data11221 <- data11221[, -1, drop = FALSE]


# 4. 设置分析数据和输出目录

gene_of_interest <- gene_names

output_dir <- output_dir_stage

data <- data11221

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


# 5. 清洗临床分期数据

data[[clinical_col]] <- ifelse(
  data[[clinical_col]] %in% c("Stage 0", "Stage X"),
  NA,
  data[[clinical_col]]
)

data <- data[!is.na(data[[clinical_col]]), , drop = FALSE]


# 6. 获取临床分组顺序

clinical_levels <- levels(factor(data[[clinical_col]]))


# 7. 循环绘制每个基因的临床相关箱线图

for (gene in gene_of_interest) {
  
  tryCatch({
    
    p <- ggplot(
      data,
      aes(
        x = .data[[clinical_col]],
        y = .data[[gene]],
        fill = .data[[clinical_col]]
      )
    ) +
      geom_boxplot(
        outlier.shape = NA,
        width = 0.6
      ) +
      geom_jitter(
        shape = 16,
        position = position_jitter(0.2),
        size = 1.5,
        alpha = 0.6
      ) +
      scale_fill_manual(
        values = stage_colors
      ) +
      labs(
        title = paste(gene, "Expression across", clinical_col),
        x = clinical_col,
        y = "Expression Level",
        fill = clinical_col
      ) +
      theme_minimal(
        base_size = 12
      ) +
      theme(
        legend.position = "top",
        plot.title = element_text(
          hjust = 0.5,
          face = "bold",
          size = 12
        ),
        axis.text.x = element_text(
          angle = 45,
          hjust = 1
        ),
        panel.grid = element_blank(),
        axis.line = element_line(
          color = "black",
          linewidth = 0.6
        )
      )
    
    my_comparisons <- combn(
      clinical_levels,
      2,
      simplify = FALSE
    )
    
    p <- p +
      stat_compare_means(
        method = diff_method,
        label = "p.format",
        comparisons = my_comparisons
      )
    
    output_file <- file.path(
      output_dir,
      paste0(gene, "_boxplot_expression_", clinical_col, ".pdf")
    )
    
    ggsave(
      output_file,
      plot = p,
      width = plot_width,
      height = plot_height
    )
    
    message(
      paste(
        "基因",
        gene,
        "在",
        clinical_col,
        "中的表达箱线图已保存:",
        output_file
      )
    )
    
  }, error = function(e) {
    
    message(
      paste(
        "处理基因",
        gene,
        "时发生错误，已跳过该基因。错误信息:",
        e$message
      )
    )
    
  })
}