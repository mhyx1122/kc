# 1. 加载必要 R 包

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(ggpubr)
  library(singscore)
  library(AUCell)
  library(GSEABase)
  library(GSVA)
  library(UCell)
})


# 2. 检查全局环境中的 Seurat 对象

if (!exists("seurat", envir = .GlobalEnv)) {
  stop("全局环境中没有 seurat 对象，请先加载 seurat 对象。")
}

seurat_obj <- get("seurat", envir = .GlobalEnv)

if (!inherits(seurat_obj, "Seurat")) {
  stop("全局环境中的 seurat 不是 Seurat 对象。")
}

if (!"cellType" %in% colnames(seurat_obj@meta.data)) {
  stop("seurat@meta.data 中没有 cellType 列，请先添加细胞类型注释。")
}


# 3. 模块一：设置比较参数

cell_types <- unique(seurat_obj@meta.data$cellType)
cell_types <- cell_types[!is.na(cell_types)]

comparisons_list <- list(
  c("Fibroblasts", "Monocytes")
)

# 如果需要多个比较组，可以按下面格式继续添加
# comparisons_list <- list(
#   c("CellTypeA", "CellTypeB"),
#   c("CellTypeC", "CellTypeD")
# )


# 4. 模块二：读取基因打分集合文件

gene_data_file <- "19PCDgenes.csv"

gene_set_table <- read.csv(
  gene_data_file,
  header = TRUE,
  row.names = NULL,
  check.names = FALSE
)

gene_set_list <- list()

for (gene_set_name in colnames(gene_set_table)) {
  
  gene_vector <- trimws(gene_set_table[[gene_set_name]])
  gene_vector <- gene_vector[gene_vector != ""]
  gene_vector <- gene_vector[!is.na(gene_vector)]
  
  gene_set_list[[gene_set_name]] <- gene_vector
}


# 5. 模块三：设置小提琴图和 UMAP 图参数

violin_width <- 8
violin_height <- 6

umap_width <- 5.5
umap_height <- 5

auc_violin_width <- 8
auc_violin_height <- 6

auc_pdf_width <- 6
auc_pdf_height <- 5


# 6. 模块四：设置 UCell 参数

ucell_max_rank <- 5000


# 7. 模块五：设置颜色参数

feature_plot_colors <- c("grey", "blue")

celltype_colors <- c(
  "#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF",
  "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#7E6148FF",
  "#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF",
  "#7AA6DCFF", "#003C67FF", "#8F7700FF", "#3B3B3BFF",
  "#A73030FF", "#374E55FF", "#DF8F44FF", "#00A1D5FF",
  "#B24745FF", "#79AF97FF", "#6A6599FF", "#80796BFF",
  "#BC3C29FF", "#0072B5FF", "#E18727FF", "#20854EFF",
  "#7876B1FF", "#6F99ADFF", "#FFDC91FF", "#EE4C97FF"
)


# 8. 模块六：创建输出目录

output_folder <- "10.1单细胞评分"

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}


# 9. 模块七：提取 Seurat 表达矩阵

analysis_seurat <- seurat_obj

expression_matrix <- tryCatch(
  {
    GetAssayData(
      object = analysis_seurat,
      assay = "RNA",
      layer = "data"
    )
  },
  error = function(e) {
    GetAssayData(
      object = analysis_seurat,
      assay = "RNA",
      slot = "data"
    )
  }
)

try({
  analysis_seurat@assays$RNA@layers$data@Dimnames <- expression_matrix@Dimnames
}, silent = TRUE)


# 10. 模块八：循环分析每一个基因集

for (gene_set_name in names(gene_set_list)) {
  
  cat(paste0("开始分析基因集：", gene_set_name, "\n"))
  
  gene_set_output_folder <- file.path(output_folder, gene_set_name)
  
  if (!dir.exists(gene_set_output_folder)) {
    dir.create(gene_set_output_folder, recursive = TRUE)
  }
  
  current_gene_set <- gene_set_list[[gene_set_name]]
  
  score_result_table <- data.frame(
    row.names = rownames(analysis_seurat@meta.data)
  )
  
  
  # 11. 模块九：singscore 单细胞评分
  
  tryCatch({
    
    expression_df <- as.data.frame(expression_matrix)
    
    ranked_expression <- singscore::rankGenes(
      as.data.frame(expression_df)
    )
    
    singscore_result <- simpleScore(
      ranked_expression,
      upSet = current_gene_set,
      knownDirection = FALSE
    )
    
    singscore_column <- paste0("singscore_", gene_set_name)
    
    analysis_seurat@meta.data[[singscore_column]] <- singscore_result$TotalScore
    
    score_result_table[[singscore_column]] <- analysis_seurat@meta.data[
      rownames(score_result_table),
      singscore_column
    ]
    
    singscore_violin_plot <- ggviolin(
      analysis_seurat@meta.data,
      x = "cellType",
      y = singscore_column,
      color = "cellType",
      add = "mean_sd",
      fill = "cellType",
      add.params = list(color = "black")
    ) +
      stat_compare_means(
        comparisons = comparisons_list,
        label = "p.signif"
      ) +
      scale_color_manual(values = celltype_colors) +
      scale_fill_manual(values = celltype_colors) +
      theme(
        axis.text.x.bottom = element_text(
          angle = 90,
          vjust = 0.5,
          hjust = 1
        )
      ) +
      NoLegend() +
      labs(x = "")
    
    ggsave(
      file.path(gene_set_output_folder, "1.singscore单细胞评分——小提琴图.pdf"),
      plot = singscore_violin_plot,
      width = violin_width,
      height = violin_height
    )
    
    singscore_umap_plot <- FeaturePlot(
      analysis_seurat,
      singscore_column,
      cols = feature_plot_colors
    )
    
    ggsave(
      file.path(gene_set_output_folder, "1.singscore单细胞评分——umap分布图.pdf"),
      plot = singscore_umap_plot,
      width = umap_width,
      height = umap_height
    )
    
    cat(paste0("singscore method succeeded for ", gene_set_name, "\n"))
    
    rm(expression_df, ranked_expression, singscore_result, singscore_violin_plot, singscore_umap_plot)
    gc()
    
  }, error = function(e) {
    cat(paste0("singscore method failed for ", gene_set_name, " with error: ", e$message, "\n"))
  })
  
  
  # 12. 模块十：AddModuleScore 单细胞评分
  
  tryCatch({
    
    add_module_features <- list(current_gene_set)
    
    add_module_prefix <- paste(
      "AddModuleScore",
      gene_set_name,
      sep = "_"
    )
    
    analysis_seurat <- AddModuleScore(
      analysis_seurat,
      features = add_module_features,
      name = add_module_prefix
    )
    
    add_module_column <- paste0(add_module_prefix, "1")
    
    score_result_table[[add_module_column]] <- analysis_seurat@meta.data[
      rownames(score_result_table),
      add_module_column
    ]
    
    add_module_violin_plot <- ggviolin(
      analysis_seurat@meta.data,
      x = "cellType",
      y = add_module_column,
      color = "cellType",
      add = "mean_sd",
      fill = "cellType",
      add.params = list(color = "black")
    ) +
      stat_compare_means(
        comparisons = comparisons_list,
        label = "p.signif"
      ) +
      scale_color_manual(values = celltype_colors) +
      scale_fill_manual(values = celltype_colors) +
      theme(
        axis.text.x.bottom = element_text(
          angle = 90,
          vjust = 0.5,
          hjust = 1
        )
      ) +
      NoLegend() +
      labs(x = "")
    
    ggsave(
      file.path(gene_set_output_folder, "2.AddModuleScor单细胞评分——小提琴图.pdf"),
      plot = add_module_violin_plot,
      width = violin_width,
      height = violin_height
    )
    
    add_module_umap_plot <- FeaturePlot(
      analysis_seurat,
      add_module_column,
      cols = feature_plot_colors
    )
    
    ggsave(
      file.path(gene_set_output_folder, "2.AddModuleScor单细胞评分——umap分布图.pdf"),
      plot = add_module_umap_plot,
      width = umap_width,
      height = umap_height
    )
    
    cat(paste0("AddModuleScore method succeeded for ", gene_set_name, "\n"))
    
    rm(add_module_features, add_module_violin_plot, add_module_umap_plot)
    gc()
    
  }, error = function(e) {
    cat(paste0("AddModuleScore method failed for ", gene_set_name, " with error: ", e$message, "\n"))
  })
  
  
  # 13. 模块十一：AUCell 单细胞评分
  
  tryCatch({
    
    aucell_rankings <- AUCell_buildRankings(
      analysis_seurat@assays$RNA@layers$data,
      plotStats = FALSE
    )
    
    aucell_gene_set <- c(
      GeneSet(
        sample(current_gene_set),
        setName = "AUCell"
      )
    )
    
    aucell_gene_collection <- GeneSetCollection(aucell_gene_set)
    
    aucell_result <- AUCell_calcAUC(
      aucell_gene_collection,
      aucell_rankings
    )
    
    pdf(
      file.path(gene_set_output_folder, "3.AUCell评分.pdf"),
      width = auc_pdf_width,
      height = auc_pdf_height
    )
    
    aucell_threshold_result <- AUCell_exploreThresholds(
      aucell_result,
      plotHist = TRUE,
      assign = TRUE
    )
    
    dev.off()
    
    aucell_score <- as.numeric(
      getAUC(aucell_result)["AUCell", ]
    )
    
    aucell_column <- paste0("AUC_", gene_set_name)
    
    analysis_seurat@meta.data[[aucell_column]] <- aucell_score
    
    score_result_table[[aucell_column]] <- analysis_seurat@meta.data[
      rownames(score_result_table),
      aucell_column
    ]
    
    aucell_umap_plot <- FeaturePlot(
      analysis_seurat,
      aucell_column,
      cols = feature_plot_colors
    )
    
    ggsave(
      file.path(gene_set_output_folder, "3.AUC单细胞评分——umap分布图.pdf"),
      plot = aucell_umap_plot,
      width = umap_width,
      height = umap_height
    )
    
    aucell_violin_plot <- ggviolin(
      analysis_seurat@meta.data,
      x = "cellType",
      y = aucell_column,
      color = "cellType",
      add = "mean_sd",
      fill = "cellType",
      add.params = list(color = "black")
    ) +
      stat_compare_means(
        comparisons = comparisons_list,
        label = "p.signif"
      ) +
      scale_color_manual(values = celltype_colors) +
      scale_fill_manual(values = celltype_colors) +
      theme(
        axis.text.x.bottom = element_text(
          angle = 90,
          vjust = 0.5,
          hjust = 1
        )
      ) +
      NoLegend() +
      labs(x = "")
    
    ggsave(
      file.path(gene_set_output_folder, "3.AUCell评分——小提琴图.pdf"),
      plot = aucell_violin_plot,
      width = auc_violin_width,
      height = auc_violin_height
    )
    
    cat(paste0("AUCell method succeeded for ", gene_set_name, "\n"))
    
    rm(
      aucell_rankings,
      aucell_gene_set,
      aucell_gene_collection,
      aucell_result,
      aucell_threshold_result,
      aucell_score,
      aucell_umap_plot,
      aucell_violin_plot
    )
    gc()
    
  }, error = function(e) {
    cat(paste0("AUCell method failed for ", gene_set_name, " with error: ", e$message, "\n"))
  })
  
  
  # 14. 模块十二：ssGSEA 单细胞评分
  
  tryCatch({
    
    ssgsea_expression_matrix <- as.matrix(
      analysis_seurat@assays$RNA@layers$data
    )
    
    ssgsea_gene_set <- as.data.frame(current_gene_set)
    
    colnames(ssgsea_gene_set) <- "ssGSEA"
    
    ssgsea_param <- ssgseaParam(
      exprData = ssgsea_expression_matrix,
      geneSets = ssgsea_gene_set
    )
    
    ssgsea_result <- gsva(
      ssgsea_param,
      verbose = TRUE
    )
    
    ssgsea_plot_data <- cbind(
      analysis_seurat@meta.data,
      t(ssgsea_result)[rownames(analysis_seurat@meta.data), ]
    )
    
    colnames(ssgsea_plot_data)[ncol(ssgsea_plot_data)] <- "ssGSEA"
    
    ssgsea_score <- ssgsea_plot_data$ssGSEA[
      match(
        rownames(analysis_seurat@meta.data),
        rownames(ssgsea_plot_data)
      )
    ]
    
    ssgsea_column <- paste0("ssGSEA_", gene_set_name)
    
    analysis_seurat@meta.data[[ssgsea_column]] <- ssgsea_score
    
    score_result_table[[ssgsea_column]] <- analysis_seurat@meta.data[
      rownames(score_result_table),
      ssgsea_column
    ]
    
    ssgsea_umap_plot <- FeaturePlot(
      analysis_seurat,
      ssgsea_column,
      cols = feature_plot_colors
    )
    
    ggsave(
      file.path(gene_set_output_folder, "4.ssGSEA单细胞评分——umap分布图.pdf"),
      plot = ssgsea_umap_plot,
      width = umap_width,
      height = umap_height
    )
    
    ssgsea_violin_plot <- ggviolin(
      analysis_seurat@meta.data,
      x = "cellType",
      y = ssgsea_column,
      color = "cellType",
      add = "mean_sd",
      fill = "cellType",
      add.params = list(color = "black")
    ) +
      stat_compare_means(
        comparisons = comparisons_list,
        label = "p.signif"
      ) +
      scale_color_manual(values = celltype_colors) +
      scale_fill_manual(values = celltype_colors) +
      theme(
        axis.text.x.bottom = element_text(
          angle = 90,
          vjust = 0.5,
          hjust = 1
        )
      ) +
      NoLegend() +
      labs(x = "")
    
    ggsave(
      file.path(gene_set_output_folder, "4.ssGSEA单细胞评分——小提琴图.pdf"),
      plot = ssgsea_violin_plot,
      width = violin_width,
      height = violin_height
    )
    
    cat(paste0("ssGSEA method succeeded for ", gene_set_name, "\n"))
    
    rm(
      ssgsea_expression_matrix,
      ssgsea_gene_set,
      ssgsea_param,
      ssgsea_result,
      ssgsea_plot_data,
      ssgsea_score,
      ssgsea_umap_plot,
      ssgsea_violin_plot
    )
    gc()
    
  }, error = function(e) {
    cat(paste0("ssGSEA method failed for ", gene_set_name, " with error: ", e$message, "\n"))
  })
  
  
  # 15. 模块十三：UCell 单细胞评分
  
  tryCatch({
    
    ucell_features <- list(current_gene_set)
    
    analysis_seurat <- AddModuleScore_UCell(
      analysis_seurat,
      features = ucell_features,
      name = "UCell",
      maxRank = ucell_max_rank
    )
    
    colnames(analysis_seurat@meta.data)[
      colnames(analysis_seurat@meta.data) == "signature_1UCell"
    ] <- paste0("UCell_", gene_set_name)
    
    ucell_column <- paste0("UCell_", gene_set_name)
    
    score_result_table[[ucell_column]] <- analysis_seurat@meta.data[
      rownames(score_result_table),
      ucell_column
    ]
    
    ucell_violin_plot <- ggviolin(
      analysis_seurat@meta.data,
      x = "cellType",
      y = ucell_column,
      color = "cellType",
      add = "mean_sd",
      fill = "cellType",
      add.params = list(color = "black")
    ) +
      stat_compare_means(
        comparisons = comparisons_list,
        label = "p.signif"
      ) +
      scale_color_manual(values = celltype_colors) +
      scale_fill_manual(values = celltype_colors) +
      theme(
        axis.text.x.bottom = element_text(
          angle = 90,
          vjust = 0.5,
          hjust = 1
        )
      ) +
      NoLegend() +
      labs(x = "")
    
    ggsave(
      file.path(gene_set_output_folder, "5.Ucell单细胞评分——小提琴图.pdf"),
      plot = ucell_violin_plot,
      width = violin_width,
      height = violin_height
    )
    
    ucell_umap_plot <- FeaturePlot(
      analysis_seurat,
      ucell_column,
      cols = feature_plot_colors
    )
    
    ggsave(
      file.path(gene_set_output_folder, "5.Ucell单细胞评分——umap分布图.pdf"),
      plot = ucell_umap_plot,
      width = umap_width,
      height = umap_height
    )
    
    score_result_table$cellType <- analysis_seurat@meta.data[
      rownames(score_result_table),
      "cellType"
    ]
    
    write.csv(
      score_result_table,
      file = file.path(gene_set_output_folder, "五种评分方法数据_score.csv"),
      row.names = TRUE
    )
    
    cat(paste0("UCell method succeeded for ", gene_set_name, "\n"))
    
    rm(ucell_features, ucell_violin_plot, ucell_umap_plot)
    gc()
    
  }, error = function(e) {
    cat(paste0("UCell method failed for ", gene_set_name, " with error: ", e$message, "\n"))
  })
  
  cat(paste0("完成基因集：", gene_set_name, "\n"))
}


# 16. 模块十四：保存更新后的 Seurat 对象到全局环境

seurat <- analysis_seurat

cat("所有基因集的五种单细胞评分分析已完成。\n")  