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

# 1. 读取评分基因集和Spatial_Data对象

# 1.1 参数设置

gene_data_file <- "基因打分集合.csv"
output_folder <- "10.1单细胞评分"

# 根据Spatial_Data$cellType中的真实名称修改比较组合
comparisons_list <- list(
  c("TCells", "B.Cells"),
  c("TCells", "Endothelial.Cells")
)

# 小提琴图保存尺寸
violin_width <- 8
violin_height <- 6

# UMAP图保存尺寸
umap_width <- 5.5
umap_height <- 5

# AUCell小提琴图保存尺寸
auc_violin_width <- 8
auc_violin_height <- 6

# AUCell阈值图保存尺寸
auc_pdf_width <- 6
auc_pdf_height <- 5

# UCell最大排名数量
max_rank <- 5000

# FeaturePlot渐变颜色
feature_plot_colors <- c("grey", "blue")

# 细胞类型颜色
color_palette <- c(
  "#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF",
  "#91D1C2FF", "#7E6148FF", "#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF",
  "#7AA6DCFF", "#003C67FF", "#8F7700FF", "#3B3B3BFF", "#A73030FF", "#374E55FF",
  "#DF8F44FF", "#00A1D5FF", "#B24745FF", "#79AF97FF", "#6A6599FF", "#80796BFF",
  "#BC3C29FF", "#0072B5FF", "#E18727FF", "#20854EFF", "#7876B1FF", "#6F99ADFF",
  "#FFDC91FF", "#EE4C97FF"
)

# 1.2 检查并读取Spatial_Data对象

if (!exists("Spatial_Data", envir = .GlobalEnv)) {
  stop("全局环境中不存在 Spatial_Data 对象。")
}

analysis_seurat <- get("Spatial_Data", envir = .GlobalEnv)

if (!inherits(analysis_seurat, "Seurat")) {
  stop("Spatial_Data 不是有效的Seurat对象。")
}

if (!"cellType" %in% colnames(analysis_seurat@meta.data)) {
  stop("Spatial_Data@meta.data 中不存在 cellType 列。")
}

# 1.3 检查评分基因集文件

if (!file.exists(gene_data_file)) {
  stop(paste0("未找到评分基因集文件：", gene_data_file))
}

gene_data <- read.csv(
  gene_data_file,
  header = TRUE,
  row.names = NULL,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (ncol(gene_data) == 0) {
  stop("评分基因集文件中没有可用的列。")
}

# 1.4 将每一列整理为一个评分基因集

gene_sets <- list()

for (gene_set_name in colnames(gene_data)) {
  current_genes <- as.character(gene_data[[gene_set_name]])
  current_genes <- current_genes[!is.na(current_genes)]
  current_genes <- trimws(current_genes)
  current_genes <- current_genes[current_genes != ""]
  current_genes <- unique(current_genes)
  
  gene_sets[[gene_set_name]] <- current_genes
}

empty_gene_sets <- names(gene_sets)[lengths(gene_sets) == 0]

if (length(empty_gene_sets) > 0) {
  stop(
    paste0(
      "以下评分基因集为空：",
      paste(empty_gene_sets, collapse = ", ")
    )
  )
}

# 1.5 检查细胞比较组合

cell_types <- unique(as.character(analysis_seurat@meta.data$cellType))

if (length(comparisons_list) == 0) {
  stop("comparisons_list不能为空，请至少设置一组比较对象。")
}

for (comparison_index in seq_along(comparisons_list)) {
  current_comparison <- comparisons_list[[comparison_index]]
  
  if (length(current_comparison) != 2) {
    stop(
      paste0(
        "comparisons_list中的第",
        comparison_index,
        "组比较必须包含两个细胞类型。"
      )
    )
  }
  
  missing_cell_types <- setdiff(current_comparison, cell_types)
  
  if (length(missing_cell_types) > 0) {
    stop(
      paste0(
        "第",
        comparison_index,
        "组比较中的以下细胞类型不存在：",
        paste(missing_cell_types, collapse = ", "),
        "。当前可用细胞类型为：",
        paste(cell_types, collapse = ", ")
      )
    )
  }
}

# 1.6 检查颜色数量

if (length(color_palette) < length(cell_types)) {
  stop(
    paste0(
      "细胞类型颜色数量不足。当前有",
      length(cell_types),
      "种细胞类型，但只提供了",
      length(color_palette),
      "个颜色。"
    )
  )
}

# 1.7 检查SCT assay

if (!"SCT" %in% Assays(analysis_seurat)) {
  stop(
    paste0(
      "Spatial_Data中不存在SCT assay。",
      "原评分代码中的AUCell和ssGSEA依赖SCT assay的data层。"
    )
  )
}

# 1.8 创建主输出目录

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

# 1.9 提取默认assay的data表达矩阵

expression_data <- GetAssayData(
  object = analysis_seurat,
  layer = "data"
)

# 2. 依次分析每一个评分基因集

for (gene_set_name in names(gene_sets)) {
  
  current_gene_set <- gene_sets[[gene_set_name]]
  
  sub_output_folder <- file.path(
    output_folder,
    gene_set_name
  )
  
  if (!dir.exists(sub_output_folder)) {
    dir.create(sub_output_folder, recursive = TRUE)
  }
  
  # 提前创建本基因集的评分结果汇总表
  score_data <- data.frame(
    row.names = rownames(analysis_seurat@meta.data)
  )
  
  # AddModuleScore和UCell共同使用的基因集格式
  module_features <- list(current_gene_set)
  
  # 3. singscore评分
  
  tryCatch(
    {
      ranked_expression <- as.data.frame(expression_data)
      ranked_expression <- singscore::rankGenes(ranked_expression)
      
      singscore_result <- singscore::simpleScore(
        ranked_expression,
        upSet = current_gene_set,
        knownDirection = FALSE
      )
      
      singscore_column <- paste0(
        "singscore_",
        gene_set_name
      )
      
      analysis_seurat@meta.data[[singscore_column]] <-
        singscore_result$TotalScore
      
      score_data[[singscore_column]] <-
        analysis_seurat@meta.data[
          rownames(score_data),
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
        scale_color_manual(values = color_palette) +
        scale_fill_manual(values = color_palette) +
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
        filename = file.path(
          sub_output_folder,
          "1.singscore单细胞评分——小提琴图.pdf"
        ),
        plot = singscore_violin_plot,
        width = violin_width,
        height = violin_height
      )
      
      singscore_feature_plot <- FeaturePlot(
        analysis_seurat,
        features = singscore_column,
        cols = feature_plot_colors
      )
      
      ggsave(
        filename = file.path(
          sub_output_folder,
          "1.singscore单细胞评分——umap分布图.pdf"
        ),
        plot = singscore_feature_plot,
        width = umap_width,
        height = umap_height
      )
      
      message(
        "singscore method succeeded for ",
        gene_set_name
      )
      
      rm(
        ranked_expression,
        singscore_result,
        singscore_violin_plot,
        singscore_feature_plot
      )
      
      gc()
    },
    error = function(e) {
      message(
        "singscore method failed for ",
        gene_set_name,
        " with error: ",
        e$message
      )
    }
  )
  
  # 4. AddModuleScore评分
  
  tryCatch(
    {
      addmodule_base_name <- paste(
        "AddModuleScore",
        gene_set_name,
        sep = "_"
      )
      
      analysis_seurat <- AddModuleScore(
        analysis_seurat,
        features = module_features,
        name = addmodule_base_name
      )
      
      addmodule_column <- paste0(
        addmodule_base_name,
        "1"
      )
      
      score_data[[addmodule_column]] <-
        analysis_seurat@meta.data[
          rownames(score_data),
          addmodule_column
        ]
      
      addmodule_violin_plot <- ggviolin(
        analysis_seurat@meta.data,
        x = "cellType",
        y = addmodule_column,
        color = "cellType",
        add = "mean_sd",
        fill = "cellType",
        add.params = list(color = "black")
      ) +
        stat_compare_means(
          comparisons = comparisons_list,
          label = "p.signif"
        ) +
        scale_color_manual(values = color_palette) +
        scale_fill_manual(values = color_palette) +
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
        filename = file.path(
          sub_output_folder,
          "2.AddModuleScor单细胞评分——小提琴图.pdf"
        ),
        plot = addmodule_violin_plot,
        width = violin_width,
        height = violin_height
      )
      
      addmodule_feature_plot <- FeaturePlot(
        analysis_seurat,
        features = addmodule_column,
        cols = feature_plot_colors
      )
      
      ggsave(
        filename = file.path(
          sub_output_folder,
          "2.AddModuleScor单细胞评分——umap分布图.pdf"
        ),
        plot = addmodule_feature_plot,
        width = umap_width,
        height = umap_height
      )
      
      message(
        "AddModuleScore method succeeded for ",
        gene_set_name
      )
      
      rm(
        addmodule_violin_plot,
        addmodule_feature_plot
      )
      
      gc()
    },
    error = function(e) {
      message(
        "AddModuleScore method failed for ",
        gene_set_name,
        " with error: ",
        e$message
      )
    }
  )
  
  # 5. AUCell评分
  
  tryCatch(
    {
      sct_expression_auc <- GetAssayData(
        object = analysis_seurat,
        assay = "SCT",
        layer = "data"
      )
      
      auc_rankings <- AUCell_buildRankings(
        sct_expression_auc,
        plotStats = FALSE
      )
      
      auc_gene_set <- c(
        GeneSet(
          sample(current_gene_set),
          setName = "AUCell"
        )
      )
      
      auc_gene_collection <- GeneSetCollection(
        auc_gene_set
      )
      
      auc_result <- AUCell_calcAUC(
        auc_gene_collection,
        auc_rankings
      )
      
      pdf(
        file = file.path(
          sub_output_folder,
          "3.AUCell评分.pdf"
        ),
        width = auc_pdf_width,
        height = auc_pdf_height
      )
      
      auc_threshold_result <- AUCell_exploreThresholds(
        auc_result,
        plotHist = TRUE,
        assign = TRUE
      )
      
      dev.off()
      
      auc_values <- as.numeric(
        getAUC(auc_result)["AUCell", ]
      )
      
      auc_column <- paste0(
        "AUC_",
        gene_set_name
      )
      
      analysis_seurat@meta.data[[auc_column]] <- auc_values
      
      score_data[[auc_column]] <-
        analysis_seurat@meta.data[
          rownames(score_data),
          auc_column
        ]
      
      auc_feature_plot <- FeaturePlot(
        analysis_seurat,
        features = auc_column,
        cols = feature_plot_colors
      )
      
      ggsave(
        filename = file.path(
          sub_output_folder,
          "3.AUC单细胞评分——umap分布图.pdf"
        ),
        plot = auc_feature_plot,
        width = umap_width,
        height = umap_height
      )
      
      auc_violin_plot <- ggviolin(
        analysis_seurat@meta.data,
        x = "cellType",
        y = auc_column,
        color = "cellType",
        add = "mean_sd",
        fill = "cellType",
        add.params = list(color = "black")
      ) +
        stat_compare_means(
          comparisons = comparisons_list,
          label = "p.signif"
        ) +
        scale_color_manual(values = color_palette) +
        scale_fill_manual(values = color_palette) +
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
        filename = file.path(
          sub_output_folder,
          "3.AUCell评分——小提琴图.pdf"
        ),
        plot = auc_violin_plot,
        width = auc_violin_width,
        height = auc_violin_height
      )
      
      message(
        "AUCell method succeeded for ",
        gene_set_name
      )
      
      rm(
        sct_expression_auc,
        auc_rankings,
        auc_gene_set,
        auc_gene_collection,
        auc_result,
        auc_threshold_result,
        auc_values,
        auc_feature_plot,
        auc_violin_plot
      )
      
      gc()
    },
    error = function(e) {
      message(
        "AUCell method failed for ",
        gene_set_name,
        " with error: ",
        e$message
      )
    }
  )
  
  # 6. ssGSEA评分
  
  tryCatch(
    {
      sct_expression_ssgsea <- as.matrix(
        GetAssayData(
          object = analysis_seurat,
          assay = "SCT",
          layer = "data"
        )
      )
      
      ssgsea_gene_set <- as.data.frame(
        current_gene_set
      )
      
      colnames(ssgsea_gene_set) <- "ssGSEA"
      
      ssgsea_parameter <- ssgseaParam(
        exprData = sct_expression_ssgsea,
        geneSets = ssgsea_gene_set
      )
      
      ssgsea_result <- gsva(
        ssgsea_parameter,
        verbose = TRUE
      )
      
      ssgsea_plot_data <- cbind(
        analysis_seurat@meta.data,
        t(ssgsea_result)[
          rownames(analysis_seurat@meta.data),
          ,
          drop = FALSE
        ]
      )
      
      colnames(ssgsea_plot_data)[ncol(ssgsea_plot_data)] <- "ssGSEA"
      
      ssgsea_vector <- ssgsea_plot_data$ssGSEA[
        match(
          rownames(analysis_seurat@meta.data),
          rownames(ssgsea_plot_data)
        )
      ]
      
      ssgsea_column <- paste0(
        "ssGSEA_",
        gene_set_name
      )
      
      analysis_seurat@meta.data[[ssgsea_column]] <-
        ssgsea_vector
      
      score_data[[ssgsea_column]] <-
        analysis_seurat@meta.data[
          rownames(score_data),
          ssgsea_column
        ]
      
      ssgsea_feature_plot <- FeaturePlot(
        analysis_seurat,
        features = ssgsea_column,
        cols = feature_plot_colors
      )
      
      ggsave(
        filename = file.path(
          sub_output_folder,
          "4.ssGSEA单细胞评分——umap分布图.pdf"
        ),
        plot = ssgsea_feature_plot,
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
        scale_color_manual(values = color_palette) +
        scale_fill_manual(values = color_palette) +
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
        filename = file.path(
          sub_output_folder,
          "4.ssGSEA单细胞评分——小提琴图.pdf"
        ),
        plot = ssgsea_violin_plot,
        width = violin_width,
        height = violin_height
      )
      
      message(
        "ssGSEA method succeeded for ",
        gene_set_name
      )
      
      rm(
        sct_expression_ssgsea,
        ssgsea_gene_set,
        ssgsea_parameter,
        ssgsea_result,
        ssgsea_plot_data,
        ssgsea_vector,
        ssgsea_feature_plot,
        ssgsea_violin_plot
      )
      
      gc()
    },
    error = function(e) {
      message(
        "ssGSEA method failed for ",
        gene_set_name,
        " with error: ",
        e$message
      )
    }
  )
  
  # 7. UCell评分并导出五种评分结果
  
  tryCatch(
    {
      analysis_seurat <- AddModuleScore_UCell(
        analysis_seurat,
        features = module_features,
        name = "UCell",
        maxRank = max_rank
      )
      
      colnames(analysis_seurat@meta.data)[
        colnames(analysis_seurat@meta.data) == "signature_1UCell"
      ] <- paste0(
        "UCell_",
        gene_set_name
      )
      
      ucell_column <- paste0(
        "UCell_",
        gene_set_name
      )
      
      score_data[[ucell_column]] <-
        analysis_seurat@meta.data[
          rownames(score_data),
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
        scale_color_manual(values = color_palette) +
        scale_fill_manual(values = color_palette) +
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
        filename = file.path(
          sub_output_folder,
          "5.Ucell单细胞评分——小提琴图.pdf"
        ),
        plot = ucell_violin_plot,
        width = violin_width,
        height = violin_height
      )
      
      ucell_feature_plot <- FeaturePlot(
        analysis_seurat,
        features = ucell_column,
        cols = feature_plot_colors
      )
      
      ggsave(
        filename = file.path(
          sub_output_folder,
          "5.Ucell单细胞评分——umap分布图.pdf"
        ),
        plot = ucell_feature_plot,
        width = umap_width,
        height = umap_height
      )
      
      score_data$cellType <-
        analysis_seurat@meta.data[
          rownames(score_data),
          "cellType"
        ]
      
      write.csv(
        score_data,
        file = file.path(
          sub_output_folder,
          "五种评分方法数据_score.csv"
        ),
        row.names = TRUE
      )
      
      message(
        "UCell method succeeded for ",
        gene_set_name
      )
      
      rm(
        ucell_violin_plot,
        ucell_feature_plot
      )
      
      gc()
    },
    error = function(e) {
      message(
        "UCell method failed for ",
        gene_set_name,
        " with error: ",
        e$message
      )
    }
  )
}

# 8. 更新全局环境中的Spatial_Data对象

assign(
  "Spatial_Data",
  analysis_seurat,
  envir = .GlobalEnv
)