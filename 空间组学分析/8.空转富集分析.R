suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(clusterProfiler)
  library(enrichplot)
})

# 1. 准备各cluster的Top基因

# 1.1 参数设置

out_dir <- "6.KEGG-GO富集分析"

# up表示上调基因，down表示下调基因
gene_direction <- "up"

# 每个cluster选取的Top基因数量
top_gene_number <- 200

# Top基因表文件名
name_top_gene_csv <- "top200gene"

# OrgDb数据库包名，保持字符形式
orgdb_object <- "org.Hs.eg.db"

# 1.2 创建输出文件夹

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# 1.3 检查并读取all_markers对象

if (!exists("all_markers", envir = .GlobalEnv)) {
  stop("全局环境中没有 all_markers 对象，请先完成差异表达分析。")
}

all_markers <- get("all_markers", envir = .GlobalEnv)

required_cols <- c("gene", "cluster", "avg_log2FC")
missing_cols <- setdiff(required_cols, colnames(all_markers))

if (length(missing_cols) > 0) {
  stop(
    paste0(
      "all_markers 缺少必要列：",
      paste(missing_cols, collapse = ", ")
    )
  )
}

# 1.4 根据基因方向筛选差异基因

if (gene_direction == "up") {
  selected_markers <- all_markers %>%
    filter(avg_log2FC > 0)
} else if (gene_direction == "down") {
  selected_markers <- all_markers %>%
    filter(avg_log2FC < 0)
} else {
  stop("gene_direction 只能设置为 up 或 down。")
}

if (nrow(selected_markers) == 0) {
  stop("筛选后的 selected_markers 为空，请检查 all_markers 或基因方向设置。")
}

# 1.5 每个cluster选取Top基因

selected_markers <- selected_markers %>%
  group_by(cluster) %>%
  slice_max(
    order_by = abs(avg_log2FC),
    n = top_gene_number,
    with_ties = FALSE
  ) %>%
  ungroup()

# 1.6 按cluster整理Top基因表

top_genes_by_cluster <- unstack(
  selected_markers,
  gene ~ cluster
)

# 1.7 保存Top基因表

top_gene_file <- file.path(
  out_dir,
  paste0(name_top_gene_csv, ".csv")
)

write.csv(
  as.data.frame(
    top_genes_by_cluster,
    stringsAsFactors = FALSE
  ),
  file = top_gene_file,
  row.names = FALSE
)

# 1.8 将各cluster的基因由SYMBOL转换为ENTREZID

top_genes_entrez <- lapply(
  top_genes_by_cluster,
  function(gene_vector) {
    gene_vector <- as.character(gene_vector)
    gene_vector <- gene_vector[!is.na(gene_vector)]
    gene_vector <- trimws(gene_vector)
    gene_vector <- gene_vector[gene_vector != ""]
    
    if (length(gene_vector) == 0) {
      return(character(0))
    }
    
    gene_df <- tryCatch(
      {
        bitr(
          gene_vector,
          fromType = "SYMBOL",
          toType = "ENTREZID",
          OrgDb = orgdb_object
        )
      },
      error = function(e) {
        NULL
      }
    )
    
    if (is.null(gene_df) || nrow(gene_df) == 0) {
      return(character(0))
    }
    
    unique(
      as.character(gene_df$ENTREZID)
    )
  }
)

# 1.9 获取可进行富集分析的cluster

available_clusters <- names(top_genes_entrez)

if (length(available_clusters) == 0) {
  stop("没有获得可用于富集分析的cluster。")
}

available_cluster_text <- paste0(
  "当前可选cluster：\n",
  paste(available_clusters, collapse = ", ")
)

prepare_summary <- paste0(
  "Top基因准备完成。\n",
  "基因方向：", gene_direction, "\n",
  "每个cluster保留Top基因数：", top_gene_number, "\n",
  "可选cluster数：", length(available_clusters), "\n",
  "可选cluster：", paste(available_clusters, collapse = ", "), "\n",
  "Top基因表已保存至：", top_gene_file, "\n",
  "当前阶段已完成Top基因筛选与ENTREZID转换。"
)

# 2. 对选定cluster进行GO和KEGG富集分析

# 2.1 参数设置

# 默认选择第一个cluster，可根据实际情况修改
selected_celltype <- available_clusters[1]

# KEGG物种代码，人类为hsa
kegg_organism_code <- "hsa"

# GO和KEGG结果保留阈值
pvalue_cutoff <- 1
qvalue_cutoff <- 1

# 富集图显示的通路数量
showCategory_go <- 5
showCategory_kegg <- 15

# 可选pvalue、p.adjust或qvalue
color_by <- "pvalue"

# GO图保存参数
w_go <- 8
h_go <- 8
name_go_pdf <- "GO_dotplot"

# KEGG图保存参数
w_kegg <- 8
h_kegg <- 8
name_kegg_pdf <- "KEGG_dotplot"

# 2.2 检查选择的cluster

if (!selected_celltype %in% available_clusters) {
  stop(
    paste0(
      "selected_celltype 不存在，可选cluster为：",
      paste(available_clusters, collapse = ", ")
    )
  )
}

gene_list <- top_genes_entrez[[selected_celltype]]

if (is.null(gene_list) || length(gene_list) == 0) {
  stop(
    paste0(
      "cluster ",
      selected_celltype,
      " 没有可用于富集分析的ENTREZID。"
    )
  )
}

# 2.3 运行GO富集分析

go_result <- enrichGO(
  gene = gene_list,
  OrgDb = orgdb_object,
  ont = "ALL",
  pvalueCutoff = pvalue_cutoff,
  qvalueCutoff = qvalue_cutoff
)

go_result_exists <- !is.null(go_result) &&
  nrow(as.data.frame(go_result)) > 0

if (go_result_exists) {
  go_result <- setReadable(
    go_result,
    OrgDb = orgdb_object,
    keyType = "ENTREZID"
  )
  
  go_result_file <- file.path(
    out_dir,
    paste0(selected_celltype, "_go_cell.csv")
  )
  
  write.csv(
    as.data.frame(go_result),
    file = go_result_file,
    row.names = FALSE
  )
  
  go_plot <- enrichplot::dotplot(
    go_result,
    showCategory = showCategory_go,
    split = "ONTOLOGY",
    label_format = 50,
    color = color_by
  ) +
    facet_grid(
      ONTOLOGY ~ .,
      scales = "free"
    )
  
  print(go_plot)
  
  go_pdf_file <- file.path(
    out_dir,
    paste0(
      selected_celltype,
      "_",
      name_go_pdf,
      ".pdf"
    )
  )
  
  ggsave(
    filename = go_pdf_file,
    plot = go_plot,
    width = w_go,
    height = h_go,
    device = "pdf"
  )
} else {
  go_result_file <- NA_character_
  go_pdf_file <- NA_character_
  
  warning(
    paste0(
      "cluster ",
      selected_celltype,
      " 未获得GO富集结果。"
    )
  )
}

# 2.4 运行KEGG富集分析

kegg_result <- enrichKEGG(
  gene = gene_list,
  organism = kegg_organism_code,
  pvalueCutoff = pvalue_cutoff,
  qvalueCutoff = qvalue_cutoff
)

kegg_result_exists <- !is.null(kegg_result) &&
  nrow(as.data.frame(kegg_result)) > 0

if (kegg_result_exists) {
  kegg_result <- setReadable(
    kegg_result,
    OrgDb = orgdb_object,
    keyType = "ENTREZID"
  )
  
  kegg_result_file <- file.path(
    out_dir,
    paste0(selected_celltype, "_kegg_cell.csv")
  )
  
  write.csv(
    as.data.frame(kegg_result),
    file = kegg_result_file,
    row.names = FALSE
  )
  
  kegg_plot <- enrichplot::dotplot(
    kegg_result,
    showCategory = showCategory_kegg,
    label_format = 50,
    color = color_by
  ) +
    theme(
      axis.text.x = element_text(
        angle = 30,
        hjust = 1
      )
    )
  
  print(kegg_plot)
  
  kegg_pdf_file <- file.path(
    out_dir,
    paste0(
      selected_celltype,
      "_",
      name_kegg_pdf,
      ".pdf"
    )
  )
  
  ggsave(
    filename = kegg_pdf_file,
    plot = kegg_plot,
    width = w_kegg,
    height = h_kegg,
    device = "pdf"
  )
} else {
  kegg_result_file <- NA_character_
  kegg_pdf_file <- NA_character_
  
  warning(
    paste0(
      "cluster ",
      selected_celltype,
      " 未获得KEGG富集结果。"
    )
  )
}

# 2.5 生成富集分析结果说明

enrichment_summary <- paste0(
  "富集分析完成。\n",
  "本次选择cluster：", selected_celltype, "\n",
  "用于富集分析的ENTREZID数量：", length(gene_list), "\n",
  "GO结果是否存在：", ifelse(go_result_exists, "是", "否"), "\n",
  "KEGG结果是否存在：", ifelse(kegg_result_exists, "是", "否"), "\n",
  "CSV和PDF结果已保存到目录：", out_dir
)

# 3. 保存分析参数记录

# 3.1 参数设置

name_params <- "enrichment_parameters"

# 3.2 生成参数记录

param_text <- paste0(
  "本次分析参数总结：\n\n",
  
  "【第1步：Top基因准备】\n",
  "- 基因方向：", gene_direction, "\n",
  "- 每个cluster选取Top基因数：", top_gene_number, "\n",
  "- Top基因表：", basename(top_gene_file), "\n",
  "- 可选cluster：", paste(available_clusters, collapse = ", "), "\n\n",
  
  "【第2步：富集分析参数】\n",
  "- OrgDb对象名：", orgdb_object, "\n",
  "- 本次选择cluster：", selected_celltype, "\n",
  "- 用于富集分析的ENTREZID数量：", length(gene_list), "\n",
  "- KEGG物种代码：", kegg_organism_code, "\n",
  "- pvalueCutoff：", pvalue_cutoff, "\n",
  "- qvalueCutoff：", qvalue_cutoff, "\n",
  "- GO显示通路数：", showCategory_go, "\n",
  "- KEGG显示通路数：", showCategory_kegg, "\n",
  "- dotplot颜色映射变量：", color_by, "\n",
  "- GO保存宽：", w_go, "英寸\n",
  "- GO保存高：", h_go, "英寸\n",
  "- KEGG保存宽：", w_kegg, "英寸\n",
  "- KEGG保存高：", h_kegg, "英寸\n",
  "- 输出目录：", out_dir, "\n\n",
  
  "写作提示词（自行组装语言，或借助AI组装）：\n",
  "1. 基于差异基因结果，按指定方向筛选各cluster的Top基因。\n",
  "2. 将筛选得到的基因由SYMBOL转换为ENTREZID，为后续富集分析做准备。\n",
  "3. 在Top基因筛选基础上，对选定cluster分别进行GO与KEGG富集分析。\n",
  "4. 比较该cluster在生物过程、分子功能、细胞组分及经典通路上的富集特征。\n",
  "5. 结合富集结果解析该区域或细胞亚群的潜在功能状态与生物学意义。"
)

# 3.3 保存参数记录

param_file <- file.path(
  out_dir,
  paste0(name_params, ".txt")
)

writeLines(
  param_text,
  con = param_file
)