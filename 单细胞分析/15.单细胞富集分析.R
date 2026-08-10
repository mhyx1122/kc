suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(clusterProfiler)
  library(enrichplot)
  library(AnnotationDbi)
})
options(timeout = 600)
# 1. 检查全局环境中是否存在 all_markers 对象

if (!exists("all_markers", envir = .GlobalEnv)) {
  stop("全局环境中没有 all_markers 对象，请先准备差异分析结果。")
}

all_markers <- get("all_markers", envir = .GlobalEnv)


# 2. 输出目录模块

# 2.1 输出目录参数
out_dir <- "6.KEGG-GO富集分析"

# 2.2 创建输出目录
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}


# 3. Top 基因准备参数模块

# 3.1 基因方向参数
# 可选："up" 或 "down"
gene_direction <- "up"

# 3.2 每个 cluster 选取 Top 基因数
top_gene_number <- 200

# 3.3 Top 基因表保存文件名
name_top_gene_csv <- "top200gene"


# 4. all_markers 检查模块

# 4.1 检查必要列
required_cols <- c("gene", "cluster", "avg_log2FC")

missing_cols <- setdiff(required_cols, colnames(all_markers))

if (length(missing_cols) > 0) {
  stop(paste0("all_markers 缺少必要列：", paste(missing_cols, collapse = ", ")))
}


# 5. 筛选 Top 基因模块

# 5.1 按方向筛选基因
if (gene_direction == "up") {
  
  selected_markers <- all_markers %>%
    dplyr::filter(avg_log2FC > 0)
  
} else if (gene_direction == "down") {
  
  selected_markers <- all_markers %>%
    dplyr::filter(avg_log2FC < 0)
  
} else {
  
  stop("gene_direction 只能设置为 up 或 down")
}

if (nrow(selected_markers) == 0) {
  stop("筛选后的 selected_markers 为空，请检查 all_markers 或 gene_direction 参数。")
}

# 5.2 每个 cluster 选择 Top 基因
selected_markers <- selected_markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(
    order_by = abs(avg_log2FC),
    n = top_gene_number,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup()

# 5.3 按 cluster 拆分基因
top_genes_by_cluster <- unstack(
  selected_markers,
  gene ~ cluster
)

# 5.4 保存 Top 基因表
write.csv(
  as.data.frame(top_genes_by_cluster, stringsAsFactors = FALSE),
  file = file.path(out_dir, paste0(name_top_gene_csv, ".csv")),
  row.names = FALSE
)


# 6. 基因 ID 转换参数模块

# 6.1 OrgDb 对象名
orgdb_object_name <- "org.Hs.eg.db"

# 6.2 KEGG 物种代码
kegg_organism_code <- "hsa"

# 6.3 加载 OrgDb 对象
if (!requireNamespace(orgdb_object_name, quietly = TRUE)) {
  stop(paste0("缺少 ", orgdb_object_name, " 包，请先安装。"))
}

OrgDb_use <- get(
  orgdb_object_name,
  envir = asNamespace(orgdb_object_name)
)


# 7. SYMBOL 转 ENTREZID 模块

top_genes_entrez <- lapply(top_genes_by_cluster, function(x) {
  
  x <- x[!is.na(x) & x != ""]
  
  if (length(x) == 0) {
    return(character(0))
  }
  
  gene_df <- tryCatch({
    bitr(
      x,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = OrgDb_use
    )
  }, error = function(e) {
    NULL
  })
  
  if (is.null(gene_df) || nrow(gene_df) == 0) {
    return(character(0))
  }
  
  unique(as.character(gene_df$ENTREZID))
})

available_clusters <- names(top_genes_entrez)


# 8. 富集分析参数模块

# 8.1 选择需要富集分析的 cluster
# 默认选择第一个 cluster
selected_celltype <- available_clusters[1]

# 如果要指定某个 cluster，可以改成类似：
# selected_celltype <- "0"
# selected_celltype <- "T cells"

# 8.2 GO / KEGG 显示通路数
showCategory_go <- 5
showCategory_kegg <- 15

# 8.3 dotplot 颜色映射变量
# 可选："pvalue"、"p.adjust"、"qvalue"
color_by <- "pvalue"

# 8.4 GO 图保存参数
w_go <- 8
h_go <- 8
name_go_pdf <- "GO_dotplot"

# 8.5 KEGG 图保存参数
w_kegg <- 8
h_kegg <- 8
name_kegg_pdf <- "KEGG_dotplot"


# 9. 富集分析前检查模块

if (is.null(selected_celltype) || selected_celltype == "") {
  stop("selected_celltype 不能为空。")
}

if (!selected_celltype %in% names(top_genes_entrez)) {
  stop(paste0("selected_celltype 不在可选 cluster 中：", selected_celltype))
}

gene_list <- top_genes_entrez[[selected_celltype]]

if (is.null(gene_list) || length(gene_list) == 0) {
  stop(paste0("cluster ", selected_celltype, " 没有可用于富集分析的 ENTREZID。"))
}


# 10. GO 富集分析模块

# 10.1 运行 GO 富集分析
go_result <- enrichGO(
  gene = gene_list,
  OrgDb = OrgDb_use,
  ont = "ALL",
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

# 10.2 保存 GO 富集结果
if (!is.null(go_result) && nrow(as.data.frame(go_result)) > 0) {
  
  go_result <- setReadable(
    go_result,
    OrgDb = OrgDb_use,
    keyType = "ENTREZID"
  )
  
  write.csv(
    as.data.frame(go_result),
    file = file.path(out_dir, paste0(selected_celltype, "_go_cell.csv")),
    row.names = FALSE
  )
  
  p_go <- dotplot(
    go_result,
    showCategory = showCategory_go,
    split = "ONTOLOGY",
    label_format = 50,
    color = color_by
  ) +
    facet_grid(ONTOLOGY ~ ., scales = "free")
  
  ggsave(
    filename = file.path(out_dir, paste0(selected_celltype, "_", name_go_pdf, ".pdf")),
    plot = p_go,
    width = w_go,
    height = h_go,
    device = "pdf"
  )
}


# 11. KEGG 富集分析模块

# 11.1 运行 KEGG 富集分析
kegg_result <- enrichKEGG(
  gene = gene_list,
  organism = kegg_organism_code,
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

# 11.2 保存 KEGG 富集结果
if (!is.null(kegg_result) && nrow(as.data.frame(kegg_result)) > 0) {
  
  kegg_result <- setReadable(
    kegg_result,
    OrgDb = OrgDb_use,
    keyType = "ENTREZID"
  )
  
  write.csv(
    as.data.frame(kegg_result),
    file = file.path(out_dir, paste0(selected_celltype, "_kegg_cell.csv")),
    row.names = FALSE
  )
  
  p_kegg <- dotplot(
    kegg_result,
    showCategory = showCategory_kegg,
    label_format = 50,
    color = color_by
  ) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  
  ggsave(
    filename = file.path(out_dir, paste0(selected_celltype, "_", name_kegg_pdf, ".pdf")),
    plot = p_kegg,
    width = w_kegg,
    height = h_kegg,
    device = "pdf"
  )
}


# 12. 参数记录模块

# 12.1 参数文件保存参数
name_params <- "enrichment_parameters"

# 12.2 生成参数记录文本
param_text <- paste0(
  "本次分析参数总结：\n",
  "【第1步：Top基因准备】\n",
  "- 基因方向：", gene_direction, "\n",
  "- 每个cluster选取Top基因数：", top_gene_number, "\n",
  "- 可选cluster：", paste(available_clusters, collapse = ", "), "\n",
  "【第2步：富集分析参数】\n",
  "- OrgDb对象名：", orgdb_object_name, "\n",
  "- 本次选择cluster：", selected_celltype, "\n",
  "- KEGG物种代码：", kegg_organism_code, "\n",
  "- GO显示通路数：", showCategory_go, "\n",
  "- KEGG显示通路数：", showCategory_kegg, "\n",
  "- dotplot颜色映射变量：", color_by, "\n",
  "- GO保存宽：", w_go, " 英寸\n",
  "- GO保存高：", h_go, " 英寸\n",
  "- KEGG保存宽：", w_kegg, " 英寸\n",
  "- KEGG保存高：", h_kegg, " 英寸\n",
  "- 输出目录：", out_dir, "\n",
  "写作提示词（自行组装语言，或借助AI组装）：\n",
  "1.基于差异基因结果，按指定方向筛选各cluster的Top基因。\n",
  "2.将筛选得到的基因由SYMBOL转换为ENTREZID，为后续富集分析做准备。\n",
  "3.在Top基因筛选基础上，对选定cluster分别进行GO与KEGG富集分析。\n",
  "4.结合富集结果可进一步解析该细胞亚群的潜在功能状态与生物学意义。"
)

# 12.3 保存参数记录
writeLines(
  param_text,
  con = file.path(out_dir, paste0(name_params, ".txt"))
)