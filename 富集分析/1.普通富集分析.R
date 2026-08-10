library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(R.utils)

# 1. 设置输入参数
input_file <- "差异分析结果.csv"
output_folder <- "富集分析"

go_orgdb <- "org.Hs.eg.db"
kegg_species <- "hsa"

kegg_show_num <- 15
go_show_num <- 5
pvalue_display <- "pvalue"

gsea_gene_set_id <- 1
gsea_colors <- c("#FFCC33", "#333399", "#FF0033")

# 2. 创建总输出文件夹
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

# 3. 读取输入数据
Genes_All <- read.csv(
  input_file,
  header = TRUE,
  check.names = FALSE,
  row.names = NULL
)

# 4. 根据输入文件列数统一列名
input_col_num <- ncol(Genes_All)

if (input_col_num == 1) {
  colnames(Genes_All)[1] <- "gene_symbol"
} else if (input_col_num == 2) {
  colnames(Genes_All)[1:2] <- c("gene_symbol", "change")
} else if (input_col_num >= 3) {
  colnames(Genes_All)[1:3] <- c("gene_symbol", "change", "logFC")
} else {
  stop("输入文件至少需要包含 gene_symbol 这一列。")
}

# 5. 设置 clusterProfiler 下载方式
R.utils::setOption("clusterProfiler.download.method", "auto")

# 6. 定义 GO + KEGG 富集分析通用流程
run_kegg_go_enrichment <- function(genes_input, output_dir, message_text) {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  cat(message_text, "\n")
  
  gene <- suppressWarnings(
    suppressMessages(
      bitr(
        genes_input$gene_symbol,
        fromType = "SYMBOL",
        toType = "ENTREZID",
        OrgDb = go_orgdb
      )
    )
  )
  
  GO <- enrichGO(
    gene$ENTREZID,
    OrgDb = go_orgdb,
    keyType = "ENTREZID",
    ont = "ALL",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    readable = TRUE
  )
  
  cat("\nGO富集分析已完成\n")
  write.csv(
    GO@result,
    file = file.path(output_dir, "GO_result.csv"),
    row.names = FALSE
  )
  
  KEGG <- enrichKEGG(
    gene$ENTREZID,
    organism = kegg_species,
    pvalueCutoff = 1,
    qvalueCutoff = 1
  )
  
  cat("\nKEGG富集分析已完成\n")
  
  KEGG <- setReadable(
    KEGG,
    OrgDb = go_orgdb,
    keyType = "ENTREZID"
  )
  
  write.csv(
    KEGG@result,
    file = file.path(output_dir, "KEGG_result.csv"),
    row.names = FALSE
  )
  
  KEGG@result$Description <- gsub(
    " - Mus musculus \\(house mouse\\)$",
    "",
    KEGG@result$Description
  )
  
  p <- suppressWarnings(
    barplot(
      GO,
      showCategory = go_show_num,
      split = "ONTOLOGY",
      color = pvalue_display,
      font.size = 8
    ) +
      facet_grid(ONTOLOGY ~ ., scale = "free")
  )
  ggsave(
    file.path(output_dir, "1.GO柱状图.pdf"),
    plot = p,
    width = 15,
    height = 16,
    units = "cm"
  )
  message("1.GO柱状图已绘制完成，并保存在文件夹中")
  
  p <- suppressWarnings(
    barplot(
      KEGG,
      showCategory = kegg_show_num,
      font.size = 8,
      color = pvalue_display,
      label_format = 40
    )
  )
  ggsave(
    file.path(output_dir, "2.KEGG柱状图.pdf"),
    plot = p,
    width = 15,
    height = 16,
    units = "cm"
  )
  message("2.KEGG柱状图已绘制完成，并保存在文件夹中")
  
  p <- dotplot(
    GO,
    showCategory = go_show_num,
    split = "ONTOLOGY",
    color = pvalue_display,
    font.size = 8
  ) +
    facet_grid(ONTOLOGY ~ ., scale = "free")
  ggsave(
    file.path(output_dir, "3.GO气泡图.pdf"),
    plot = p,
    width = 15,
    height = 16,
    units = "cm"
  )
  message("3.GO气泡图已绘制完成，并保存在文件夹中")
  
  p <- dotplot(
    KEGG,
    showCategory = kegg_show_num,
    font.size = 8,
    color = pvalue_display
  )
  ggsave(
    file.path(output_dir, "4.KEGG气泡图.pdf"),
    plot = p,
    width = 15,
    height = 16,
    units = "cm"
  )
  message("4.KEGG气泡图已绘制完成，并保存在文件夹中")
  
  return(list(
    GO = GO,
    KEGG = KEGG,
    Genes_Input = genes_input,
    gene = gene
  ))
}

# 7. 初始化结果列表
result410 <- list()

# 8. 根据输入列数准备 All genes 数据
if (input_col_num == 1) {
  Genes_All_for_enrich <- Genes_All
} else {
  Genes_All_for_enrich <- Genes_All[Genes_All$change != "NOT", ]
}

# 9. 运行 All genes 的 GO + KEGG 富集分析
all_dir <- file.path(output_folder, "All_Gene")

result410$All_Gene <- tryCatch({
  run_kegg_go_enrichment(
    genes_input = Genes_All_for_enrich,
    output_dir = all_dir,
    message_text = "开始KEGG-GO全部基因或全部差异基因富集分析······"
  )
}, error = function(e) {
  message("All genes 富集分析出错：", e$message)
  NULL
})

# 10. 如果输入文件有 change 列，则额外运行 UP 和 DOWN 分析
if (input_col_num >= 2) {
  
  up_dir <- file.path(output_folder, "UP_Gene")
  down_dir <- file.path(output_folder, "DOWN_Gene")
  
  Genes_UP <- Genes_All[Genes_All$change != "NOT", ]
  Genes_UP <- Genes_UP[Genes_UP$change != "DOWN", ]
  
  result410$UP_Gene <- tryCatch({
    run_kegg_go_enrichment(
      genes_input = Genes_UP,
      output_dir = up_dir,
      message_text = "开始KEGG-GO上调基因富集分析······"
    )
  }, error = function(e) {
    message("UP genes 富集分析出错：", e$message)
    NULL
  })
  
  Genes_DOWN <- Genes_All[Genes_All$change != "NOT", ]
  Genes_DOWN <- Genes_DOWN[Genes_DOWN$change != "UP", ]
  
  result410$DOWN_Gene <- tryCatch({
    run_kegg_go_enrichment(
      genes_input = Genes_DOWN,
      output_dir = down_dir,
      message_text = "开始KEGG-GO下调基因富集分析······"
    )
  }, error = function(e) {
    message("DOWN genes 富集分析出错：", e$message)
    NULL
  })
}

# 11. 如果输入文件有 logFC 列，则额外运行 GSEA 分析
if (input_col_num >= 3) {
  
  gsea_dir <- file.path(output_folder, "GSEA")
  
  if (!dir.exists(gsea_dir)) {
    dir.create(gsea_dir, recursive = TRUE)
  }
  
  result410$GSEA <- tryCatch({
    
    Genes_GSEA <- Genes_All
    colnames(Genes_GSEA)[1:3] <- c("gene_symbol", "change", "logFC")
    
    gene <- suppressWarnings(
      suppressMessages(
        bitr(
          Genes_GSEA$gene_symbol,
          fromType = "SYMBOL",
          toType = "ENTREZID",
          OrgDb = go_orgdb
        )
      )
    )
    
    info <- Genes_GSEA[, c("gene_symbol", "logFC")]
    names(info) <- c("SYMBOL", "Log2FoldChange")
    
    info_merge <- merge(info, gene, by = "SYMBOL")
    
    GSEA_input <- info_merge$Log2FoldChange
    names(GSEA_input) <- info_merge$ENTREZID
    
    cat("GSEA-KEGG分析开始······\n")
    
    GSEA_input <- sort(GSEA_input, decreasing = TRUE)
    
    KEGG_ges <- suppressWarnings(
      gseKEGG(
        GSEA_input,
        organism = kegg_species,
        pvalueCutoff = 1,
        eps = 0
      )
    )
    
    KEGG_ges <- setReadable(
      KEGG_ges,
      OrgDb = go_orgdb,
      keyType = "ENTREZID"
    )
    
    KEGG_ges_result <- KEGG_ges@result
    
    write.csv(
      KEGG_ges_result,
      file = file.path(gsea_dir, "GSEA-KEGG.csv")
    )
    
    p <- ridgeplot(
      KEGG_ges,
      showCategory = as.numeric(kegg_show_num),
      fill = pvalue_display,
      decreasing = TRUE
    ) +
      theme(
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8)
      )
    
    ggsave(
      file.path(gsea_dir, "1.GSEA-KEGG山峦图.pdf"),
      plot = p,
      width = 15,
      height = 15,
      units = "cm"
    )
    message("GSEA-KEGG山峦图绘制完成，并保存在文件夹中")
    
    pp <- dotplot(
      KEGG_ges,
      showCategory = kegg_show_num,
      font.size = 8,
      color = pvalue_display
    )
    
    ggsave(
      file.path(gsea_dir, "2.GSEA-KEGG气泡图.pdf"),
      plot = pp,
      width = 15,
      height = 15,
      units = "cm"
    )
    message("GSEA-KEGG气泡图绘制完成，并保存在文件夹中")
    
    p5 <- gseaplot2(
      KEGG_ges,
      geneSetID = gsea_gene_set_id,
      color = gsea_colors,
      pvalue_table = T,
      pvalue_table_columns = c("pvalue", "NES"),
      ES_geom = "line"
    )
    
    ggsave(
      file.path(gsea_dir, "3.GSEA-KEGG主图.pdf"),
      plot = p5,
      width = 30,
      height = 30,
      units = "cm"
    )
    message("GSEA-KEGG主图绘制完成，并保存在文件夹中")
    
    cat("\nGSEA-GO分析开始······\n")
    
    GO_ges <- suppressWarnings(
      gseGO(
        geneList = GSEA_input,
        OrgDb = go_orgdb,
        ont = "ALL",
        pvalueCutoff = 1,
        eps = 0
      )
    )
    
    GO_ges <- setReadable(
      GO_ges,
      OrgDb = go_orgdb,
      keyType = "ENTREZID"
    )
    
    GO_ges_result <- GO_ges@result
    
    write.csv(
      GO_ges_result,
      file = file.path(gsea_dir, "GSEA-GO.csv")
    )
    
    p <- ridgeplot(
      GO_ges,
      showCategory = go_show_num * 3,
      fill = pvalue_display,
      decreasing = TRUE
    ) +
      theme(
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8)
      )
    
    ggsave(
      file.path(gsea_dir, "4.GSEA-GO山峦图.pdf"),
      plot = p,
      width = 15,
      height = 15,
      units = "cm"
    )
    message("GSEA-GO山峦图绘制完成，并保存在文件夹中")
    
    pp <- dotplot(
      GO_ges,
      showCategory = go_show_num,
      split = "ONTOLOGY",
      color = pvalue_display,
      font.size = 8
    ) +
      facet_grid(ONTOLOGY ~ ., scale = "free")
    
    ggsave(
      file.path(gsea_dir, "5.GSEA-GO气泡图.pdf"),
      plot = pp,
      width = 15,
      height = 15,
      units = "cm"
    )
    message("GSEA-GO气泡图绘制完成，并保存在文件夹中")
    
    p5 <- gseaplot2(
      GO_ges,
      geneSetID = gsea_gene_set_id,
      color = gsea_colors,
      pvalue_table = TRUE,
      ES_geom = "line"
    )
    
    ggsave(
      file.path(gsea_dir, "6.GSEA-GO主图.pdf"),
      plot = p5,
      width = 30,
      height = 30,
      units = "cm"
    )
    message("GSEA-GO主图绘制完成，并保存在文件夹中")
    
    list(
      GO_ges = GO_ges,
      KEGG_ges = KEGG_ges,
      KEGG_ges_result = KEGG_ges_result,
      GO_ges_result = GO_ges_result
    )
    
  }, error = function(e) {
    message("GSEA 分析出错：", e$message)
    NULL
  })
}

# 12. 输出完成提示
cat("富集分析流程已完成，结果保存在：", output_folder, "\n")

