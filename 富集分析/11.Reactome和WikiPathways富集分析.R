library(clusterProfiler)
library(enrichplot)
library(ReactomePA)
library(ggplot2)
library(R.utils)

# 1. 设置输入参数
input_file <- "基因文件.csv"
out_dir <- "Reactome和WikiPathways"
species <- "human"
color_by <- "pvalue"
show_n <- 15
label_format <- 50
palette_reactome <- "#E64B35FF,#4DBBD5FF,#00A087FF,#3C5488FF,#F39B7FFF,#8491B4FF,#91D1C2FF,#7E6148FF,#0073C2FF,#EFC000FF,#868686FF,#CD534CFF,#7AA6DCFF,#003C67FF,#8F7700FF,#3B3B3BFF,#A73030FF,#374E55FF,#DF8F44FF,#00A1D5FF,#B24745FF,#79AF97FF,#6A6599FF,#80796BFF,#BC3C29FF,#0072B5FF,#E18727FF,#20854EFF,#7876B1FF,#6F99ADFF,#FFDC91FF,#EE4C97FF"
palette_wp <- "#E64B35FF,#4DBBD5FF,#00A087FF,#3C5488FF,#F39B7FFF,#8491B4FF,#91D1C2FF,#7E6148FF,#0073C2FF,#EFC000FF,#868686FF,#CD534CFF,#7AA6DCFF,#003C67FF,#8F7700FF,#3B3B3BFF,#A73030FF,#374E55FF,#DF8F44FF,#00A1D5FF,#B24745FF,#79AF97FF,#6A6599FF,#80796BFF,#BC3C29FF,#0072B5FF,#E18727FF,#20854EFF,#7876B1FF,#6F99ADFF,#FFDC91FF,#EE4C97FF"
w_reactome <- 8
h_reactome <- 8
name_reactome <- "all_reactome_cell"
w_wp <- 8
h_wp <- 8
name_wp <- "all_wikipathways_cell"
name_params <- "Reactome_WikiPathways_parameters"

# 2. 创建输出文件夹
if (!dir.exists(out_dir)) {
  dir.create(out_dir,recursive = TRUE)
}

# 3. 根据物种设置分析参数
if (species == "human") {
  reactome_organism <- "human"
  wp_organism <- "Homo sapiens"
  orgdb_name <- "org.Hs.eg.db"
  species_cn <- "人"
} else if (species == "mouse") {
  reactome_organism <- "mouse"
  wp_organism <- "Mus musculus"
  orgdb_name <- "org.Mm.eg.db"
  species_cn <- "小鼠"
} else if (species == "rat") {
  reactome_organism <- "rat"
  wp_organism <- "Rattus norvegicus"
  orgdb_name <- "org.Rn.eg.db"
  species_cn <- "大鼠"
} else {
  stop("不支持的物种类型")
}

# 4. 加载物种注释包
if (!requireNamespace(orgdb_name,quietly = TRUE)) {
  stop(
    paste0(
      "缺少注释包：",orgdb_name,
      "。请先安装，例如 BiocManager::install('",orgdb_name,"')"
    )
  )
}
orgdb <- get(orgdb_name,envir = asNamespace(orgdb_name))

# 5. 读取输入文件
deg_df <- read.csv(
  input_file,
  header = TRUE,
  check.names = FALSE,
  row.names = NULL
)

# 6. 检查 gene_symbol 列
if (!"gene_symbol" %in% colnames(deg_df)) {
  stop("输入文件中必须包含 gene_symbol 列。")
}

# 7. 基因 SYMBOL 转换为 ENTREZID
gene_id_df <- suppressWarnings(
  suppressMessages(
    bitr(
      deg_df$gene_symbol,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = orgdb
    )
  )
)
gene_vector <- unique(gene_id_df$ENTREZID)
if (length(gene_vector) == 0) {
  stop("未成功转换出 ENTREZID，请检查 gene_symbol 和物种是否匹配。")
}

# 8. 设置 clusterProfiler 下载方式
R.utils::setOption("clusterProfiler.download.method","auto")

# 9. Reactome 富集分析
reactome_res <- enrichPathway(
  gene = gene_vector,
  organism = reactome_organism,
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  readable = FALSE
)
reactome_res <- setReadable(
  reactome_res,
  OrgDb = orgdb,
  keyType = "ENTREZID"
)

# 10. 保存 Reactome 富集结果表
write.csv(
  as.data.frame(reactome_res),
  file = file.path(out_dir,"all_reactome_cell.csv"),
  row.names = FALSE
)

# 11. WikiPathways 富集分析
wp_res <- enrichWP(
  gene = gene_vector,
  organism = wp_organism,
  pvalueCutoff = 1,
  qvalueCutoff = 1
)
wp_res <- setReadable(
  wp_res,
  OrgDb = orgdb,
  keyType = "ENTREZID"
)

# 12. 保存 WikiPathways 富集结果表
write.csv(
  as.data.frame(wp_res),
  file = file.path(out_dir,"all_wikipathways_cell.csv"),
  row.names = FALSE
)

# 13. 处理 Reactome 颜色集合
pal_reactome <- unlist(strsplit(palette_reactome,","))
pal_reactome <- trimws(pal_reactome)
pal_reactome <- pal_reactome[pal_reactome != ""]

# 14. 绘制 Reactome 气泡图
p_reactome <- dotplot(
  reactome_res,
  showCategory = show_n,
  label_format = label_format,
  color = color_by
) +
  theme(axis.text.x = element_text(angle = 30,hjust = 1))
if (length(pal_reactome) >= 2) {
  p_reactome <- p_reactome + scale_color_gradientn(colours = pal_reactome)
}

# 15. 保存 Reactome 气泡图
ggsave(
  filename = file.path(out_dir,paste0(name_reactome,".pdf")),
  plot = p_reactome,
  width = w_reactome,
  height = h_reactome,
  device = "pdf"
)

# 16. 处理 WikiPathways 颜色集合
pal_wp <- unlist(strsplit(palette_wp,","))
pal_wp <- trimws(pal_wp)
pal_wp <- pal_wp[pal_wp != ""]

# 17. 绘制 WikiPathways 气泡图
p_wp <- dotplot(
  wp_res,
  showCategory = show_n,
  label_format = label_format,
  color = color_by
) +
  theme(axis.text.x = element_text(angle = 30,hjust = 1))
if (length(pal_wp) >= 2) {
  p_wp <- p_wp + scale_color_gradientn(colours = pal_wp)
}

# 18. 保存 WikiPathways 气泡图
ggsave(
  filename = file.path(out_dir,paste0(name_wp,".pdf")),
  plot = p_wp,
  width = w_wp,
  height = h_wp,
  device = "pdf"
)

# 19. 生成参数记录文本
param_text <- paste0(
  "本次分析参数总结：\n",
  "- 输入文件：",basename(input_file),"\n",
  "- 物种：",species_cn,"\n",
  "- Reactome organism：",reactome_organism,"\n",
  "- WikiPathways organism：",wp_organism,"\n",
  "- OrgDb：",orgdb_name,"\n",
  "- 显示前 n 个通路：",show_n,"\n",
  "- 颜色映射变量：",color_by,"\n",
  "- 标签最大字符数：",label_format,"\n",
  "说明：\n",
  "1. 使用 gene_symbol 列进行 SYMBOL -> ENTREZID 转换。\n",
  "2. 分别进行了 Reactome 和 WikiPathways 富集分析。\n",
  "3. 可根据 pvalue 或 p.adjust 进行着色显示，并支持手动修改颜色集合。"
)

# 20. 保存参数记录文件
if (is.null(name_params) || trimws(name_params) == "") {
  name_params <- "Reactome_WikiPathways_parameters"
}
writeLines(
  param_text,
  con = file.path(out_dir,paste0(name_params,".txt"))
)

# 21. 输出完成提示
message("Reactome 和 WikiPathways 富集分析完成，结果已保存到文件夹：",out_dir)