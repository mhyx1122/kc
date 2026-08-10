# 1. 加载必要的 R 包

suppressPackageStartupMessages({
  library(Seurat)
  library(SCopeLoomR)
  library(AUCell)
  library(SCENIC)
  library(arrow)
  library(RcisTarget)
  library(foreach)
})


# 2. 清理内存和连接

gc()

closeAllConnections()

options(max.connections = 5000)


# 3. 检查 Seurat 对象

if (!exists("seurat", envir = .GlobalEnv)) {
  stop("全局环境中没有找到 seurat 对象。")
}

if (!inherits(seurat, "Seurat")) {
  stop("全局环境中的 seurat 不是 Seurat 对象。")
}

if (!"RNA" %in% names(seurat@assays)) {
  stop("seurat 对象中不存在 RNA assay。")
}

if (!"cellType" %in% colnames(seurat@meta.data)) {
  stop("seurat@meta.data 中不存在 cellType 列。")
}


# 4. 设置基因过滤参数

# 可选：
# filter_method <- "pct"
# filter_method <- "hvg"
filter_method <- "hvg"

# filter_method = "pct" 时使用
min_pct <- 0.05

# filter_method = "hvg" 时使用
top_hvg_n <- 2000


# 5. 设置物种参数

# 人类使用 hgnc
# 小鼠使用 mgi
org <- "hgnc"


# 6. 设置并行核心数

nCores410 <- 1


# 7. 设置 SCENIC 数据库文件

db_prefix <- ifelse(
  org == "hgnc",
  "hg38",
  "mm10"
)

mm_dbs <- list(
  "500bp" = paste0(
    db_prefix,
    "_500bp.feather"
  ),
  "10kb" = paste0(
    db_prefix,
    "_10kb.feather"
  )
)


# 8. 创建结果文件夹

output_dir <- "SCENIC"

if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}


# 9. 设置细胞身份

Idents(seurat) <- "cellType"

seurat_object <- seurat

if (!"celltype" %in% colnames(seurat_object@meta.data)) {
  seurat_object@meta.data$celltype <- Idents(seurat_object)
}


# 10. 提取 RNA counts 表达矩阵

exprMat_raw <- as.matrix(
  GetAssayData(
    object = seurat_object,
    assay = "RNA",
    layer = "counts"
  )
)


# 11. 对基因进行预筛选

if (filter_method == "pct") {
  
  gene_keep_idx <- rowMeans(
    exprMat_raw != 0
  ) >= min_pct
  
  exprMat_sub <- exprMat_raw[
    gene_keep_idx,
    ,
    drop = FALSE
  ]
  
} else if (filter_method == "hvg") {
  
  seurat_tmp <- seurat_object
  
  DefaultAssay(seurat_tmp) <- "RNA"
  
  seurat_tmp <- FindVariableFeatures(
    object = seurat_tmp,
    selection.method = "vst",
    nfeatures = top_hvg_n,
    verbose = FALSE
  )
  
  hvg_genes <- VariableFeatures(seurat_tmp)
  
  hvg_genes <- intersect(
    hvg_genes,
    rownames(exprMat_raw)
  )
  
  exprMat_sub <- exprMat_raw[
    hvg_genes,
    ,
    drop = FALSE
  ]
  
} else {
  stop("filter_method 只能是 'pct' 或 'hvg'。")
}

if (nrow(exprMat_sub) == 0) {
  stop("基因过滤后没有剩余基因，请调整过滤参数。")
}


# 12. 整理细胞注释信息

cell_info <- seurat_object@meta.data

cell_info_subset <- cell_info[
  ,
  c(
    "cellType",
    "nCount_RNA",
    "nFeature_RNA"
  )
]

colnames(cell_info_subset) <- c(
  "CellType",
  "nGene",
  "nUMI"
)

cellInfo <- cell_info_subset


# 13. 初始化 SCENIC 参数

scenicOptions <- tryCatch(
  {
    initializeScenic(
      org = org,
      dbDir = "SCENIC_database",
      dbs = mm_dbs,
      datasetTitle = "expanalysis",
      nCores = nCores410
    )
  },
  error = function(e) {
    
    cat("APP作者v：cgxr410\n")
    
    motifAnnotations_hgnc <<- motifAnnotations
    motifAnnotations_mgi <<- motifAnnotations
    
    initializeScenic(
      org = org,
      dbDir = "SCENIC_database",
      dbs = mm_dbs,
      datasetTitle = "expanalysis",
      nCores = nCores410
    )
  }
)


# 14. 使用 SCENIC 自带方法进一步过滤基因

genesKept <- geneFiltering(
  exprMat_sub,
  scenicOptions
)

exprMat_filtered <- exprMat_sub[
  genesKept,
  ,
  drop = FALSE
]

if (nrow(exprMat_filtered) == 0) {
  stop(
    paste0(
      "SCENIC 自带 geneFiltering 后没有剩余基因，",
      "请调整预筛选参数。"
    )
  )
}


# 15. 计算基因相关性

runCorrelation(
  exprMat_filtered,
  scenicOptions
)


# 16. 使用 GENIE3 构建基因调控网络

runGenie3(
  exprMat_filtered,
  scenicOptions
)


# 17. 将共表达网络转换为调控模块

runSCENIC_1_coexNetwork2modules(
  scenicOptions
)


# 18. 基于转录因子结合基序创建 Regulon

runSCENIC_2_createRegulons(
  scenicOptions
)


# 19. 计算每个细胞的 Regulon 活性

runSCENIC_3_scoreCells(
  scenicOptions,
  exprMat_raw
)


# 20. 对 Regulon 活性进行二值化

runSCENIC_4_aucell_binarize(
  scenicOptions
)

message("SCENIC 分析完成，结果已输出到 SCENIC 文件夹。")