suppressPackageStartupMessages({
  library(Seurat)
  library(patchwork)
  library(ggplot2)
  library(SeuratWrappers)
})

options(future.globals.maxSize = 500 * 1024^3)


# 1. 检查全局环境中是否存在 seurat 对象
if (!exists("seurat", envir = .GlobalEnv)) {
  stop("全局环境中未找到对象 seurat。请先生成或加载 seurat。")
}

srt0 <- get("seurat", envir = .GlobalEnv)


# 2. 输出目录模块

# 2.1 输出目录参数
scRNA_default_dir <- "1.3去批次对比图"

# 2.2 创建输出目录
if (!dir.exists(scRNA_default_dir)) {
  dir.create(scRNA_default_dir, recursive = TRUE)
}


# 3. Integration 参数模块

# 3.1 UMAP 使用的维度
scRNA_dims_n <- 20
dims_use <- 1:scRNA_dims_n

# 3.2 选择去批次方法
# 可选："CCA"、"RPCA"、"Harmony"、"FastMNN"、"JointPCA"
scRNA_methods <- c("Harmony", "FastMNN")

# 3.3 split 使用的 metadata 列
split_col <- "orig.ident"

# 3.4 绘图颜色
scRNA_palette_text <- "#E64B35FF,#4DBBD5FF,#00A087FF,#3C5488FF,#F39B7FFF,#8491B4FF,#91D1C2FF,#7E6148FF,#0073C2FF,#EFC000FF,#868686FF,#CD534CFF"

scRNA_palette <- unlist(strsplit(scRNA_palette_text, ",", fixed = TRUE))
scRNA_palette <- trimws(scRNA_palette)
scRNA_palette <- scRNA_palette[nchar(scRNA_palette) > 0]

if (length(scRNA_palette) == 0) {
  scRNA_palette <- NULL
}


# 4. 保存参数模块

# 4.1 拼图保存参数
scRNA_w_panel <- 14
scRNA_h_panel <- 8
scRNA_name_panel <- "integration_compare"


# 5. 参数检查模块

# 5.1 检查去批次方法
if (length(scRNA_methods) < 1) {
  stop("请至少选择一个去批次方法。")
}

# 5.2 检查 orig.ident
if (!split_col %in% colnames(srt0[[]])) {
  stop(paste0("metadata 中没有列：", split_col))
}

# 5.3 检查 PCA
if (!"pca" %in% Reductions(srt0)) {
  stop("seurat 对象中没有 reduction pca，请先 RunPCA。")
}


# 6. 去批次方法规格模块

method_spec <- list(
  CCA = list(
    method_fun = CCAIntegration,
    new_reduction = "integrated.cca",
    needs_orig = TRUE
  ),
  RPCA = list(
    method_fun = RPCAIntegration,
    new_reduction = "integrated.rpca",
    needs_orig = TRUE
  ),
  Harmony = list(
    method_fun = HarmonyIntegration,
    new_reduction = "harmony",
    needs_orig = TRUE
  ),
  FastMNN = list(
    method_fun = FastMNNIntegration,
    new_reduction = "integrated.mnn",
    needs_orig = FALSE
  ),
  JointPCA = list(
    method_fun = JointPCAIntegration,
    new_reduction = "integrated.Join",
    needs_orig = TRUE
  )
)

unknown_methods <- setdiff(scRNA_methods, names(method_spec))

if (length(unknown_methods) > 0) {
  stop(paste0(
    "存在未知去批次方法：",
    paste(unknown_methods, collapse = ", ")
  ))
}


# 7. split RNA assay 模块

srt <- srt0

srt[["RNA"]] <- split(
  srt[["RNA"]],
  f = srt[[split_col]][, 1]
)


# 8. 多方法 IntegrateLayers 模块

ok_methods <- character(0)
fail_methods <- character(0)
fail_msgs <- character(0)

for (i in seq_along(scRNA_methods)) {
  
  m <- scRNA_methods[[i]]
  spec <- method_spec[[m]]
  
  res <- tryCatch({
    
    if (isTRUE(spec$needs_orig)) {
      
      IntegrateLayers(
        object = srt,
        method = spec$method_fun,
        orig.reduction = "pca",
        new.reduction = spec$new_reduction,
        features = VariableFeatures(srt),
        verbose = FALSE
      )
      
    } else {
      
      IntegrateLayers(
        object = srt,
        method = spec$method_fun,
        new.reduction = spec$new_reduction,
        features = VariableFeatures(srt),
        verbose = FALSE
      )
    }
    
  }, error = function(e) {
    e
  })
  
  if (inherits(res, "error")) {
    
    fail_methods <- c(fail_methods, m)
    fail_msgs <- c(fail_msgs, conditionMessage(res))
    
  } else {
    
    srt <- res
    ok_methods <- c(ok_methods, m)
  }
}


# 9. JoinLayers 模块

srt <- JoinLayers(srt)


# 10. 计算 PCA 基线 UMAP 模块

srt <- RunUMAP(
  srt,
  dims = dims_use,
  reduction = "pca",
  verbose = FALSE
)


# 11. 计算各去批次方法对应的 UMAP 模块

for (m in ok_methods) {
  
  red <- method_spec[[m]]$new_reduction
  umap_tmp <- paste0("umap.", red, ".tmp")
  
  srt <- tryCatch({
    
    RunUMAP(
      srt,
      dims = dims_use,
      reduction = red,
      reduction.name = umap_tmp,
      verbose = FALSE
    )
    
  }, error = function(e) {
    srt
  })
}


# 12. 生成 PCA 基线图模块

p_base <- DimPlot(
  srt,
  reduction = "umap",
  group.by = "orig.ident"
) +
  ggtitle("UMAP (PCA baseline)")

if (!is.null(scRNA_palette)) {
  p_base <- p_base + scale_color_manual(values = scRNA_palette)
}


# 13. 生成各方法 UMAP 图模块

p_methods <- list()

for (m in ok_methods) {
  
  red <- method_spec[[m]]$new_reduction
  umap_tmp <- paste0("umap.", red, ".tmp")
  
  if (!umap_tmp %in% Reductions(srt)) {
    next
  }
  
  p <- DimPlot(
    srt,
    reduction = umap_tmp,
    group.by = "orig.ident"
  ) +
    ggtitle(paste0("UMAP (", red, ")"))
  
  if (!is.null(scRNA_palette)) {
    p <- p + scale_color_manual(values = scRNA_palette)
  }
  
  p_methods[[m]] <- p
}


# 14. 拼接对比图模块

all_plots <- c(
  list(p_base),
  unname(p_methods)
)

panel_ncol <- if (length(all_plots) >= 4) 3 else 2

p_panel <- wrap_plots(
  all_plots,
  ncol = panel_ncol
)


# 15. 保存对比图模块

ggsave(
  filename = file.path(scRNA_default_dir, paste0(scRNA_name_panel, ".pdf")),
  plot = p_panel,
  width = scRNA_w_panel,
  height = scRNA_h_panel,
  device = "pdf"
)


# 16. 保存运行信息模块

integration_info_text <- paste0(
  "本次去批次分析运行信息：\n",
  "- split_by：", split_col, "\n",
  "- selected_methods：", paste(scRNA_methods, collapse = ", "), "\n",
  "- ok_methods：", ifelse(length(ok_methods) > 0, paste(ok_methods, collapse = ", "), "无"), "\n",
  "- fail_methods：", ifelse(length(fail_methods) > 0, paste(fail_methods, collapse = ", "), "无"), "\n",
  "- fail_msgs：", ifelse(length(fail_msgs) > 0, paste(fail_msgs, collapse = " | "), "无"), "\n",
  "- reductions_in_seurat：", paste(Reductions(srt), collapse = ", "), "\n",
  "- umap_dims：1:", scRNA_dims_n, "\n",
  "- palette：", ifelse(is.null(scRNA_palette), "NULL", paste(scRNA_palette, collapse = ", ")), "\n"
)

writeLines(
  integration_info_text,
  con = file.path(scRNA_default_dir, "integration_run_info.txt")
)


# 17. 写回 seurat 对象

seurat <- srt