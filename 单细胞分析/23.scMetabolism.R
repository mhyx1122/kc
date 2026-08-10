library(scMetabolism)
library(ggplot2)
library(rsvd)
library(scCustomize)

metabolism_type = "KEGG"

# v5 转换到 V4
seurat_V4 <- Convert_Assay(seurat_object = seurat, convert_to = "V3", assay = "RNA")

countexp.Seurat <- suppressWarnings(sc.metabolism.Seurat(
  obj = seurat_V4,
  method = "AUCell",
  imputation = FALSE,
  ncores = 1,
  metabolism.type = metabolism_type
))

metabolism.matrix <- t(countexp.Seurat@assays$METABOLISM$score)
colnames(metabolism.matrix) <- make.names(colnames(metabolism.matrix), unique = TRUE)

if (
  identical(
    gsub("[.-]", "", rownames(seurat@meta.data)),
    gsub("[.-]", "", rownames(metabolism.matrix))
  )
) {
  message("两者行名顺序一致，可以合并，已合并到 meta.data 中。")
} else {
  warning("两者行名顺序不一致，不能合并，请检查结果。")
}

# 返回修改后的 Seurat 对象，而不是修改全局变量
seurat@meta.data[, colnames(metabolism.matrix)] <<- metabolism.matrix
message("可用的代谢通路名称：")
print(colnames(metabolism.matrix) )