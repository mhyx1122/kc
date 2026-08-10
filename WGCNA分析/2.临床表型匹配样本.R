library(WGCNA)

# 1. 设置输入参数
trait_file <- "表型数据.csv"

# 2. 读取表型数据
traitData <- read.csv(
  trait_file,
  row.names = 1,
  header = TRUE
)

# 3. 查看表型数据维度和列名
print(dim(traitData))
print(names(traitData))

# 4. 按表达矩阵样本顺序匹配表型数据
fpkmSamples <- rownames(datExpr0)
traitSamples <- rownames(traitData)
traitRows <- match(fpkmSamples,traitSamples)
datTraits <- traitData[traitRows,,drop = FALSE]
print(rownames(datTraits))

# 5. 对非离群样本重新聚类
sampleTree2 <- hclust(
  dist(datExpr0),
  method = method
)

# 6. 将表型数据转换为颜色矩阵
traitColors <- numbers2colors(
  datTraits,
  signed = FALSE
)

# 7. 绘制并保存样本聚类图和表型热图
pdf(
  file = file.path(output_dir,"3_sample_map.pdf"),
  width = tree_width,
  height = tree_hight
)
plotDendroAndColors(
  dendro = sampleTree2,
  colors = traitColors,
  groupLabels = names(datTraits),
  main = "Sample dendrogram and trait heatmap",
  ylab = "",
  cex.axis = cex.axis,
  cex.main = cex.main,
  cex.dendroLabels = cex,
  marAll = BianJu
)
dev.off()

# 8. 输出完成提示
message("数据处理完成。图像已保存为 PDF 文件：",file.path(output_dir,"3_sample_map.pdf"))