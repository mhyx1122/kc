library(WGCNA)

# 1. 设置输入参数
minModuleSize <- 25
deepSplit <- 2
fig6_w <- 8
fig6_h <- 6
Pmethod <- "pearson"
width8 <- 10
height8 <- 8
ZLlim <- -1
ZHliM <- 1
par8 <- c(6,10,3,3)

# 2. 动态剪切基因模块
dynamicMods <- cutreeDynamic(
  dendro = geneTree,
  distM = dissTOM,
  deepSplit = deepSplit,
  pamRespectsDendro = FALSE,
  minClusterSize = minModuleSize
)
table(dynamicMods)

# 3. 将模块编号转换为模块颜色
dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)
print(table(dynamicColors))

# 4. 绘制并保存基因聚类树和模块颜色图
pdf(
  file = file.path(output_dir,"6_Tree.pdf"),
  width = fig6_w,
  height = fig6_h
)
plotDendroAndColors(
  geneTree,
  dynamicColors,
  "Dynamic Tree Cut",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Gene dendrogram and module colors"
)
dev.off()

# 5. 提取模块特征基因
MEList <- moduleEigengenes(
  datExpr0,
  colors = dynamicColors
)
MEs <- MEList$eigengenes

# 6. 计算模块特征基因之间的不相似度并聚类
MEDiss <- 1 - cor(MEs)
METree <- hclust(
  as.dist(MEDiss),
  method = "average"
)

# 7. 绘制并保存模块特征基因聚类图
pdf(
  file = file.path(output_dir,"7_module.pdf"),
  width = fig6_w,
  height = fig6_h
)
plot(
  METree,
  main = "Clustering of module eigengenes",
  xlab = "",
  sub = ""
)
dev.off()

# 8. 计算模块与性状之间的相关性和P值
nGenes <- ncol(datExpr0)
nSamples <- nrow(datExpr0)
moduleTraitCor <- cor(
  MEs,
  datTraits,
  use = "p",
  method = Pmethod
)
moduleTraitPvalue <- corPvalueStudent(
  moduleTraitCor,
  nSamples
)

# 9. 绘制并保存全部模块-性状相关性热图
pdf(
  file = file.path(output_dir,"8_Module.pdf"),
  width = width8,
  height = height8
)
textMatrix <- paste(
  signif(moduleTraitCor,2),
  "\n(",
  signif(moduleTraitPvalue,1),
  ")",
  sep = ""
)
dim(textMatrix) <- dim(moduleTraitCor)
par(mar = par8)
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = names(datTraits),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.5,
  zlim = c(ZLlim,ZHliM),
  main = paste("Module-trait relationships")
)
dev.off()

# 10. 判断每个模块是否至少有一个性状相关P值小于0.05
significantModuleRows <- apply(
  moduleTraitPvalue,
  1,
  function(p) any(p < 0.05)
)

# 11. 提取显著相关模块的相关性矩阵和P值矩阵
filteredModuleTraitCor <- moduleTraitCor[significantModuleRows,,drop = FALSE]
filteredModuleTraitPvalue <- moduleTraitPvalue[significantModuleRows,,drop = FALSE]

# 12. 生成显著模块热图文字矩阵
textMatrix <- paste(
  signif(filteredModuleTraitCor,2),
  "\n(",
  signif(filteredModuleTraitPvalue,1),
  ")",
  sep = ""
)
dim(textMatrix) <- dim(filteredModuleTraitCor)

# 13. 绘制并保存只保留显著模块的热图
pdf(
  file = file.path(output_dir,"8_Module(只保留显著).pdf"),
  width = width8,
  height = height8
)
par(mar = par8)
labeledHeatmap(
  Matrix = filteredModuleTraitCor,
  xLabels = names(datTraits),
  yLabels = names(MEs)[significantModuleRows],
  ySymbols = names(MEs)[significantModuleRows],
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.5,
  zlim = c(ZLlim,ZHliM),
  main = "Module-trait relationships"
)
dev.off()

# 14. 输出提示信息
onlyP <- "只显示显著相关（P＜0.05）的模块，当模块较多时可选择用此图"
print(onlyP)
print(onlyP)
print(onlyP)

# 15. 输出完成提示
message("分析完成！所有图像已保存。")