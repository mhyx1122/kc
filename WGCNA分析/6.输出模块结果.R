library(WGCNA)

# 1. 设置输入参数
nSelectgene <- 1000

# 2. 创建输出文件夹
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# 3. 导出所有基因所在模块
moduleColors <- dynamicColors
probes <- colnames(datExpr0)
geneInfo0 <- data.frame(
  probes = probes,
  moduleColor = moduleColors
)
geneOrder <- order(geneInfo0$moduleColor)
geneInfo <- geneInfo0[geneOrder,]
write.table(
  geneInfo,
  file = file.path(output_dir,"all_genes.csv"),
  sep = ",",
  row.names = FALSE,
  quote = FALSE
)

# 4. 输出每个模块的基因
module_gene_dir <- file.path(output_dir,"module_gene")
dir.create(module_gene_dir)
for (mod in 1:nrow(table(moduleColors))) {
  modules <- names(table(moduleColors))[mod]
  probes <- colnames(datExpr0)
  inModule <- moduleColors == modules
  modGenes <- probes[inModule]
  filePath <- file.path(module_gene_dir,paste0("module_",modules,".csv"))
  write.table(
    modGenes,
    file = filePath,
    sep = ",",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
}

# 5. 计算模块成员关系 MM
modNames <- substring(names(MEs),3)
geneModuleMembership <- as.data.frame(
  cor(datExpr0,MEs,use = "p")
)
MMPvalue <- as.data.frame(
  corPvalueStudent(
    as.matrix(geneModuleMembership),
    nSamples
  )
)
names(geneModuleMembership) <- paste("MM",modNames,sep = "")
names(MMPvalue) <- paste("p.MM",modNames,sep = "")

# 6. 计算基因性状显著性 GS
traitNames <- names(datTraits)
geneTraitSignificance <- as.data.frame(
  cor(datExpr0,datTraits,use = "p")
)
GSPvalue <- as.data.frame(
  corPvalueStudent(
    as.matrix(geneTraitSignificance),
    nSamples
  )
)
names(geneTraitSignificance) <- paste("GS.",traitNames,sep = "")
names(GSPvalue) <- paste("p.GS.",traitNames,sep = "")

# 7. 批量输出性状和模块散点图
picDir <- file.path(output_dir,"module_trait")
dir.create(picDir)
for (trait in traitNames) {
  traitColumn <- match(trait,traitNames)
  for (module in modNames) {
    column <- match(module,modNames)
    moduleGenes <- moduleColors == module
    if (nrow(geneModuleMembership[moduleGenes,]) > 1) {
      pdfFile <- paste("9_",trait,"_",module,".pdf",sep = "")
      outPdf <- file.path(picDir,pdfFile)
      pdf(file = outPdf,width = 7,height = 7)
      par(mfrow = c(1,1))
      verboseScatterplot(
        abs(geneModuleMembership[moduleGenes,column]),
        abs(geneTraitSignificance[moduleGenes,traitColumn]),
        xlab = paste("Module Membership in",module,"module"),
        ylab = paste("Gene significance for ",trait),
        main = paste("Module membership vs. gene significance\n"),
        cex.main = 1.2,
        cex.lab = 1.2,
        cex.axis = 1.2,
        col = module
      )
      dev.off()
    }
  }
}

# 8. 整理并输出 GS 和 MM 数据
probes <- colnames(datExpr0)
geneInfo0 <- data.frame(
  probes = probes,
  moduleColor = moduleColors
)
for (Tra in 1:ncol(geneTraitSignificance)) {
  oldNames <- names(geneInfo0)
  geneInfo0 <- data.frame(
    geneInfo0,
    geneTraitSignificance[,Tra],
    GSPvalue[,Tra]
  )
  names(geneInfo0) <- c(
    oldNames,
    names(geneTraitSignificance)[Tra],
    names(GSPvalue)[Tra]
  )
}
for (mod in 1:ncol(geneModuleMembership)) {
  oldNames <- names(geneInfo0)
  geneInfo0 <- data.frame(
    geneInfo0,
    geneModuleMembership[,mod],
    MMPvalue[,mod]
  )
  names(geneInfo0) <- c(
    oldNames,
    names(geneModuleMembership)[mod],
    names(MMPvalue)[mod]
  )
}
geneOrder <- order(geneInfo0$moduleColor)
geneInfo <- geneInfo0[geneOrder,]
write.table(
  geneInfo,
  file = file.path(output_dir,"GS_MM.xls"),
  sep = "\t",
  row.names = FALSE
)

# 9. 随机抽取部分基因并绘制 TOM 网络热图
nSelect <- nSelectgene
set.seed(10)
select <- sample(nGenes,size = nSelect)
selectTOM <- dissTOM[select,select]
selectTree <- hclust(
  as.dist(selectTOM),
  method = "average"
)
selectColors <- moduleColors[select]
plotDiss <- selectTOM^7
diag(plotDiss) <- NA
png(
  file = file.path(output_dir,"9_selectgenemap.png"),
  width = 9,
  height = 9,
  units = "in",
  res = 300
)
TOMplot(
  plotDiss,
  selectTree,
  selectColors,
  main = "Network heatmap plot, selected genes"
)
dev.off()

# 10. 绘制模块特征基因聚类树
pdf(
  file = file.path(output_dir,"10_dendrogram.pdf"),
  width = 8,
  height = 8
)
par(cex = 1.0)
plotEigengeneNetworks(
  MEs,
  "Eigengene dendrogram",
  marDendro = c(0,4,2,0),
  plotHeatmaps = FALSE
)
dev.off()

# 11. 绘制模块特征基因邻接热图
pdf(
  file = file.path(output_dir,"11_modulemap.pdf"),
  width = 8,
  height = 8
)
par(cex = 1.0)
plotEigengeneNetworks(
  MEs,
  "Eigengene adjacency heatmap",
  marHeatmap = c(3,4,2,2),
  plotDendrograms = FALSE,
  xLabelsAngle = 90
)
dev.off()

# 12. 绘制模块特征基因聚类树联合热图
pdf(
  file = file.path(output_dir,"12_dendrogrammap.pdf"),
  width = 9,
  height = 14
)
par(cex = 0.9)
plotEigengeneNetworks(
  MEs,
  "",
  marDendro = c(0,4,1,2),
  marHeatmap = c(3,4,1,2),
  cex.lab = 0.8,
  xLabelsAngle = 90
)
dev.off()

# 13. 创建 Cytoscape 输入文件夹
cytoDir <- file.path(output_dir,"CytoscapeInput")
dir.create(cytoDir)

# 14. 批量输出 Cytoscape 输入文件
num_modules <- nrow(table(moduleColors))
cat("Number of modules:",num_modules,"\n")
for (mod in 1:num_modules) {
  modules <- names(table(moduleColors))[mod]
  probes <- colnames(datExpr0)
  inModule <- moduleColors == modules
  modProbes <- probes[inModule]
  modGenes <- modProbes
  if (length(modProbes) == 0) {
    cat("modProbes is empty for module",modules,"\n")
    next
  }
  modTOM <- TOM[inModule,inModule]
  cat("Module:",modules,"Type of modTOM:",class(modTOM),"\n")
  cat("Module:",modules,"Dimensions of modTOM:",dim(modTOM),"\n")
  if (is.matrix(modTOM) || is.data.frame(modTOM)) {
    dimnames(modTOM) <- list(modProbes,modProbes)
  } else {
    cat("modTOM is not a matrix or data frame for module",modules,"\n")
    next
  }
  edges_File <- paste("CytoscapeInput-edges-",modules,".txt",sep = "")
  nodes_File <- paste("CytoscapeInput-nodes-",modules,".txt",sep = "")
  outEdge <- file.path(cytoDir,edges_File)
  outNode <- file.path(cytoDir,nodes_File)
  cyt <- exportNetworkToCytoscape(
    modTOM,
    edgeFile = outEdge,
    nodeFile = outNode,
    weighted = TRUE,
    threshold = 0.02,
    nodeNames = modProbes,
    altNodeNames = modGenes,
    nodeAttr = moduleColors[inModule]
  )
}

# 15. 输出完成提示
prompt <- "有bug请联系作者VX：cgxr410解决"
print(prompt)
print(prompt)
print(prompt)
message("成功绘制热图，已随机选择 ",nSelectgene," 个基因。")