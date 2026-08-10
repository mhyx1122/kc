library(WGCNA)

# 1. 设置输入参数
Gree_width <- 12
Gree_hight <- 9

# 2. 构建邻接矩阵
adjacency <- adjacency(
  datExpr0,
  power = softPower
)

# 3. 计算拓扑重叠矩阵 TOM
TOM <- TOMsimilarity(adjacency)

# 4. 计算 TOM 不相似矩阵
dissTOM <- 1 - TOM

# 5. 还原 TOM 矩阵
TOM <- 1 - dissTOM

# 6. 基于 TOM 不相似矩阵构建基因聚类树
geneTree <- hclust(
  as.dist(dissTOM),
  method = method
)

# 7. 保存基因聚类树图
pdf(
  file = file.path(output_dir,"5_gene.pdf"),
  width = Gree_width,
  height = Gree_hight
)
plot(
  geneTree,
  xlab = "",
  sub = "",
  main = "Gene clustering on TOM-based dissimilarity",
  labels = FALSE,
  hang = 0.04
)
dev.off()

# 8. 输出完成提示
message("TOM矩阵构建完成，基因聚类树图已保存到：",file.path(output_dir,"5_gene.pdf"))