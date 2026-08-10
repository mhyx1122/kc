# 1. 加载必要的 R 包

suppressPackageStartupMessages({
  library(monocle3)
})


# 2. 检查 cds 对象

if (!exists("cds", envir = .GlobalEnv)) {
  stop("全局环境中没有找到 cds 对象，请先完成 Monocle3 轨迹学习分析。")
}


# 3. 设置拟时序起点选择方式

# FALSE：在轨迹图中手动选择起点
# TRUE：根据指定细胞类型自动确定起点
auto <- FALSE


# 4. 手动选择起点并计算拟时序

if (auto == FALSE) {
  
  cds <- order_cells(cds)
  
} else {
  
  # 5. 设置自动选择起点的参数
  
  # 用于确定起始细胞的注释列
  select_column <- "cellType"
  
  # 指定作为轨迹起点的细胞类型
  auto_select <- "TCells"
  
  
  # 6. 提取指定细胞类型对应的细胞
  
  cell_ids <- which(
    colData(cds)[, select_column] == auto_select
  )
  
  
  # 7. 获取每个细胞在轨迹主图中对应的最近节点
  
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$
    pr_graph_cell_proj_closest_vertex
  
  closest_vertex <- as.matrix(
    closest_vertex[colnames(cds), ]
  )
  
  
  # 8. 选择指定细胞类型中出现频率最高的轨迹节点作为起点
  
  root_pr_nodes <- igraph::V(
    principal_graph(cds)[["UMAP"]]
  )$name[
    as.numeric(
      names(
        which.max(
          table(
            closest_vertex[cell_ids, ]
          )
        )
      )
    )
  ]
  
  
  # 9. 使用自动确定的起始节点计算拟时序
  
  cds <- order_cells(
    cds,
    root_pr_nodes = root_pr_nodes
  )
}