# 1. 加载必要 R 包

suppressPackageStartupMessages({
  library(CellChat)
  library(ggplot2)
})


# 2. 检查全局环境中的 cellchat 对象

if (!exists("cellchat", envir = .GlobalEnv)) {
  stop("全局环境中没有 cellchat 对象，请先运行 CellChat 分析。")
}

cellchat_obj <- get("cellchat", envir = .GlobalEnv)


# 3. 创建输出文件夹

# 3.1 总输出目录
out_dir <- "细胞通讯分析"

# 3.2 分模块输出目录
summary_dir <- file.path(out_dir, "细胞通讯汇总图")
single_dir <- file.path(out_dir, "单信号作图")
multi_dir <- file.path(out_dir, "多细胞多信号作图")

# 3.3 创建目录
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

if (!dir.exists(summary_dir)) {
  dir.create(summary_dir, recursive = TRUE)
}

if (!dir.exists(single_dir)) {
  dir.create(single_dir, recursive = TRUE)
}

if (!dir.exists(multi_dir)) {
  dir.create(multi_dir, recursive = TRUE)
}


# 4. 读取 CellChat 对象中的基础信息

# 4.1 读取信号通路
pathways <- cellchat_obj@netP$pathways

if (length(pathways) == 0) {
  stop("cellchat_obj@netP$pathways 为空，请确认已经运行 computeCommunProbPathway() 和 aggregateNet()。")
}

# 4.2 读取细胞类型
cell_types <- levels(cellchat_obj@idents)

if (length(cell_types) == 0) {
  stop("cellchat_obj@idents 中没有细胞类型信息。")
}

# 4.3 读取细胞群大小
group_size <- as.numeric(table(cellchat_obj@idents))

# 4.4 打印基础信息
cat("已读取 cellchat 对象。\n")
cat("细胞群数量：", length(cell_types), "\n")
cat("信号通路数量：", length(pathways), "\n")
cat("细胞类型：", paste(cell_types, collapse = ", "), "\n")
cat("信号通路：", paste(pathways, collapse = ", "), "\n")


# 5. 汇总作图板块：数据库分类图

# 5.1 数据库分类图参数
w_db_category <- 10
h_db_category <- 7
name_db_category <- "showDatabaseCategory"

# 5.2 绘制数据库分类图
p_db_category <- showDatabaseCategory(cellchat_obj@DB)

# 5.3 显示数据库分类图
print(p_db_category)

# 5.4 保存数据库分类图
file_db_category <- file.path(
  out_dir,
  paste0(name_db_category, ".pdf")
)

pdf(
  file = file_db_category,
  width = w_db_category,
  height = h_db_category
)

print(p_db_category)

dev.off()

cat("数据库分类图保存成功：", file_db_category, "\n")


# 6. 汇总作图板块：通讯网络图，按连接数量

# 6.1 通讯网络连接数量图参数
w_net_count <- 8
h_net_count <- 8
name_net_count <- "1.细胞-细胞通讯网络(连接数量)"

palette_count_text <- "#E64B35FF,#4DBBD5FF,#00A087FF,#3C5488FF,#F39B7FFF,#8491B4FF,#91D1C2FF,#7E6148FF,#0073C2FF,#EFC000FF,#868686FF,#CD534CFF,#7AA6DCFF,#003C67FF,#8F7700FF,#3B3B3BFF,#A73030FF,#374E55FF,#DF8F44FF,#00A1D5FF,#B24745FF,#79AF97FF,#6A6599FF,#80796BFF,#BC3C29FF,#0072B5FF,#E18727FF,#20854EFF,#7876B1FF,#6F99ADFF,#FFDC91FF,#EE4C97FF"

cluster_colors_count <- unlist(strsplit(palette_count_text, ","))
cluster_colors_count <- trimws(cluster_colors_count)
cluster_colors_count <- cluster_colors_count[cluster_colors_count != ""]

if (length(cluster_colors_count) < length(cell_types)) {
  stop("palette_count_text 中颜色数量不足以覆盖所有 cellType 分组，请补充颜色。")
}

cluster_colors_count <- setNames(
  cluster_colors_count[seq_along(cell_types)],
  cell_types
)

# 6.2 显示通讯网络连接数量图
par(mfrow = c(1, 1), xpd = TRUE)

netVisual_circle(
  cellchat_obj@net$count,
  vertex.weight = group_size,
  weight.scale = TRUE,
  color.use = cluster_colors_count,
  label.edge = TRUE,
  title.name = "Number of interactions"
)

# 6.3 保存通讯网络连接数量图
file_net_count <- file.path(
  summary_dir,
  paste0(name_net_count, ".pdf")
)

pdf(
  file = file_net_count,
  width = w_net_count,
  height = h_net_count
)

par(mfrow = c(1, 1), xpd = TRUE)

netVisual_circle(
  cellchat_obj@net$count,
  vertex.weight = group_size,
  weight.scale = TRUE,
  color.use = cluster_colors_count,
  label.edge = TRUE,
  title.name = "Number of interactions"
)

dev.off()

cat("通讯网络连接数量图保存成功：", file_net_count, "\n")


# 7. 汇总作图板块：通讯网络图，按连接强度

# 7.1 通讯网络连接强度图参数
w_net_weight <- 8
h_net_weight <- 8
name_net_weight <- "2.细胞-细胞通讯网络(连接强度)"

palette_weight_text <- "#E64B35FF,#4DBBD5FF,#00A087FF,#3C5488FF,#F39B7FFF,#8491B4FF,#91D1C2FF,#7E6148FF,#0073C2FF,#EFC000FF,#868686FF,#CD534CFF,#7AA6DCFF,#003C67FF,#8F7700FF,#3B3B3BFF,#A73030FF,#374E55FF,#DF8F44FF,#00A1D5FF,#B24745FF,#79AF97FF,#6A6599FF,#80796BFF,#BC3C29FF,#0072B5FF,#E18727FF,#20854EFF,#7876B1FF,#6F99ADFF,#FFDC91FF,#EE4C97FF"

cluster_colors_weight <- unlist(strsplit(palette_weight_text, ","))
cluster_colors_weight <- trimws(cluster_colors_weight)
cluster_colors_weight <- cluster_colors_weight[cluster_colors_weight != ""]

if (length(cluster_colors_weight) < length(cell_types)) {
  stop("palette_weight_text 中颜色数量不足以覆盖所有 cellType 分组，请补充颜色。")
}

cluster_colors_weight <- setNames(
  cluster_colors_weight[seq_along(cell_types)],
  cell_types
)

# 7.2 显示通讯网络连接强度图
par(mfrow = c(1, 1), xpd = TRUE)

netVisual_circle(
  cellchat_obj@net$weight,
  vertex.weight = group_size,
  weight.scale = TRUE,
  color.use = cluster_colors_weight,
  label.edge = TRUE,
  title.name = "Interaction weights/strength"
)

# 7.3 保存通讯网络连接强度图
file_net_weight <- file.path(
  summary_dir,
  paste0(name_net_weight, ".pdf")
)

pdf(
  file = file_net_weight,
  width = w_net_weight,
  height = h_net_weight
)

par(mfrow = c(1, 1), xpd = TRUE)

netVisual_circle(
  cellchat_obj@net$weight,
  vertex.weight = group_size,
  weight.scale = TRUE,
  color.use = cluster_colors_weight,
  label.edge = TRUE,
  title.name = "Interaction weights/strength"
)

dev.off()

cat("通讯网络连接强度图保存成功：", file_net_weight, "\n")


# 8. 汇总作图板块：单个细胞群发送的信号

# 8.1 单个细胞群发送信号图参数
w_sender_signal <- 12
h_sender_signal <- 10
name_sender_signal <- "3.单个细胞群发送的信号"

signal_rows <- 3
signal_cols <- 4

palette_sender_text <- "#E64B35FF,#4DBBD5FF,#00A087FF,#3C5488FF,#F39B7FFF,#8491B4FF,#91D1C2FF,#7E6148FF,#0073C2FF,#EFC000FF,#868686FF,#CD534CFF,#7AA6DCFF,#003C67FF,#8F7700FF,#3B3B3BFF,#A73030FF,#374E55FF,#DF8F44FF,#00A1D5FF,#B24745FF,#79AF97FF,#6A6599FF,#80796BFF,#BC3C29FF,#0072B5FF,#E18727FF,#20854EFF,#7876B1FF,#6F99ADFF,#FFDC91FF,#EE4C97FF"

cluster_colors_sender <- unlist(strsplit(palette_sender_text, ","))
cluster_colors_sender <- trimws(cluster_colors_sender)
cluster_colors_sender <- cluster_colors_sender[cluster_colors_sender != ""]

if (length(cluster_colors_sender) < length(cell_types)) {
  stop("palette_sender_text 中颜色数量不足以覆盖所有 cellType 分组，请补充颜色。")
}

cluster_colors_sender <- setNames(
  cluster_colors_sender[seq_along(cell_types)],
  cell_types
)

signal_matrix <- cellchat_obj@net$weight

# 8.2 显示单个细胞群发送信号图
par(mfrow = c(signal_rows, signal_cols), xpd = TRUE)

for (i in seq_len(nrow(signal_matrix))) {
  temp_matrix <- matrix(
    0,
    nrow = nrow(signal_matrix),
    ncol = ncol(signal_matrix),
    dimnames = dimnames(signal_matrix)
  )
  
  temp_matrix[i, ] <- signal_matrix[i, ]
  
  netVisual_circle(
    temp_matrix,
    vertex.weight = group_size,
    weight.scale = TRUE,
    color.use = cluster_colors_sender,
    label.edge = TRUE,
    edge.weight.max = max(signal_matrix),
    title.name = rownames(signal_matrix)[i]
  )
}

# 8.3 保存单个细胞群发送信号图
file_sender_signal <- file.path(
  summary_dir,
  paste0(name_sender_signal, ".pdf")
)

pdf(
  file = file_sender_signal,
  width = w_sender_signal,
  height = h_sender_signal
)

par(mfrow = c(signal_rows, signal_cols), xpd = TRUE)

for (i in seq_len(nrow(signal_matrix))) {
  temp_matrix <- matrix(
    0,
    nrow = nrow(signal_matrix),
    ncol = ncol(signal_matrix),
    dimnames = dimnames(signal_matrix)
  )
  
  temp_matrix[i, ] <- signal_matrix[i, ]
  
  netVisual_circle(
    temp_matrix,
    vertex.weight = group_size,
    weight.scale = TRUE,
    color.use = cluster_colors_sender,
    label.edge = TRUE,
    edge.weight.max = max(signal_matrix),
    title.name = rownames(signal_matrix)[i]
  )
}

dev.off()

cat("单个细胞群发送信号图保存成功：", file_sender_signal, "\n")


# 9. 单信号作图板块：设置单信号参数

# 9.1 单信号基础参数
single_signal <- pathways[1]

single_target_cells <- cell_types[seq_len(min(2, length(cell_types)))]

single_heatmap_color <- "Reds"
single_font_size <- 10
single_role_network_width <- 8
single_role_network_height <- 2.5

palette_single_text <- "#E64B35FF,#4DBBD5FF,#00A087FF,#3C5488FF,#F39B7FFF,#8491B4FF,#91D1C2FF,#7E6148FF,#0073C2FF,#EFC000FF,#868686FF,#CD534CFF,#7AA6DCFF,#003C67FF,#8F7700FF,#3B3B3BFF,#A73030FF,#374E55FF,#DF8F44FF,#00A1D5FF,#B24745FF,#79AF97FF,#6A6599FF,#80796BFF,#BC3C29FF,#0072B5FF,#E18727FF,#20854EFF,#7876B1FF,#6F99ADFF,#FFDC91FF,#EE4C97FF"

cluster_colors_single <- unlist(strsplit(palette_single_text, ","))
cluster_colors_single <- trimws(cluster_colors_single)
cluster_colors_single <- cluster_colors_single[cluster_colors_single != ""]

if (length(cluster_colors_single) < length(cell_types)) {
  stop("palette_single_text 中颜色数量不足以覆盖所有 cellType 分组，请补充颜色。")
}

cluster_colors_single <- setNames(
  cluster_colors_single[seq_along(cell_types)],
  cell_types
)

single_target_cell_indices <- match(single_target_cells, cell_types)
single_target_cell_indices <- single_target_cell_indices[!is.na(single_target_cell_indices)]

if (length(single_target_cell_indices) == 0) {
  stop("single_target_cells 没有匹配到任何细胞类型，请检查参数。")
}

cat("当前单信号通路：", single_signal, "\n")
cat("当前单信号层级图关注细胞：", paste(single_target_cells, collapse = ", "), "\n")


# 10. 单信号作图板块：特定信号圈图

# 10.1 特定信号圈图参数
w_single_circle <- 15
h_single_circle <- 8
name_single_circle <- "3.细胞通讯_特定信号"

# 10.2 显示特定信号圈图
par(mfrow = c(1, 1), xpd = TRUE)

netVisual_aggregate(
  cellchat_obj,
  color.use = cluster_colors_single,
  signaling = single_signal,
  layout = "circle"
)

# 10.3 保存特定信号圈图
file_single_circle <- file.path(
  single_dir,
  paste0(name_single_circle, ".pdf")
)

pdf(
  file = file_single_circle,
  width = w_single_circle,
  height = h_single_circle
)

par(mfrow = c(1, 1), xpd = TRUE)

netVisual_aggregate(
  cellchat_obj,
  color.use = cluster_colors_single,
  signaling = single_signal,
  layout = "circle"
)

dev.off()

cat("特定信号圈图保存成功：", file_single_circle, "\n")


# 11. 单信号作图板块：特定信号弦图

# 11.1 特定信号弦图参数
w_single_chord <- 15
h_single_chord <- 8
name_single_chord <- "4.细胞通讯_特定信号弦图"

# 11.2 显示特定信号弦图
par(mfrow = c(1, 1), xpd = TRUE)

netVisual_aggregate(
  cellchat_obj,
  color.use = cluster_colors_single,
  signaling = single_signal,
  layout = "chord"
)

if (requireNamespace("circlize", quietly = TRUE)) {
  circlize::circos.clear()
}

# 11.3 保存特定信号弦图
file_single_chord <- file.path(
  single_dir,
  paste0(name_single_chord, ".pdf")
)

pdf(
  file = file_single_chord,
  width = w_single_chord,
  height = h_single_chord
)

par(mfrow = c(1, 1), xpd = TRUE)

netVisual_aggregate(
  cellchat_obj,
  color.use = cluster_colors_single,
  signaling = single_signal,
  layout = "chord"
)

if (requireNamespace("circlize", quietly = TRUE)) {
  circlize::circos.clear()
}

dev.off()

cat("特定信号弦图保存成功：", file_single_chord, "\n")


# 12. 单信号作图板块：特定信号热图

# 12.1 特定信号热图参数
w_single_heatmap <- 15
h_single_heatmap <- 8
name_single_heatmap <- "5.细胞通讯_特定信号热图"

# 12.2 绘制特定信号热图
p_single_heatmap <- netVisual_heatmap(
  cellchat_obj,
  color.use = cluster_colors_single,
  signaling = single_signal,
  color.heatmap = single_heatmap_color
)

# 12.3 显示特定信号热图
print(p_single_heatmap)

# 12.4 保存特定信号热图
file_single_heatmap <- file.path(
  single_dir,
  paste0(name_single_heatmap, ".pdf")
)

pdf(
  file = file_single_heatmap,
  width = w_single_heatmap,
  height = h_single_heatmap
)

print(p_single_heatmap)

dev.off()

cat("特定信号热图保存成功：", file_single_heatmap, "\n")


# 13. 单信号作图板块：单信号层级图

# 13.1 单信号层级图参数
w_hierarchy_single <- 25
h_hierarchy_single <- 8
name_hierarchy_single <- "9.细胞通讯_单信号下的层级图"

# 13.2 显示单信号层级图
par(mfrow = c(1, 1), xpd = TRUE)

netVisual_aggregate(
  cellchat_obj,
  color.use = cluster_colors_single,
  vertex.receiver = single_target_cell_indices,
  signaling = single_signal,
  layout = "hierarchy"
)

# 13.3 保存单信号层级图
file_hierarchy_single <- file.path(
  single_dir,
  paste0(name_hierarchy_single, ".pdf")
)

pdf(
  file = file_hierarchy_single,
  width = w_hierarchy_single,
  height = h_hierarchy_single
)

par(mfrow = c(1, 1), xpd = TRUE)

netVisual_aggregate(
  cellchat_obj,
  color.use = cluster_colors_single,
  vertex.receiver = single_target_cell_indices,
  signaling = single_signal,
  layout = "hierarchy"
)

dev.off()

cat("单信号层级图保存成功：", file_hierarchy_single, "\n")


# 14. 单信号作图板块：信号角色网络图

# 14.1 信号角色网络图参数
w_single_role_network <- 8
h_single_role_network <- 8
name_single_role_network <- "6.1细胞群在该信号中的角色"

# 14.2 计算中心性
cellchat_single_centrality <- netAnalysis_computeCentrality(
  cellchat_obj,
  slot.name = "netP"
)

# 14.3 显示信号角色网络图
par(mfrow = c(1, 1), xpd = TRUE)

netAnalysis_signalingRole_network(
  cellchat_single_centrality,
  color.use = cluster_colors_single,
  signaling = single_signal,
  width = single_role_network_width,
  height = single_role_network_height,
  font.size = single_font_size
)

# 14.4 保存信号角色网络图
file_single_role_network <- file.path(
  single_dir,
  paste0(name_single_role_network, ".pdf")
)

pdf(
  file = file_single_role_network,
  width = w_single_role_network,
  height = h_single_role_network
)

par(mfrow = c(1, 1), xpd = TRUE)

netAnalysis_signalingRole_network(
  cellchat_single_centrality,
  color.use = cluster_colors_single,
  signaling = single_signal,
  width = single_role_network_width,
  height = single_role_network_height,
  font.size = single_font_size
)

dev.off()

cat("信号角色网络图保存成功：", file_single_role_network, "\n")


# 15. 单信号作图板块：信号角色散点图

# 15.1 信号角色散点图参数
w_single_role_scatter <- 6
h_single_role_scatter <- 5
name_single_role_scatter <- "6.2网络中信号角色在指定通路中的权重"

# 15.2 绘制信号角色散点图
p_single_role_scatter <- netAnalysis_signalingRole_scatter(
  cellchat_single_centrality,
  color.use = cluster_colors_single,
  signaling = single_signal
)

# 15.3 显示信号角色散点图
print(p_single_role_scatter)

# 15.4 保存信号角色散点图
file_single_role_scatter <- file.path(
  single_dir,
  paste0(name_single_role_scatter, ".pdf")
)

pdf(
  file = file_single_role_scatter,
  width = w_single_role_scatter,
  height = h_single_role_scatter
)

print(p_single_role_scatter)

dev.off()

cat("信号角色散点图保存成功：", file_single_role_scatter, "\n")


# 16. 单信号作图板块：提取单信号配体-受体对

# 16.1 提取单信号配体-受体对
single_lr_df <- extractEnrichedLR(
  cellchat_obj,
  signaling = single_signal,
  geneLR.return = FALSE
)

if (is.null(single_lr_df) || nrow(single_lr_df) == 0) {
  stop("当前 single_signal 下没有提取到配体-受体对。")
}

# 16.2 整理配体-受体标签
if ("interaction_name_2" %in% colnames(single_lr_df)) {
  single_lr_labels <- as.character(single_lr_df$interaction_name_2)
} else if ("interaction_name" %in% colnames(single_lr_df)) {
  single_lr_labels <- as.character(single_lr_df$interaction_name)
} else {
  single_lr_labels <- apply(
    single_lr_df,
    1,
    function(x) paste(x, collapse = " | ")
  )
}

single_lr_labels[is.na(single_lr_labels) | single_lr_labels == ""] <- paste0(
  "LR_",
  seq_along(single_lr_labels)
)[is.na(single_lr_labels) | single_lr_labels == ""]

# 16.3 选择单个配体-受体对
# 如需指定某个配体-受体对，可手动修改 single_lr
single_lr <- single_lr_labels[1]

single_lr_index <- match(single_lr, single_lr_labels)

if (is.na(single_lr_index)) {
  single_lr_index <- 1
}

single_lr_selected <- single_lr_df[single_lr_index, , drop = FALSE]

cat("当前单信号配体-受体对：", single_lr, "\n")


# 17. 单信号作图板块：单信号 LR 圈图

# 17.1 单信号 LR 圈图参数
w_single_lr_circle <- 15
h_single_lr_circle <- 8
name_single_lr_circle <- "10.1单信号下_配受体对水平的细胞通讯_圈图"

# 17.2 显示单信号 LR 圈图
par(mfrow = c(1, 1), xpd = TRUE)

netVisual_individual(
  cellchat_obj,
  color.use = cluster_colors_single,
  signaling = single_signal,
  pairLR.use = single_lr_selected,
  layout = "circle"
)

# 17.3 保存单信号 LR 圈图
file_single_lr_circle <- file.path(
  single_dir,
  paste0(name_single_lr_circle, ".pdf")
)

pdf(
  file = file_single_lr_circle,
  width = w_single_lr_circle,
  height = h_single_lr_circle
)

par(mfrow = c(1, 1), xpd = TRUE)

netVisual_individual(
  cellchat_obj,
  color.use = cluster_colors_single,
  signaling = single_signal,
  pairLR.use = single_lr_selected,
  layout = "circle"
)

dev.off()

cat("单信号 LR 圈图保存成功：", file_single_lr_circle, "\n")


# 18. 单信号作图板块：单信号 LR 弦图

# 18.1 单信号 LR 弦图参数
w_single_lr_chord <- 15
h_single_lr_chord <- 8
name_single_lr_chord <- "10.2单信号下_配受体对水平的细胞通讯_弦图"

# 18.2 显示单信号 LR 弦图
par(mfrow = c(1, 1), xpd = TRUE)

netVisual_individual(
  cellchat_obj,
  color.use = cluster_colors_single,
  signaling = single_signal,
  pairLR.use = single_lr_selected,
  layout = "chord"
)

if (requireNamespace("circlize", quietly = TRUE)) {
  circlize::circos.clear()
}

# 18.3 保存单信号 LR 弦图
file_single_lr_chord <- file.path(
  single_dir,
  paste0(name_single_lr_chord, ".pdf")
)

pdf(
  file = file_single_lr_chord,
  width = w_single_lr_chord,
  height = h_single_lr_chord
)

par(mfrow = c(1, 1), xpd = TRUE)

netVisual_individual(
  cellchat_obj,
  color.use = cluster_colors_single,
  signaling = single_signal,
  pairLR.use = single_lr_selected,
  layout = "chord"
)

if (requireNamespace("circlize", quietly = TRUE)) {
  circlize::circos.clear()
}

dev.off()

cat("单信号 LR 弦图保存成功：", file_single_lr_chord, "\n")


# 19. 单信号作图板块：单信号 LR 层级图

# 19.1 单信号 LR 层级图参数
w_single_lr_hierarchy <- 15
h_single_lr_hierarchy <- 8
name_single_lr_hierarchy <- "10.3单信号下_配受体对水平的细胞通讯_层级图"

# 19.2 显示单信号 LR 层级图
par(mfrow = c(1, 1), xpd = TRUE)

netVisual_individual(
  cellchat_obj,
  color.use = cluster_colors_single,
  signaling = single_signal,
  pairLR.use = single_lr_selected,
  layout = "hierarchy",
  vertex.receiver = single_target_cell_indices
)

# 19.3 保存单信号 LR 层级图
file_single_lr_hierarchy <- file.path(
  single_dir,
  paste0(name_single_lr_hierarchy, ".pdf")
)

pdf(
  file = file_single_lr_hierarchy,
  width = w_single_lr_hierarchy,
  height = h_single_lr_hierarchy
)

par(mfrow = c(1, 1), xpd = TRUE)

netVisual_individual(
  cellchat_obj,
  color.use = cluster_colors_single,
  signaling = single_signal,
  pairLR.use = single_lr_selected,
  layout = "hierarchy",
  vertex.receiver = single_target_cell_indices
)

dev.off()

cat("单信号 LR 层级图保存成功：", file_single_lr_hierarchy, "\n")


# 20. 单信号作图板块：配体受体贡献排名

# 20.1 配体受体贡献排名图参数
w_single_contribution <- 8
h_single_contribution <- 8
name_single_contribution <- "7.单信号：信号通路-配体受体贡献排名"

# 20.2 绘制配体受体贡献排名图
p_single_contribution <- netAnalysis_contribution(
  cellchat_obj,
  signaling = single_signal
)

# 20.3 显示配体受体贡献排名图
print(p_single_contribution)

# 20.4 保存配体受体贡献排名图
file_single_contribution <- file.path(
  single_dir,
  paste0(name_single_contribution, ".pdf")
)

pdf(
  file = file_single_contribution,
  width = w_single_contribution,
  height = h_single_contribution
)

print(p_single_contribution)

dev.off()

cat("单信号配体受体贡献排名图保存成功：", file_single_contribution, "\n")


# 21. 多细胞多信号作图板块：设置多信号参数

# 21.1 多信号基础参数
multi_signal <- pathways[seq_len(min(3, length(pathways)))]

target_cells <- cell_types[seq_len(min(2, length(cell_types)))]
receiver_cells <- cell_types[seq_len(min(2, length(cell_types)))]
signal_source <- cell_types[1]

palette_multi_text <- "#E64B35FF,#4DBBD5FF,#00A087FF,#3C5488FF,#F39B7FFF,#8491B4FF,#91D1C2FF,#7E6148FF,#0073C2FF,#EFC000FF,#868686FF,#CD534CFF,#7AA6DCFF,#003C67FF,#8F7700FF,#3B3B3BFF,#A73030FF,#374E55FF,#DF8F44FF,#00A1D5FF,#B24745FF,#79AF97FF,#6A6599FF,#80796BFF,#BC3C29FF,#0072B5FF,#E18727FF,#20854EFF,#7876B1FF,#6F99ADFF,#FFDC91FF,#EE4C97FF"

cluster_colors_multi <- unlist(strsplit(palette_multi_text, ","))
cluster_colors_multi <- trimws(cluster_colors_multi)
cluster_colors_multi <- cluster_colors_multi[cluster_colors_multi != ""]

if (length(cluster_colors_multi) < length(cell_types)) {
  stop("palette_multi_text 中颜色数量不足以覆盖所有 cellType 分组，请补充颜色。")
}

cluster_colors_multi <- setNames(
  cluster_colors_multi[seq_along(cell_types)],
  cell_types
)

target_cell_indices <- match(target_cells, cell_types)
target_cell_indices <- target_cell_indices[!is.na(target_cell_indices)]

if (length(target_cell_indices) == 0) {
  stop("target_cells 没有匹配到任何细胞类型，请检查参数。")
}

if (length(multi_signal) == 0) {
  stop("multi_signal 为空，请至少选择一个信号通路。")
}

cat("当前多信号通路：", paste(multi_signal, collapse = ", "), "\n")
cat("当前多信号层级图关注细胞：", paste(target_cells, collapse = ", "), "\n")
cat("当前气泡图发送细胞：", paste(signal_source, collapse = ", "), "\n")
cat("当前气泡图接收细胞：", paste(receiver_cells, collapse = ", "), "\n")


# 22. 多细胞多信号作图板块：多信号圈图

# 22.1 多信号圈图参数
w_multi_circle <- 15
h_multi_circle <- 8
name_multi_circle <- "8.1多信号_圈图"

# 22.2 显示多信号圈图
par(mfrow = c(1, 1), xpd = TRUE)

netVisual_aggregate(
  cellchat_obj,
  color.use = cluster_colors_multi,
  signaling = multi_signal,
  layout = "circle"
)

# 22.3 保存多信号圈图
file_multi_circle <- file.path(
  multi_dir,
  paste0(name_multi_circle, ".pdf")
)

pdf(
  file = file_multi_circle,
  width = w_multi_circle,
  height = h_multi_circle
)

par(mfrow = c(1, 1), xpd = TRUE)

netVisual_aggregate(
  cellchat_obj,
  color.use = cluster_colors_multi,
  signaling = multi_signal,
  layout = "circle"
)

dev.off()

cat("多信号圈图保存成功：", file_multi_circle, "\n")


# 23. 多细胞多信号作图板块：多信号弦图

# 23.1 多信号弦图参数
w_multi_chord <- 15
h_multi_chord <- 8
name_multi_chord <- "8.2多信号_弦图"

# 23.2 显示多信号弦图
par(mfrow = c(1, 1), xpd = TRUE)

netVisual_aggregate(
  cellchat_obj,
  color.use = cluster_colors_multi,
  signaling = multi_signal,
  layout = "chord"
)

if (requireNamespace("circlize", quietly = TRUE)) {
  circlize::circos.clear()
}

# 23.3 保存多信号弦图
file_multi_chord <- file.path(
  multi_dir,
  paste0(name_multi_chord, ".pdf")
)

pdf(
  file = file_multi_chord,
  width = w_multi_chord,
  height = h_multi_chord
)

par(mfrow = c(1, 1), xpd = TRUE)

netVisual_aggregate(
  cellchat_obj,
  color.use = cluster_colors_multi,
  signaling = multi_signal,
  layout = "chord"
)

if (requireNamespace("circlize", quietly = TRUE)) {
  circlize::circos.clear()
}

dev.off()

cat("多信号弦图保存成功：", file_multi_chord, "\n")


# 24. 多细胞多信号作图板块：多信号贡献排名

# 24.1 多信号贡献排名图参数
w_multi_contribution <- 8
h_multi_contribution <- 8
name_multi_contribution <- "7.多信号：信号通路-配体受体贡献排名"

# 24.2 绘制多信号贡献排名图
p_multi_contribution <- netAnalysis_contribution(
  cellchat_obj,
  signaling = multi_signal
)

# 24.3 显示多信号贡献排名图
print(p_multi_contribution)

# 24.4 保存多信号贡献排名图
file_multi_contribution <- file.path(
  multi_dir,
  paste0(name_multi_contribution, ".pdf")
)

pdf(
  file = file_multi_contribution,
  width = w_multi_contribution,
  height = h_multi_contribution
)

print(p_multi_contribution)

dev.off()

cat("多信号贡献排名图保存成功：", file_multi_contribution, "\n")


# 25. 多细胞多信号作图板块：多信号基因表达图

# 25.1 多信号基因表达图参数
w_multi_gene_expr <- 8
h_multi_gene_expr <- 15
name_multi_gene_expr <- "8.多信号：信号通路中基因在细胞亚群中的表达量"

# 25.2 绘制多信号基因表达图
p_multi_gene_expr <- plotGeneExpression(
  cellchat_obj,
  color.use = cluster_colors_multi,
  signaling = multi_signal
)

# 25.3 显示多信号基因表达图
suppressWarnings(print(p_multi_gene_expr))

# 25.4 保存多信号基因表达图
file_multi_gene_expr <- file.path(
  multi_dir,
  paste0(name_multi_gene_expr, ".pdf")
)

pdf(
  file = file_multi_gene_expr,
  width = w_multi_gene_expr,
  height = h_multi_gene_expr
)

suppressWarnings(print(p_multi_gene_expr))

dev.off()

cat("多信号基因表达图保存成功：", file_multi_gene_expr, "\n")


# 26. 多细胞多信号作图板块：层级图，多信号

# 26.1 层级图参数
w_hierarchy_multi <- 25
h_hierarchy_multi <- 8
name_hierarchy_multi <- "9.细胞通讯_多信号下的层级图"

# 26.2 显示多信号层级图
par(mfrow = c(1, 1), xpd = TRUE)

netVisual_aggregate(
  cellchat_obj,
  color.use = cluster_colors_multi,
  vertex.receiver = target_cell_indices,
  signaling = multi_signal,
  layout = "hierarchy"
)

# 26.3 保存多信号层级图
file_hierarchy_multi <- file.path(
  multi_dir,
  paste0(name_hierarchy_multi, ".pdf")
)

pdf(
  file = file_hierarchy_multi,
  width = w_hierarchy_multi,
  height = h_hierarchy_multi
)

par(mfrow = c(1, 1), xpd = TRUE)

netVisual_aggregate(
  cellchat_obj,
  color.use = cluster_colors_multi,
  vertex.receiver = target_cell_indices,
  signaling = multi_signal,
  layout = "hierarchy"
)

dev.off()

cat("多信号层级图保存成功：", file_hierarchy_multi, "\n")


# 27. 多细胞多信号作图板块：提取多信号配体-受体对

# 27.1 提取多信号配体-受体对
multi_lr_df <- extractEnrichedLR(
  cellchat_obj,
  signaling = multi_signal,
  geneLR.return = FALSE
)

if (is.null(multi_lr_df) || nrow(multi_lr_df) == 0) {
  stop("当前 multi_signal 下没有提取到配体-受体对。")
}

# 27.2 整理多信号配体-受体标签
if ("interaction_name_2" %in% colnames(multi_lr_df)) {
  multi_lr_labels <- as.character(multi_lr_df$interaction_name_2)
} else if ("interaction_name" %in% colnames(multi_lr_df)) {
  multi_lr_labels <- as.character(multi_lr_df$interaction_name)
} else {
  multi_lr_labels <- apply(
    multi_lr_df,
    1,
    function(x) paste(x, collapse = " | ")
  )
}

multi_lr_labels[is.na(multi_lr_labels) | multi_lr_labels == ""] <- paste0(
  "LR_",
  seq_along(multi_lr_labels)
)[is.na(multi_lr_labels) | multi_lr_labels == ""]

# 27.3 选择多信号下的单个配体-受体对
# 如需指定某个配体-受体对，可手动修改 multi_lr
multi_lr <- multi_lr_labels[1]

multi_lr_index <- match(multi_lr, multi_lr_labels)

if (is.na(multi_lr_index)) {
  multi_lr_index <- 1
}

multi_lr_selected <- multi_lr_df[multi_lr_index, , drop = FALSE]

cat("当前多信号配体-受体对：", multi_lr, "\n")


# 28. 多细胞多信号作图板块：多信号 LR 圈图

# 28.1 多信号 LR 圈图参数
w_multi_lr_circle <- 15
h_multi_lr_circle <- 8
name_multi_lr_circle <- "10.1多信号下_配受体对水平的细胞通讯_圈图"

# 28.2 显示多信号 LR 圈图
par(mfrow = c(1, 1), xpd = TRUE)

netVisual_individual(
  cellchat_obj,
  color.use = cluster_colors_multi,
  signaling = multi_signal,
  pairLR.use = multi_lr_selected,
  layout = "circle"
)

# 28.3 保存多信号 LR 圈图
file_multi_lr_circle <- file.path(
  multi_dir,
  paste0(name_multi_lr_circle, ".pdf")
)

pdf(
  file = file_multi_lr_circle,
  width = w_multi_lr_circle,
  height = h_multi_lr_circle
)

par(mfrow = c(1, 1), xpd = TRUE)

netVisual_individual(
  cellchat_obj,
  color.use = cluster_colors_multi,
  signaling = multi_signal,
  pairLR.use = multi_lr_selected,
  layout = "circle"
)

dev.off()

cat("多信号 LR 圈图保存成功：", file_multi_lr_circle, "\n")


# 29. 多细胞多信号作图板块：多信号 LR 弦图

# 29.1 多信号 LR 弦图参数
w_multi_lr_chord <- 15
h_multi_lr_chord <- 8
name_multi_lr_chord <- "10.2多信号下_配受体对水平的细胞通讯_弦图"

# 29.2 显示多信号 LR 弦图
par(mfrow = c(1, 1), xpd = TRUE)

netVisual_individual(
  cellchat_obj,
  color.use = cluster_colors_multi,
  signaling = multi_signal,
  pairLR.use = multi_lr_selected,
  layout = "chord"
)

if (requireNamespace("circlize", quietly = TRUE)) {
  circlize::circos.clear()
}

# 29.3 保存多信号 LR 弦图
file_multi_lr_chord <- file.path(
  multi_dir,
  paste0(name_multi_lr_chord, ".pdf")
)

pdf(
  file = file_multi_lr_chord,
  width = w_multi_lr_chord,
  height = h_multi_lr_chord
)

par(mfrow = c(1, 1), xpd = TRUE)

netVisual_individual(
  cellchat_obj,
  color.use = cluster_colors_multi,
  signaling = multi_signal,
  pairLR.use = multi_lr_selected,
  layout = "chord"
)

if (requireNamespace("circlize", quietly = TRUE)) {
  circlize::circos.clear()
}

dev.off()

cat("多信号 LR 弦图保存成功：", file_multi_lr_chord, "\n")


# 30. 多细胞多信号作图板块：多信号 LR 层级图

# 30.1 多信号 LR 层级图参数
w_multi_lr_hierarchy <- 15
h_multi_lr_hierarchy <- 8
name_multi_lr_hierarchy <- "10.3多信号下_配受体对水平的细胞通讯_层级图"

# 30.2 显示多信号 LR 层级图
par(mfrow = c(1, 1), xpd = TRUE)

netVisual_individual(
  cellchat_obj,
  color.use = cluster_colors_multi,
  signaling = multi_signal,
  pairLR.use = multi_lr_selected,
  layout = "hierarchy",
  vertex.receiver = target_cell_indices
)

# 30.3 保存多信号 LR 层级图
file_multi_lr_hierarchy <- file.path(
  multi_dir,
  paste0(name_multi_lr_hierarchy, ".pdf")
)

pdf(
  file = file_multi_lr_hierarchy,
  width = w_multi_lr_hierarchy,
  height = h_multi_lr_hierarchy
)

par(mfrow = c(1, 1), xpd = TRUE)

netVisual_individual(
  cellchat_obj,
  color.use = cluster_colors_multi,
  signaling = multi_signal,
  pairLR.use = multi_lr_selected,
  layout = "hierarchy",
  vertex.receiver = target_cell_indices
)

dev.off()

cat("多信号 LR 层级图保存成功：", file_multi_lr_hierarchy, "\n")


# 31. 多细胞多信号作图板块：多信号气泡图

# 31.1 多信号气泡图参数
w_multi_bubble <- 8
h_multi_bubble <- 8
name_multi_bubble <- "11.多个信号、配体对通路水平号的细胞通讯可视化"

if (length(signal_source) == 0) {
  stop("signal_source 为空，请至少设置一个发送细胞。")
}

if (length(receiver_cells) == 0) {
  stop("receiver_cells 为空，请至少设置一个接收细胞。")
}

# 31.2 提取用于气泡图的多信号配体-受体对
multi_signal_lr_for_bubble <- extractEnrichedLR(
  cellchat_obj,
  signaling = multi_signal
)

# 31.3 绘制多信号气泡图
p_multi_bubble <- netVisual_bubble(
  cellchat_obj,
  sources.use = signal_source,
  targets.use = receiver_cells,
  pairLR.use = multi_signal_lr_for_bubble,
  remove.isolate = FALSE
)

# 31.4 显示多信号气泡图
print(p_multi_bubble)

# 31.5 保存多信号气泡图
file_multi_bubble <- file.path(
  multi_dir,
  paste0(name_multi_bubble, ".pdf")
)

pdf(
  file = file_multi_bubble,
  width = w_multi_bubble,
  height = h_multi_bubble
)

print(p_multi_bubble)

dev.off()

cat("多信号气泡图保存成功：", file_multi_bubble, "\n")


# 32. 保存参数记录

# 32.1 参数文件名
name_params <- "cellchat_plot_parameters"

# 32.2 整理参数文本
param_text <- paste0(
  "当前作图参数总结：
",
  "- cellchat 来源：全局环境变量 cellchat
",
  "- 可展示的信号通路：", paste(pathways, collapse = ", "), "
",
"- 细胞类型：", paste(cell_types, collapse = ", "), "
",
"- 单信号通路 single_signal：", single_signal, "
",
"- 单信号配体受体对 single_lr：", single_lr, "
",
"- 单信号层级图关注细胞 single_target_cells：", paste(single_target_cells, collapse = ", "), "
",
"- 多信号通路 multi_signal：", paste(multi_signal, collapse = ", "), "
",
"- 多信号配体受体对 multi_lr：", multi_lr, "
",
"- 多信号层级图关注细胞 target_cells：", paste(target_cells, collapse = ", "), "
",
"- 气泡图发送细胞 signal_source：", paste(signal_source, collapse = ", "), "
",
"- 气泡图接收细胞 receiver_cells：", paste(receiver_cells, collapse = ", "), "
",
"- 单信号热图配色 single_heatmap_color：", single_heatmap_color, "
",
"- 单信号角色网络图字体大小 single_font_size：", single_font_size, "
",
"- 单信号角色网络图内部宽度参数：", single_role_network_width, "
",
"- 单信号角色网络图内部高度参数：", single_role_network_height, "
",
"- 汇总图颜色 palette_count_text：", palette_count_text, "
",
"- 汇总图颜色 palette_weight_text：", palette_weight_text, "
",
"- 汇总图颜色 palette_sender_text：", palette_sender_text, "
",
"- 单信号图颜色 palette_single_text：", palette_single_text, "
",
"- 多信号图颜色 palette_multi_text：", palette_multi_text, "
",
"- 输出目录：", out_dir
)

# 32.3 保存参数文件
param_file <- file.path(
  out_dir,
  paste0(name_params, ".txt")
)

writeLines(
  param_text,
  con = param_file
)

cat("参数保存成功：", param_file, "\n")


# 33. 完成提示

cat("CellChat 总作图流程运行完成。\n")
cat("汇总图输出目录：", summary_dir, "\n")
cat("单信号图输出目录：", single_dir, "\n")
cat("多细胞多信号图输出目录：", multi_dir, "\n")