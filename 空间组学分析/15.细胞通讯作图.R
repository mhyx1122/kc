suppressPackageStartupMessages({
  library(CellChat)
  library(ggplot2)
})

# 1. 读取CellChat对象并准备公共信息

# 1.1 参数设置

out_dir <- "细胞通讯分析"

default_palette <- c(
  "#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF",
  "#91D1C2FF", "#7E6148FF", "#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF",
  "#7AA6DCFF", "#003C67FF", "#8F7700FF", "#3B3B3BFF", "#A73030FF", "#374E55FF",
  "#DF8F44FF", "#00A1D5FF", "#B24745FF", "#79AF97FF", "#6A6599FF", "#80796BFF",
  "#BC3C29FF", "#0072B5FF", "#E18727FF", "#20854EFF", "#7876B1FF", "#6F99ADFF",
  "#FFDC91FF", "#EE4C97FF"
)

# 1.2 检查并读取cellchat对象

if (!exists("cellchat", envir = .GlobalEnv)) {
  stop("全局环境中没有 cellchat 对象，请先完成CellChat通讯分析。")
}

cellchat_obj <- get("cellchat", envir = .GlobalEnv)

if (!inherits(cellchat_obj, "CellChat")) {
  stop("全局环境中的cellchat不是有效的CellChat对象。")
}

if (is.null(cellchat_obj@net$count) || is.null(cellchat_obj@net$weight)) {
  stop("cellchat对象中不存在net$count或net$weight，请先完成aggregateNet()。")
}

if (is.null(cellchat_obj@netP$pathways)) {
  stop("cellchat对象中不存在netP$pathways，请先完成computeCommunProbPathway()。")
}

pathways <- cellchat_obj@netP$pathways
cell_types <- levels(cellchat_obj@idents)
group_size <- as.numeric(table(cellchat_obj@idents))

if (length(pathways) == 0) {
  stop("cellchat对象中没有可用于作图的信号通路。")
}

if (length(cell_types) == 0) {
  stop("cellchat对象中没有有效的细胞分组信息。")
}

if (length(default_palette) < length(cell_types)) {
  stop(
    paste0(
      "默认颜色数量不足。当前有",
      length(cell_types),
      "个细胞群，但仅提供了",
      length(default_palette),
      "个颜色。"
    )
  )
}

# 1.3 创建输出文件夹

summary_output_folder <- file.path(out_dir, "细胞通讯汇总图")
single_output_folder <- file.path(out_dir, "单信号作图")
multi_output_folder <- file.path(out_dir, "多细胞多信号作图")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(summary_output_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(single_output_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(multi_output_folder, recursive = TRUE, showWarnings = FALSE)

cat(
  "当前细胞群：", paste(cell_types, collapse = ", "), "\n",
  "当前可用信号通路：", paste(pathways, collapse = ", "), "\n",
  sep = ""
)

# 2. 绘制数据库分类图

# 2.1 参数设置

w_db_category <- 10
h_db_category <- 7
name_db_category <- "showDatabaseCategory"

# 2.2 绘制数据库分类图

db_category_plot <- showDatabaseCategory(
  cellchat_obj@DB
)

print(db_category_plot)

# 2.3 保存数据库分类图

db_category_file <- file.path(
  out_dir,
  paste0(name_db_category, ".pdf")
)

ggsave(
  filename = db_category_file,
  plot = db_category_plot,
  width = w_db_category,
  height = h_db_category,
  device = "pdf"
)

# 3. 绘制细胞通讯网络连接数量图

# 3.1 参数设置

palette_count <- default_palette

w_net_count <- 8
h_net_count <- 8
name_net_count <- "1.细胞-细胞通讯网络(连接数量)"

# 3.2 准备细胞颜色

if (length(palette_count) < length(cell_types)) {
  stop("palette_count中的颜色数量不足以覆盖所有细胞群。")
}

cluster_colors_count <- setNames(
  palette_count[seq_along(cell_types)],
  cell_types
)

# 3.3 保存连接数量网络图

net_count_file <- file.path(
  summary_output_folder,
  paste0(name_net_count, ".pdf")
)

pdf(
  file = net_count_file,
  width = w_net_count,
  height = h_net_count
)

par(xpd = TRUE)

netVisual_circle(
  cellchat_obj@net$count,
  vertex.weight = group_size,
  weight.scale = TRUE,
  color.use = cluster_colors_count,
  label.edge = TRUE,
  title.name = "Number of interactions"
)

dev.off()

# 4. 绘制细胞通讯网络连接强度图

# 4.1 参数设置

palette_weight <- default_palette

w_net_weight <- 8
h_net_weight <- 8
name_net_weight <- "2.细胞-细胞通讯网络(连接强度)"

# 4.2 准备细胞颜色

if (length(palette_weight) < length(cell_types)) {
  stop("palette_weight中的颜色数量不足以覆盖所有细胞群。")
}

cluster_colors_weight <- setNames(
  palette_weight[seq_along(cell_types)],
  cell_types
)

# 4.3 保存连接强度网络图

net_weight_file <- file.path(
  summary_output_folder,
  paste0(name_net_weight, ".pdf")
)

pdf(
  file = net_weight_file,
  width = w_net_weight,
  height = h_net_weight
)

par(xpd = TRUE)

netVisual_circle(
  cellchat_obj@net$weight,
  vertex.weight = group_size,
  weight.scale = TRUE,
  color.use = cluster_colors_weight,
  label.edge = TRUE,
  title.name = "Interaction weights/strength"
)

dev.off()

# 5. 绘制每个细胞群发送的信号

# 5.1 参数设置

palette_sender <- default_palette

signal_rows <- 3
signal_cols <- 4

w_sender_signal <- 12
h_sender_signal <- 10
name_sender_signal <- "3.单个细胞群发送的信号"

# 5.2 准备细胞颜色和通讯矩阵

if (length(palette_sender) < length(cell_types)) {
  stop("palette_sender中的颜色数量不足以覆盖所有细胞群。")
}

cluster_colors_sender <- setNames(
  palette_sender[seq_along(cell_types)],
  cell_types
)

signal_matrix <- cellchat_obj@net$weight

# 5.3 保存每个细胞群发送的信号图

sender_signal_file <- file.path(
  summary_output_folder,
  paste0(name_sender_signal, ".pdf")
)

pdf(
  file = sender_signal_file,
  width = w_sender_signal,
  height = h_sender_signal
)

par(
  mfrow = c(signal_rows, signal_cols),
  xpd = TRUE
)

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

# 6. 设置单信号作图参数并选择配体受体对

# 6.1 参数设置

# 可从pathways中选择一个信号通路
single_signal <- pathways[1]

# 设置为NULL时，默认选择该信号通路中的第一个配体受体对
# 也可设置为具体显示名称，例如：
# single_lr_label <- "MIF - (CD74+CXCR4)"
single_lr_label <- NULL

# 单信号层级图关注的接收细胞
single_target_cells <- cell_types[
  seq_len(min(2, length(cell_types)))
]

# 单信号热图配色，可选：
# Reds、Blues、Greens、Purples、Oranges、Greys
single_heatmap_color <- "Reds"

# 信号角色网络图内部参数
single_font_size <- 10
single_role_network_width <- 8
single_role_network_height <- 2.5

palette_single <- default_palette

# 6.2 检查单信号参数

if (!single_signal %in% pathways) {
  stop(
    paste0(
      "single_signal不在当前可用信号通路中。可选值为：",
      paste(pathways, collapse = ", ")
    )
  )
}

missing_single_target_cells <- setdiff(
  single_target_cells,
  cell_types
)

if (length(missing_single_target_cells) > 0) {
  stop(
    paste0(
      "以下single_target_cells不存在：",
      paste(missing_single_target_cells, collapse = ", ")
    )
  )
}

if (length(single_target_cells) == 0) {
  stop("single_target_cells不能为空。")
}

if (length(palette_single) < length(cell_types)) {
  stop("palette_single中的颜色数量不足以覆盖所有细胞群。")
}

cluster_colors_single <- setNames(
  palette_single[seq_along(cell_types)],
  cell_types
)

single_target_indices <- match(
  single_target_cells,
  cell_types
)

single_target_indices <- single_target_indices[
  !is.na(single_target_indices)
]

# 6.3 提取单信号中的配体受体对

single_lr_df <- extractEnrichedLR(
  cellchat_obj,
  signaling = single_signal,
  geneLR.return = FALSE
)

if (is.null(single_lr_df) || nrow(single_lr_df) == 0) {
  stop(
    paste0(
      "信号通路 ",
      single_signal,
      " 中没有可用的配体受体对。"
    )
  )
}

# 6.4 生成配体受体显示名称

if ("interaction_name_2" %in% colnames(single_lr_df)) {
  single_lr_labels <- as.character(
    single_lr_df$interaction_name_2
  )
} else if ("interaction_name" %in% colnames(single_lr_df)) {
  single_lr_labels <- as.character(
    single_lr_df$interaction_name
  )
} else {
  single_lr_labels <- apply(
    single_lr_df,
    1,
    function(x) {
      paste(x, collapse = " | ")
    }
  )
}

empty_single_labels <- is.na(single_lr_labels) |
  single_lr_labels == ""

single_lr_labels[empty_single_labels] <- paste0(
  "LR_",
  which(empty_single_labels)
)

# 6.5 根据名称选择配体受体对

if (
  is.null(single_lr_label) ||
  length(single_lr_label) == 0 ||
  trimws(single_lr_label) == ""
) {
  single_lr_index <- 1
  single_lr_label <- single_lr_labels[1]
} else {
  single_lr_index <- match(
    single_lr_label,
    single_lr_labels
  )
  
  if (is.na(single_lr_index)) {
    warning(
      paste0(
        "未找到single_lr_label：",
        single_lr_label,
        "，将使用第一个配体受体对：",
        single_lr_labels[1]
      )
    )
    
    single_lr_index <- 1
    single_lr_label <- single_lr_labels[1]
  }
}

selected_single_lr <- single_lr_df[
  single_lr_index,
  ,
  drop = FALSE
]

cat(
  "单信号通路：", single_signal, "\n",
  "单信号配体受体对：", single_lr_label, "\n",
  "单信号层级图关注细胞：",
  paste(single_target_cells, collapse = ", "),
  "\n",
  sep = ""
)

# 7. 绘制单信号圈图

# 7.1 参数设置

w_single_circle <- 15
h_single_circle <- 8
name_single_circle <- "3.细胞通讯_特定信号"

# 7.2 保存单信号圈图

single_circle_file <- file.path(
  single_output_folder,
  paste0(name_single_circle, ".pdf")
)

pdf(
  file = single_circle_file,
  width = w_single_circle,
  height = h_single_circle
)

par(
  mfrow = c(1, 1),
  xpd = TRUE
)

netVisual_aggregate(
  cellchat_obj,
  color.use = cluster_colors_single,
  signaling = single_signal,
  layout = "circle"
)

dev.off()

# 8. 绘制单信号弦图

# 8.1 参数设置

w_single_chord <- 15
h_single_chord <- 8
name_single_chord <- "4.细胞通讯_特定信号弦图"

# 8.2 保存单信号弦图

single_chord_file <- file.path(
  single_output_folder,
  paste0(name_single_chord, ".pdf")
)

pdf(
  file = single_chord_file,
  width = w_single_chord,
  height = h_single_chord
)

par(
  mfrow = c(1, 1),
  xpd = TRUE
)

netVisual_aggregate(
  cellchat_obj,
  color.use = cluster_colors_single,
  signaling = single_signal,
  layout = "chord"
)

dev.off()

# 9. 绘制单信号热图

# 9.1 参数设置

w_single_heatmap <- 15
h_single_heatmap <- 8
name_single_heatmap <- "5.细胞通讯_特定信号热图"

# 9.2 生成单信号热图

single_heatmap_plot <- netVisual_heatmap(
  cellchat_obj,
  color.use = cluster_colors_single,
  signaling = single_signal,
  color.heatmap = single_heatmap_color
)

print(single_heatmap_plot)

# 9.3 保存单信号热图

single_heatmap_file <- file.path(
  single_output_folder,
  paste0(name_single_heatmap, ".pdf")
)

pdf(
  file = single_heatmap_file,
  width = w_single_heatmap,
  height = h_single_heatmap
)

print(single_heatmap_plot)

dev.off()

# 10. 绘制单信号层级图

# 10.1 参数设置

w_hierarchy_single <- 25
h_hierarchy_single <- 8
name_hierarchy_single <- "9.细胞通讯_单信号下的层级图"

# 10.2 保存单信号层级图

hierarchy_single_file <- file.path(
  single_output_folder,
  paste0(name_hierarchy_single, ".pdf")
)

pdf(
  file = hierarchy_single_file,
  width = w_hierarchy_single,
  height = h_hierarchy_single
)

par(
  mfrow = c(1, 1),
  xpd = TRUE
)

tryCatch(
  {
    netVisual_aggregate(
      cellchat_obj,
      color.use = cluster_colors_single,
      vertex.receiver = single_target_indices,
      signaling = single_signal,
      layout = "hierarchy"
    )
  },
  error = function(e) {
    print(
      ggplot() +
        theme_void() +
        ggtitle(
          paste(
            "层级图绘制失败：",
            e$message
          )
        )
    )
  }
)

dev.off()

# 11. 计算信号中心性并绘制角色网络图

# 11.1 参数设置

w_single_role_network <- 8
h_single_role_network <- 8
name_single_role_network <- "6.1细胞群在该信号中的角色"

# 11.2 计算信号通路网络中心性

cellchat_centrality <- netAnalysis_computeCentrality(
  cellchat_obj,
  slot.name = "netP"
)

# 11.3 保存信号角色网络图

single_role_network_file <- file.path(
  single_output_folder,
  paste0(name_single_role_network, ".pdf")
)

pdf(
  file = single_role_network_file,
  width = w_single_role_network,
  height = h_single_role_network
)

par(xpd = TRUE)

netAnalysis_signalingRole_network(
  cellchat_centrality,
  color.use = cluster_colors_single,
  signaling = single_signal,
  width = single_role_network_width,
  height = single_role_network_height,
  font.size = single_font_size
)

dev.off()

# 12. 绘制信号角色散点图

# 12.1 参数设置

w_single_role_scatter <- 6
h_single_role_scatter <- 5
name_single_role_scatter <- "6.2网络中信号角色在指定通路中的权重"

# 12.2 生成信号角色散点图

single_role_scatter_plot <- netAnalysis_signalingRole_scatter(
  cellchat_centrality,
  color.use = cluster_colors_single,
  signaling = single_signal
)

print(single_role_scatter_plot)

# 12.3 保存信号角色散点图

single_role_scatter_file <- file.path(
  single_output_folder,
  paste0(name_single_role_scatter, ".pdf")
)

ggsave(
  filename = single_role_scatter_file,
  plot = single_role_scatter_plot,
  width = w_single_role_scatter,
  height = h_single_role_scatter,
  device = "pdf"
)

# 13. 绘制单信号下指定配体受体对圈图

# 13.1 参数设置

w_single_lr_circle <- 15
h_single_lr_circle <- 8
name_single_lr_circle <- "10.1单信号下_配受体对水平的细胞通讯_圈图"

# 13.2 保存配体受体对圈图

single_lr_circle_file <- file.path(
  single_output_folder,
  paste0(name_single_lr_circle, ".pdf")
)

pdf(
  file = single_lr_circle_file,
  width = w_single_lr_circle,
  height = h_single_lr_circle
)

par(xpd = TRUE)

netVisual_individual(
  cellchat_obj,
  color.use = cluster_colors_single,
  signaling = single_signal,
  pairLR.use = selected_single_lr,
  layout = "circle"
)

dev.off()

# 14. 绘制单信号下指定配体受体对弦图

# 14.1 参数设置

w_single_lr_chord <- 15
h_single_lr_chord <- 8
name_single_lr_chord <- "10.2单信号下_配受体对水平的细胞通讯_弦图"

# 14.2 保存配体受体对弦图

single_lr_chord_file <- file.path(
  single_output_folder,
  paste0(name_single_lr_chord, ".pdf")
)

pdf(
  file = single_lr_chord_file,
  width = w_single_lr_chord,
  height = h_single_lr_chord
)

par(xpd = TRUE)

netVisual_individual(
  cellchat_obj,
  color.use = cluster_colors_single,
  signaling = single_signal,
  pairLR.use = selected_single_lr,
  layout = "chord"
)

dev.off()

# 15. 绘制单信号下指定配体受体对层级图

# 15.1 参数设置

w_single_lr_hierarchy <- 15
h_single_lr_hierarchy <- 8
name_single_lr_hierarchy <- "10.3单信号下_配受体对水平的细胞通讯_层级图"

# 15.2 保存配体受体对层级图

single_lr_hierarchy_file <- file.path(
  single_output_folder,
  paste0(name_single_lr_hierarchy, ".pdf")
)

pdf(
  file = single_lr_hierarchy_file,
  width = w_single_lr_hierarchy,
  height = h_single_lr_hierarchy
)

par(xpd = TRUE)

tryCatch(
  {
    netVisual_individual(
      cellchat_obj,
      color.use = cluster_colors_single,
      signaling = single_signal,
      pairLR.use = selected_single_lr,
      layout = "hierarchy",
      vertex.receiver = single_target_indices
    )
  },
  error = function(e) {
    print(
      ggplot() +
        theme_void() +
        ggtitle(
          paste(
            "层级图绘制失败：",
            e$message
          )
        )
    )
  }
)

dev.off()

# 16. 绘制单信号配体受体贡献排名

# 16.1 参数设置

w_single_contribution <- 8
h_single_contribution <- 8
name_single_contribution <- "7.单信号：信号通路-配体受体贡献排名"

# 16.2 生成贡献排名图

single_contribution_plot <- netAnalysis_contribution(
  cellchat_obj,
  signaling = single_signal
)

print(single_contribution_plot)

# 16.3 保存贡献排名图

single_contribution_file <- file.path(
  single_output_folder,
  paste0(name_single_contribution, ".pdf")
)

pdf(
  file = single_contribution_file,
  width = w_single_contribution,
  height = h_single_contribution
)

print(single_contribution_plot)

dev.off()

# 17. 设置多信号作图参数并选择配体受体对

# 17.1 参数设置

# 默认选择前三个信号通路，可手动修改
multi_signal <- pathways[
  seq_len(min(3, length(pathways)))
]

# 设置为NULL时，默认选择多信号结果中的第一个配体受体对
multi_lr_label <- NULL

# 多信号层级图关注的接收细胞
target_cells <- cell_types[
  seq_len(min(2, length(cell_types)))
]

# 气泡图发送细胞
signal_source <- cell_types[1]

# 气泡图接收细胞
receiver_cells <- cell_types[
  seq_len(min(2, length(cell_types)))
]

palette_multi <- default_palette

# 17.2 检查多信号和细胞参数

if (length(multi_signal) == 0) {
  stop("multi_signal不能为空。")
}

missing_multi_signals <- setdiff(
  multi_signal,
  pathways
)

if (length(missing_multi_signals) > 0) {
  stop(
    paste0(
      "以下multi_signal不在可用通路中：",
      paste(missing_multi_signals, collapse = ", ")
    )
  )
}

missing_target_cells <- setdiff(
  target_cells,
  cell_types
)

if (length(missing_target_cells) > 0) {
  stop(
    paste0(
      "以下target_cells不存在：",
      paste(missing_target_cells, collapse = ", ")
    )
  )
}

missing_signal_source <- setdiff(
  signal_source,
  cell_types
)

if (length(missing_signal_source) > 0) {
  stop(
    paste0(
      "以下signal_source不存在：",
      paste(missing_signal_source, collapse = ", ")
    )
  )
}

missing_receiver_cells <- setdiff(
  receiver_cells,
  cell_types
)

if (length(missing_receiver_cells) > 0) {
  stop(
    paste0(
      "以下receiver_cells不存在：",
      paste(missing_receiver_cells, collapse = ", ")
    )
  )
}

if (length(target_cells) == 0) {
  stop("target_cells不能为空。")
}

if (length(signal_source) == 0) {
  stop("signal_source不能为空。")
}

if (length(receiver_cells) == 0) {
  stop("receiver_cells不能为空。")
}

if (length(palette_multi) < length(cell_types)) {
  stop("palette_multi中的颜色数量不足以覆盖所有细胞群。")
}

cluster_colors_multi <- setNames(
  palette_multi[seq_along(cell_types)],
  cell_types
)

target_cell_indices <- match(
  target_cells,
  cell_types
)

target_cell_indices <- target_cell_indices[
  !is.na(target_cell_indices)
]

# 17.3 提取多信号中的配体受体对

multi_lr_df <- extractEnrichedLR(
  cellchat_obj,
  signaling = multi_signal,
  geneLR.return = FALSE
)

if (is.null(multi_lr_df) || nrow(multi_lr_df) == 0) {
  stop("所选multi_signal中没有可用的配体受体对。")
}

# 17.4 生成多信号配体受体显示名称

if ("interaction_name_2" %in% colnames(multi_lr_df)) {
  multi_lr_labels <- as.character(
    multi_lr_df$interaction_name_2
  )
} else if ("interaction_name" %in% colnames(multi_lr_df)) {
  multi_lr_labels <- as.character(
    multi_lr_df$interaction_name
  )
} else {
  multi_lr_labels <- apply(
    multi_lr_df,
    1,
    function(x) {
      paste(x, collapse = " | ")
    }
  )
}

empty_multi_labels <- is.na(multi_lr_labels) |
  multi_lr_labels == ""

multi_lr_labels[empty_multi_labels] <- paste0(
  "LR_",
  which(empty_multi_labels)
)

# 17.5 根据名称选择多信号配体受体对

if (
  is.null(multi_lr_label) ||
  length(multi_lr_label) == 0 ||
  trimws(multi_lr_label) == ""
) {
  multi_lr_index <- 1
  multi_lr_label <- multi_lr_labels[1]
} else {
  multi_lr_index <- match(
    multi_lr_label,
    multi_lr_labels
  )
  
  if (is.na(multi_lr_index)) {
    warning(
      paste0(
        "未找到multi_lr_label：",
        multi_lr_label,
        "，将使用第一个配体受体对：",
        multi_lr_labels[1]
      )
    )
    
    multi_lr_index <- 1
    multi_lr_label <- multi_lr_labels[1]
  }
}

selected_multi_lr <- multi_lr_df[
  multi_lr_index,
  ,
  drop = FALSE
]

cat(
  "多信号通路：", paste(multi_signal, collapse = ", "), "\n",
  "多信号配体受体对：", multi_lr_label, "\n",
  "层级图关注细胞：", paste(target_cells, collapse = ", "), "\n",
  "气泡图发送细胞：", paste(signal_source, collapse = ", "), "\n",
  "气泡图接收细胞：", paste(receiver_cells, collapse = ", "), "\n",
  sep = ""
)

# 18. 绘制多信号圈图

# 18.1 参数设置

w_multi_circle <- 15
h_multi_circle <- 8
name_multi_circle <- "8.1多信号_圈图"

# 18.2 保存多信号圈图

multi_circle_file <- file.path(
  multi_output_folder,
  paste0(name_multi_circle, ".pdf")
)

pdf(
  file = multi_circle_file,
  width = w_multi_circle,
  height = h_multi_circle
)

par(
  mfrow = c(1, 1),
  xpd = TRUE
)

netVisual_aggregate(
  cellchat_obj,
  color.use = cluster_colors_multi,
  signaling = multi_signal,
  layout = "circle"
)

dev.off()

# 19. 绘制多信号弦图

# 19.1 参数设置

w_multi_chord <- 15
h_multi_chord <- 8
name_multi_chord <- "8.2多信号_弦图"

# 19.2 保存多信号弦图

multi_chord_file <- file.path(
  multi_output_folder,
  paste0(name_multi_chord, ".pdf")
)

pdf(
  file = multi_chord_file,
  width = w_multi_chord,
  height = h_multi_chord
)

par(
  mfrow = c(1, 1),
  xpd = TRUE
)

netVisual_aggregate(
  cellchat_obj,
  color.use = cluster_colors_multi,
  signaling = multi_signal,
  layout = "chord"
)

dev.off()

# 20. 绘制多信号贡献排名图

# 20.1 参数设置

w_multi_contribution <- 8
h_multi_contribution <- 8
name_multi_contribution <- "7.多信号：信号通路-配体受体贡献排名"

# 20.2 生成多信号贡献排名图

multi_contribution_plot <- netAnalysis_contribution(
  cellchat_obj,
  signaling = multi_signal
)

print(multi_contribution_plot)

# 20.3 保存多信号贡献排名图

multi_contribution_file <- file.path(
  multi_output_folder,
  paste0(name_multi_contribution, ".pdf")
)

pdf(
  file = multi_contribution_file,
  width = w_multi_contribution,
  height = h_multi_contribution
)

print(multi_contribution_plot)

dev.off()

# 21. 绘制多信号相关基因表达图

# 21.1 参数设置

w_multi_gene_expr <- 8
h_multi_gene_expr <- 15
name_multi_gene_expr <- "8.多信号：信号通路中基因在细胞亚群中的表达量"

# 21.2 生成多信号基因表达图

multi_gene_expr_plot <- suppressWarnings(
  plotGeneExpression(
    cellchat_obj,
    color.use = cluster_colors_multi,
    signaling = multi_signal
  )
)

print(multi_gene_expr_plot)

# 21.3 保存多信号基因表达图

multi_gene_expr_file <- file.path(
  multi_output_folder,
  paste0(name_multi_gene_expr, ".pdf")
)

pdf(
  file = multi_gene_expr_file,
  width = w_multi_gene_expr,
  height = h_multi_gene_expr
)

suppressWarnings(
  print(multi_gene_expr_plot)
)

dev.off()

# 22. 绘制多信号层级图

# 22.1 参数设置

w_hierarchy_multi <- 25
h_hierarchy_multi <- 8
name_hierarchy_multi <- "9.细胞通讯_多信号下的层级图"

# 22.2 保存多信号层级图

hierarchy_multi_file <- file.path(
  multi_output_folder,
  paste0(name_hierarchy_multi, ".pdf")
)

pdf(
  file = hierarchy_multi_file,
  width = w_hierarchy_multi,
  height = h_hierarchy_multi
)

par(
  mfrow = c(1, 1),
  xpd = TRUE
)

tryCatch(
  {
    netVisual_aggregate(
      cellchat_obj,
      color.use = cluster_colors_multi,
      vertex.receiver = target_cell_indices,
      signaling = multi_signal,
      layout = "hierarchy"
    )
  },
  error = function(e) {
    print(
      ggplot() +
        theme_void() +
        ggtitle(
          paste(
            "层级图绘制失败：",
            e$message
          )
        )
    )
  }
)

dev.off()

# 23. 绘制多信号下指定配体受体对圈图

# 23.1 参数设置

w_multi_lr_circle <- 15
h_multi_lr_circle <- 8
name_multi_lr_circle <- "10.1多信号下_配受体对水平的细胞通讯_圈图"

# 23.2 保存多信号配体受体对圈图

multi_lr_circle_file <- file.path(
  multi_output_folder,
  paste0(name_multi_lr_circle, ".pdf")
)

pdf(
  file = multi_lr_circle_file,
  width = w_multi_lr_circle,
  height = h_multi_lr_circle
)

par(xpd = TRUE)

netVisual_individual(
  cellchat_obj,
  color.use = cluster_colors_multi,
  signaling = multi_signal,
  pairLR.use = selected_multi_lr,
  layout = "circle"
)

dev.off()

# 24. 绘制多信号下指定配体受体对弦图

# 24.1 参数设置

w_multi_lr_chord <- 15
h_multi_lr_chord <- 8
name_multi_lr_chord <- "10.2多信号下_配受体对水平的细胞通讯_弦图"

# 24.2 保存多信号配体受体对弦图

multi_lr_chord_file <- file.path(
  multi_output_folder,
  paste0(name_multi_lr_chord, ".pdf")
)

pdf(
  file = multi_lr_chord_file,
  width = w_multi_lr_chord,
  height = h_multi_lr_chord
)

par(xpd = TRUE)

netVisual_individual(
  cellchat_obj,
  color.use = cluster_colors_multi,
  signaling = multi_signal,
  pairLR.use = selected_multi_lr,
  layout = "chord"
)

dev.off()

# 25. 绘制多信号下指定配体受体对层级图

# 25.1 参数设置

w_multi_lr_hierarchy <- 15
h_multi_lr_hierarchy <- 8
name_multi_lr_hierarchy <- "10.3多信号下_配受体对水平的细胞通讯_层级图"

# 25.2 保存多信号配体受体对层级图

multi_lr_hierarchy_file <- file.path(
  multi_output_folder,
  paste0(name_multi_lr_hierarchy, ".pdf")
)

pdf(
  file = multi_lr_hierarchy_file,
  width = w_multi_lr_hierarchy,
  height = h_multi_lr_hierarchy
)

par(xpd = TRUE)

tryCatch(
  {
    netVisual_individual(
      cellchat_obj,
      color.use = cluster_colors_multi,
      signaling = multi_signal,
      pairLR.use = selected_multi_lr,
      layout = "hierarchy",
      vertex.receiver = target_cell_indices
    )
  },
  error = function(e) {
    print(
      ggplot() +
        theme_void() +
        ggtitle(
          paste(
            "层级图绘制失败：",
            e$message
          )
        )
    )
  }
)

dev.off()

# 26. 绘制多信号配体受体气泡图

# 26.1 参数设置

w_multi_bubble <- 8
h_multi_bubble <- 8
name_multi_bubble <- "11.多个信号、配体对通路水平号的细胞通讯可视化"

# 26.2 提取多信号配体受体对

multi_signal_lr_for_bubble <- extractEnrichedLR(
  cellchat_obj,
  signaling = multi_signal
)

if (
  is.null(multi_signal_lr_for_bubble) ||
  nrow(multi_signal_lr_for_bubble) == 0
) {
  stop("所选multi_signal中没有可用于气泡图的配体受体对。")
}

# 26.3 生成多信号气泡图

multi_bubble_plot <- netVisual_bubble(
  cellchat_obj,
  sources.use = signal_source,
  targets.use = receiver_cells,
  pairLR.use = multi_signal_lr_for_bubble,
  remove.isolate = FALSE
)

print(multi_bubble_plot)

# 26.4 保存多信号气泡图

multi_bubble_file <- file.path(
  multi_output_folder,
  paste0(name_multi_bubble, ".pdf")
)

pdf(
  file = multi_bubble_file,
  width = w_multi_bubble,
  height = h_multi_bubble
)

print(multi_bubble_plot)

dev.off()

# 27. 保存作图参数记录

# 27.1 参数设置

name_params <- "cellchat_plot_parameters"

# 27.2 生成参数记录

param_text <- paste0(
  "本次CellChat作图参数总结：\n",
  "- cellchat来源：全局环境变量cellchat\n",
  "- 细胞群：", paste(cell_types, collapse = ", "), "\n",
  "- 可展示的信号通路：", paste(pathways, collapse = ", "), "\n",
  "- 输出目录：", out_dir, "\n\n",
  
  "【汇总作图参数】\n",
  "- 细胞群发送信号排版：", signal_rows, "行 × ", signal_cols, "列\n",
  "- 汇总图目录：", summary_output_folder, "\n\n",
  
  "【单信号作图参数】\n",
  "- single_signal：", single_signal, "\n",
  "- single_lr_label：", single_lr_label, "\n",
  "- single_target_cells：", paste(single_target_cells, collapse = ", "), "\n",
  "- single_heatmap_color：", single_heatmap_color, "\n",
  "- single_font_size：", single_font_size, "\n",
  "- single_role_network_width：", single_role_network_width, "\n",
  "- single_role_network_height：", single_role_network_height, "\n",
  "- 单信号图目录：", single_output_folder, "\n\n",
  
  "【多信号作图参数】\n",
  "- multi_signal：", paste(multi_signal, collapse = ", "), "\n",
  "- multi_lr_label：", multi_lr_label, "\n",
  "- target_cells：", paste(target_cells, collapse = ", "), "\n",
  "- signal_source：", paste(signal_source, collapse = ", "), "\n",
  "- receiver_cells：", paste(receiver_cells, collapse = ", "), "\n",
  "- 多信号图目录：", multi_output_folder, "\n\n",
  
  "说明：\n",
  "1. 汇总作图包括数据库分类、通讯数量、通讯强度和各细胞群发送信号。\n",
  "2. 单信号作图包括圈图、弦图、热图、层级图、信号角色图及配体受体对图。\n",
  "3. 多信号作图包括圈图、弦图、贡献排名、基因表达、层级图和气泡图。\n",
  "4. 所有细胞颜色均由代码中的颜色向量控制，不使用随机颜色。"
)

# 27.3 保存参数文件

param_file <- file.path(
  out_dir,
  paste0(name_params, ".txt")
)

writeLines(
  param_text,
  con = param_file
)