# 1. 加载必要的 R 包

suppressPackageStartupMessages({
  library(CellChat)
  library(NMF)
  library(ggalluvial)
})


# 2. 检查 cellchat 对象

if (!exists("cellchat", envir = .GlobalEnv)) {
  stop("全局环境中没有找到 cellchat 对象，请先完成 CellChat 分析。")
}

cellchat <- get("cellchat", envir = .GlobalEnv)

if (!inherits(cellchat, "CellChat")) {
  stop("全局环境中的 cellchat 不是 CellChat 对象。")
}


# 3. 设置输出文件夹

out_dir <- "细胞通讯模式"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}


# 4. 识别细胞通讯模式

# 模式类型：
# "outgoing" 表示发送模式
# "incoming" 表示接收模式
pattern_mod <- "outgoing"

# 根据上一阶段 selectK() 的结果设置模式数量
nPatterns <- 6

# 模式识别图字体大小
font_pattern <- 8

# 模式识别图保存参数
pattern_width <- 15
pattern_height <- 12
pattern_file <- "2.确定模式数量.pdf"

message("正在识别细胞通讯模式，请耐心等待运行完毕")

pdf(
  file = file.path(out_dir, pattern_file),
  width = pattern_width,
  height = pattern_height
)

cellchat_pattern <- identifyCommunicationPatterns(
  cellchat,
  pattern = pattern_mod,
  k = nPatterns,
  font.size = font_pattern
)

dev.off()

message(
  "细胞通讯模式识别完成，结果图已保存至：",
  file.path(out_dir, pattern_file)
)


# 5. 绘制细胞通讯模式桑基图

font_river <- 2.5
font_river_title <- 12

river_width <- 15
river_height <- 12
river_file <- "3.通讯模式桑基图.pdf"

message("正在生成细胞通讯模式桑基图")

p_river <- netAnalysis_river(
  cellchat_pattern,
  pattern = pattern_mod,
  font.size = font_river,
  font.size.title = font_river_title
)

pdf(
  file = file.path(out_dir, river_file),
  width = river_width,
  height = river_height
)

print(p_river)

dev.off()

message(
  "细胞通讯模式桑基图已保存至：",
  file.path(out_dir, river_file)
)


# 6. 绘制细胞通讯模式泡泡图

font_dot <- 10
font_dot_title <- 12

dot_min <- 1
dot_max <- 3

dot_width <- 15
dot_height <- 12
dot_file <- "4.通讯模式泡泡图.pdf"

message("正在生成细胞通讯模式泡泡图")

p_dot <- netAnalysis_dot(
  cellchat_pattern,
  pattern = pattern_mod,
  font.size = font_dot,
  dot.size = c(dot_min, dot_max),
  font.size.title = font_dot_title
)

pdf(
  file = file.path(out_dir, dot_file),
  width = dot_width,
  height = dot_height
)

print(p_dot)

dev.off()

message(
  "细胞通讯模式泡泡图已保存至：",
  file.path(out_dir, dot_file)
)


# 7. 保存本次分析参数

params_file <- "communication_pattern_parameters.txt"

param_text <- paste0(
  "本次通讯模式分析参数总结：\n",
  "- 模式类型：", pattern_mod, "\n",
  "- 模式数量：", nPatterns, "\n",
  "- 模式识别图字体大小：", font_pattern, "\n",
  "- 桑基图字体大小：", font_river, "\n",
  "- 桑基图标题字体大小：", font_river_title, "\n",
  "- 泡泡图字体大小：", font_dot, "\n",
  "- 泡泡图标题字体大小：", font_dot_title, "\n",
  "- 泡泡图点大小范围：", dot_min, " ~ ", dot_max, "\n",
  "- 输出文件夹：", out_dir, "\n",
  "- 结果对象：cellchat_pattern\n",
  "\n",
  "写作提示词（自行组装语言，或借助AI组装）：\n",
  "1. 基于 CellChat 识别细胞通讯中的主要发送或接收模式。\n",
  "2. 通过桑基图展示细胞群、信号通路与通讯模式之间的对应关系。\n",
  "3. 通过泡泡图展示不同通讯模式在信号通路层面的贡献特征。\n"
)

writeLines(
  text = param_text,
  con = file.path(out_dir, params_file)
)

cat(param_text)

message(
  "分析参数已保存至：",
  file.path(out_dir, params_file)
)
