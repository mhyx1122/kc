# 1. 加载必要的 R 包

suppressPackageStartupMessages({
  library(CellChat)
  library(NMF)
})


# 2. 提示分析开始

message("正在识别细胞通讯模式，请耐心等待运行完毕")


# 3. 创建结果保存文件夹

output_folder <- "细胞通讯模式"

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}


# 4. 识别细胞通讯模式数量并保存结果

pdf(
  file = file.path(output_folder, "1.识别细胞模式数量.pdf"),
  width = 12,
  height = 6
)

print(
  selectK(
    cellchat,
    pattern = c("outgoing", "incoming")
  )
)

dev.off()