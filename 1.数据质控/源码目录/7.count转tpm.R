# 1. 设置输入文件路径
counts_geneexp_file <- "选择的表达矩阵文件.csv"
counts_idType <- "Symbol"  # 基因标识符类型
counts_org <- "hsa"  # 生物体选择，hsa为人类，mmus为小鼠

# 2. 读取表达矩阵文件
data <- read.csv(counts_geneexp_file, row.names = 1)

# 3. 输出检查信息
print("如果报错，请先检查idType（基因名类型）和org（物种缩写）是不是写错了")

# 4. 调用 count2tpm 函数将 counts 转换为 TPM 格式
# 请确保你已经加载了适当的包，并且 `count2tpm` 函数已正确定义
tpm <- count2tpm(data, idType = counts_idType, org = counts_org)

# 5. 保存转换后的数据
write.csv(tpm, "counts转换为TPM数据结果/tpm转换后的数据.csv", row.names = TRUE)

# 6. 假设你需要保存参考基因长度的基因组信息（如 grch38）
# 这里的 `anno_grch38` 应该是一个包含基因长度的基因组信息表格
write.csv(anno_grch38, "counts转换为TPM数据结果/grch38基因组（参考的基因长度信息）.csv", row.names = TRUE)