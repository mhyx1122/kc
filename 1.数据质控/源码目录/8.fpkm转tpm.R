# 1. 设置输入文件路径
fpkm_data_file <- "选择的FPKM表达量表格文件.csv"

# 2. 读取FPKM数据文件
expMatrix <- read.csv(fpkm_data_file, header=TRUE, row.names=1)

# 3. 定义FPKM到TPM的转换函数
fpkmToTpm <- function(fpkm) {
  return(exp(log(fpkm) - log(sum(fpkm)) + log(1e6)))
}

# 4. 应用转换函数将 FPKM 转换为 TPM
tpms <- apply(expMatrix, 2, fpkmToTpm)

# 5. 转换后的数据框
tpms <- as.data.frame(tpms)

# 6. 打印每个样本的总和，检查转换是否成功
message("作者提示：如果每个样本的总和加起来是100万（10的6次方，在R中显示为1e+06），代表转换成功！")
message("计算结果如下：")
print(colSums(tpms))

# 7. 创建文件夹用于保存结果（如果不存在的话）
if (!dir.exists("FPKM转TPM结果")) {
  dir.create("FPKM转TPM结果")
}

# 8. 保存转换后的结果为CSV文件
write.csv(tpms, file = "FPKM转TPM结果/tpm_genes.csv", row.names = TRUE)

# 9. 输出处理完成信息
cat("数据转换完成！结果文件已保存。\n")