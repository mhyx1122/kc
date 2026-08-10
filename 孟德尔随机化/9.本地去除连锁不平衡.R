suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
})

# 1、读取输入文件夹中的CSV文件

# 输入文件夹路径
folder_input <- "exposure_data"

if (!dir.exists(folder_input)) {
  stop("输入文件夹不存在：", folder_input)
}

csv_files <- list.files(
  path = folder_input,
  pattern = "\\.csv$",
  full.names = TRUE
)

if (length(csv_files) == 0) {
  stop("输入文件夹中没有找到CSV文件：", folder_input)
}


# 2、设置Clumping参数

# Clumping窗口大小，单位为kb
clump_kb <- 10000

# Clumping的R²阈值
clump_r2 <- 0.001

# Clumping的P值阈值
clump_p1 <- 1

# LD参考人群，可设置为EUR、SAS、EAS、AFR、AMR或legacy
pop <- "EUR"


# 3、创建结果保存文件夹

# 结果保存文件夹
folder_path <- "clump_file"

if (!dir.exists(folder_path)) {
  dir.create(folder_path, recursive = TRUE)
}


# 4、依次读取每个CSV文件并执行Clumping

for (i in seq_along(csv_files)) {
  file_path <- csv_files[i]
  
  message(
    "[", i, "/", length(csv_files), "] 正在处理：",
    basename(file_path)
  )
  
  exposure_datacg <- data.table::fread(file_path)
  
  exposure_data <- TwoSampleMR::clump_data(
    dat = exposure_datacg,
    clump_kb = clump_kb,
    clump_r2 = clump_r2,
    clump_p1 = clump_p1,
    pop = pop,
    bfile = NULL,
    plink_bin = NULL
  )
  
  output_file <- file.path(
    folder_path,
    paste0("clumped_", basename(file_path))
  )
  
  write.csv(
    exposure_data,
    file = output_file,
    row.names = FALSE
  )
  
  message(
    "Clumping完成，保留SNP数量：",
    nrow(exposure_data)
  )
}


# 5、输出完成信息

message("所有CSV文件均已处理完成。")
message(
  "结果保存位置：",
  normalizePath(folder_path, winslash = "/", mustWork = FALSE)
)