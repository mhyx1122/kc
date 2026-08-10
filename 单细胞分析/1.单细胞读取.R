library(Seurat)
library(Matrix)
library(data.table)

# 1. 选择要读取的数据类型
# 可选值："10X"、"H5"、"counts"
data_type <- "10X"

if (data_type == "10X") {
  
  # 2. 读取 10X 数据模块
  
  # 2.1 10X 数据模块参数
  min_cells <- 3
  min_features <- 200
  data_dir_10x <- paste0(getwd(), "/10Xdata")
  
  # 2.2 获取所有样本目录
  sample_dirs <- list.files(data_dir_10x)
  sample_paths <- file.path(data_dir_10x, sample_dirs)
  
  # 2.3 逐个样本读取 10X 数据并创建 Seurat 对象
  seurat_list <- list()
  
  for (sample_path in sample_paths) {
    sample_name <- basename(sample_path)
    
    counts_matrix <- Read10X(data.dir = sample_path)
    
    seurat_obj <- CreateSeuratObject(
      counts = counts_matrix,
      project = sample_name,
      min.cells = min_cells,
      min.features = min_features
    )
    
    seurat_obj$orig.ident <- sample_name
    
    cat("已读取：")
    print(table(seurat_obj$orig.ident))
    
    seurat_list[[sample_name]] <- seurat_obj
  }
  
  # 2.4 合并所有 10X 样本
  seurat <- merge(
    seurat_list[[1]],
    y = seurat_list[-1],
    add.cell.ids = names(seurat_list),
    project = "Combined_Seurat"
  )
  
  cat("共读取如下样本，每个样本的细胞数为：")
  print(table(seurat$orig.ident))
  
} else if (data_type == "H5") {
  
  # 3. 读取 H5 数据模块
  
  # 3.1 H5 数据模块参数
  min_cells <- 3
  min_features <- 200
  data_dir_h5 <- paste0(getwd(), "/h5data")
  
  # 3.2 获取所有 h5 文件
  h5_files <- list.files(data_dir_h5, pattern = "\\.h5$", full.names = TRUE)
  
  # 3.3 逐个 h5 文件读取数据并创建 Seurat 对象
  seurat_list <- list()
  
  for (h5_file in h5_files) {
    counts_matrix <- Read10X_h5(h5_file, use.names = TRUE)
    
    sample_name <- gsub(".h5", "", basename(h5_file))
    
    seurat_obj <- CreateSeuratObject(
      counts = counts_matrix,
      project = sample_name,
      min.cells = min_cells,
      min.features = min_features
    )
    
    cat("已读取：")
    print(table(seurat_obj$orig.ident))
    
    seurat_list[[sample_name]] <- seurat_obj
  }
  
  # 3.4 合并所有 H5 样本
  seurat <- merge(
    seurat_list[[1]],
    y = seurat_list[-1],
    add.cell.ids = names(seurat_list),
    project = "Combined_Seurat"
  )
  
  cat("共读取如下样本，每个样本的细胞数为：")
  print(table(seurat$orig.ident))
  
} else if (data_type == "counts") {
  
  # 4. 读取 counts 矩阵模块
  
  # 4.1 counts 数据模块参数
  min_cells <- 3
  min_features <- 200
  count_dir <- "counts"
  
  # 4.2 获取所有 csv 格式的 counts 文件
  count_files <- list.files(count_dir, pattern = "\\.csv$", full.names = TRUE)
  
  # 4.3 逐个 counts 文件读取数据并创建 Seurat 对象
  seurat_list <- list()
  
  for (count_file in count_files) {
    count_data <- fread(count_file, header = TRUE)
    count_data <- as.data.frame(count_data)
    
    count_data <- count_data[!duplicated(count_data[, 1]), ]
    rownames(count_data) <- count_data[, 1]
    count_data <- count_data[, -1]
    
    sparse_matrix <- as(as.matrix(count_data), "dgCMatrix")
    
    sample_name <- gsub("\\.csv$", "", basename(count_file))
    
    seurat_obj <- CreateSeuratObject(
      counts = sparse_matrix,
      project = sample_name,
      min.cells = min_cells,
      min.features = min_features
    )
    
    cat("已读取：")
    print(table(seurat_obj$orig.ident))
    
    seurat_list[[basename(count_file)]] <- seurat_obj
  }
  
  # 4.4 合并所有 counts 样本
  seurat <- merge(
    seurat_list[[1]],
    y = seurat_list[-1],
    add.cell.ids = names(seurat_list),
    project = "Combined_Seurat"
  )
  
  cat("共读取如下样本，每个样本的细胞数为：")
  print(table(seurat$orig.ident))
  
} else {
  
  # 5. 数据类型错误提示
  stop("未知的数据类型，请选择 '10X', 'H5' 或 'counts'.")
}

# 6. 查看 Seurat 对象概况
summary(seurat)