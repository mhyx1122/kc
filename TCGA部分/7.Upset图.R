# 1. 加载必要 R 包

library(UpSetR)
library(dplyr)
library(tidyr)
library(ggsci)


# 2. 设置输入和输出参数

input_file <- "your_input.csv"

output_dir <- "5.基因取交集"

color_point <- "#e76c54"

pdf_width <- 8

pdf_height <- 6


# 3. 创建输出目录

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


# 4. 读取输入数据

mcc_data <- read.csv(
  input_file,
  header = TRUE,
  check.names = FALSE
)


# 5. 整理每一列中的元素和所属集合

all_elements <- data.frame(
  Element = character(),
  Hub = character()
)

for (hub in names(mcc_data)) {
  
  column_data <- na.omit(mcc_data[[hub]])
  
  column_data <- column_data[column_data != ""]
  
  temp <- data.frame(
    Element = column_data,
    Hub = hub
  )
  
  all_elements <- rbind(
    all_elements,
    temp
  )
}

all_elements <- na.omit(all_elements)


# 6. 转换为 UpSetR 需要的 0/1 矩阵格式

all_elements_upset <- all_elements %>%
  mutate(
    incidence = 1
  ) %>%
  spread(
    Hub,
    incidence,
    fill = 0
  )


# 7. 根据交集组合数量生成颜色

intersections <- apply(
  all_elements_upset[, -1],
  1,
  function(x) {
    paste(
      which(x == 1),
      collapse = ","
    )
  }
)

unique_intersections <- unique(intersections)

randomColor123 <- function() {
  paste0(
    "#",
    paste0(
      sample(
        c(0:9, letters[1:6]),
        6,
        replace = TRUE
      ),
      collapse = ""
    )
  )
}

randomColors123 <- replicate(
  length(unique_intersections),
  randomColor123()
)

intersection_colors <- randomColors123


# 8. 设置集合柱状图颜色

set_colors <- c(
  ggsci::pal_npg()(9),
  ggsci::pal_jco()(9),
  ggsci::pal_jama()(7),
  ggsci::pal_nejm()(8)
)

set_colors <- rep(
  set_colors,
  length.out = length(names(mcc_data))
)


# 9. 保存 UpSet 图为 PDF

upset_pdf_path <- file.path(
  output_dir,
  "upset_plot_colored_elements.pdf"
)

pdf(
  upset_pdf_path,
  width = pdf_width,
  height = pdf_height
)

upset(
  all_elements_upset,
  sets = names(mcc_data),
  order.by = "degree",
  matrix.color = color_point,
  main.bar.color = intersection_colors,
  sets.bar.color = set_colors,
  sets.x.label = "Set Size",
  point.size = 3,
  line.size = 1,
  mb.ratio = c(0.6, 0.4),
  shade.color = "gray80",
  shade.alpha = 0.3,
  matrix.dot.alpha = 0.6,
  show.numbers = "yes",
  number.angles = 0,
  group.by = "degree",
  text.scale = 1.2,
  set_size.angles = 0
)

dev.off()


# 10. 保存 UpSet 矩阵数据

all_elements_upset_path <- file.path(
  output_dir,
  "all_elements_upset.csv"
)

write.csv(
  all_elements_upset,
  all_elements_upset_path,
  row.names = FALSE
)


# 11. 计算每个元素出现于多少个集合

int <- all_elements_upset

rownames(int) <- int[[1]]

int$row_sums <- rowSums(
  int[, -1]
)

all_upset_path <- file.path(
  output_dir,
  "all__upset.csv"
)

write.csv(
  int,
  all_upset_path,
  row.names = FALSE
)


# 12. 输出结果路径

cat("UpSet 图已保存：", upset_pdf_path, "\n")

cat("UpSet 矩阵数据已保存：", all_elements_upset_path, "\n")

cat("带 row_sums 的结果已保存：", all_upset_path, "\n")