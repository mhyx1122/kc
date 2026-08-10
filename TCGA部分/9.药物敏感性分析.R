# 1. 加载必要 R 包

library(limma)
library(ggpubr)
library(pRRophetic)
library(ggplot2)


# 2. 加载 pRRophetic 内置药物数据库

data(cgp2016ExprRma)
data(PANCANCER_IC_Tue_Aug_9_15_28_57_2016)

allDrugs <- unique(drugData2016$Drug.name)


# 3. 模块一：目标基因与全部药物敏感性分析

input_file <- "after_group_TCGA.csv"

gene_Drugnames <- c("RALA")

row_mean_filter <- 0.5

cor_method <- "pearson"
# 可选：
# cor_method <- "spearman"

pFilter <- 0.05

corFilter <- 0.3

plot_width <- 5

plot_height <- 5

output_folder <- "药物敏感性分析结果"

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

dataKK145 <- read.csv(
  input_file,
  header = TRUE,
  sep = ",",
  check.names = FALSE
)

gene <- gene_Drugnames

rt <- dataKK145

rt <- as.matrix(rt)

rownames(rt) <- rt[, 1]

exp <- rt[, 2:ncol(rt)]

dimnames <- list(
  rownames(exp),
  colnames(exp)
)

data <- matrix(
  as.numeric(as.matrix(exp)),
  nrow = nrow(exp),
  dimnames = dimnames
)

data <- avereps(data)

data <- data[
  rowMeans(data) > row_mean_filter,
]

data <- t(data)

rownames(data) <- substr(
  rownames(data),
  1,
  12
)

data <- t(avereps(data))

geneExp <- as.data.frame(
  t(data[gene, , drop = FALSE])
)

geneExp$Type <- ifelse(
  geneExp[, gene] > median(geneExp[, gene]),
  "High",
  "Low"
)

result_df <- data.frame(
  Drug = character(),
  Pvalue = numeric(),
  Correlation = numeric(),
  stringsAsFactors = FALSE
)

for (drug in allDrugs) {
  
  possibleError <- tryCatch({
    senstivity <- pRRopheticPredict(
      data,
      drug,
      selection = 1,
      dataset = "cgp2016"
    )
  }, error = function(e) e)
  
  if (inherits(possibleError, "error")) {
    next
  }
  
  senstivity <- senstivity[!is.nan(senstivity)]

  senstivity[
    senstivity > quantile(senstivity, 0.99)
  ] <- quantile(senstivity, 0.99)
  
  sameSample <- intersect(
    row.names(geneExp),
    names(senstivity)
  )
  
  geneExp1 <- geneExp[
    sameSample,
    "Type",
    drop = FALSE
  ]
  
  geneExp2 <- geneExp[
    sameSample,
    gene,
    drop = FALSE
  ]
  
  senstivity <- senstivity[sameSample]
  
  rt <- cbind(
    geneExp1,
    senstivity
  )
  
  rt$Type <- factor(
    rt$Type,
    levels = c("Low", "High")
  )
  
  type <- levels(
    factor(rt[, "Type"])
  )
  
  comp <- combn(type, 2)
  
  my_comparisons <- list()
  
  for (i in 1:ncol(comp)) {
    my_comparisons[[i]] <- comp[, i]
  }
  
  test <- wilcox.test(
    senstivity ~ Type,
    data = rt
  )
  
  diffPvalue <- test$p.value
  
  x <- as.numeric(
    geneExp2[, gene]
  )
  
  y <- as.numeric(
    rt$senstivity
  )
  
  corT <- cor.test(
    x,
    y,
    method = cor_method
  )
  
  cor <- corT$estimate
  
  pvalue <- corT$p.value
  
  result_df <- rbind(
    result_df,
    data.frame(
      Drug = drug,
      Pvalue = diffPvalue,
      Correlation = cor
    )
  )
  
  if (diffPvalue < pFilter) {
    
    boxplot <- ggboxplot(
      rt,
      x = "Type",
      y = "senstivity",
      fill = "Type",
      xlab = gene,
      ylab = paste0(
        drug,
        " sensitivity (IC50)"
      ),
      palette = c("#0066FF", "#FF0000")
    ) +
      labs(
        fill = gene
      ) +
      stat_compare_means(
        comparisons = my_comparisons
      )
    
    drugname <- gsub(
      "/",
      "_",
      drug
    )
    
    pdf(
      file = paste0(
        output_folder,
        "/drugSensitivity_Boxplot_",
        drugname,
        ".pdf"
      ),
      width = plot_width,
      height = plot_height
    )
    
    print(boxplot)
    
    dev.off()
  }
  
  if (abs(cor) > corFilter & pvalue < pFilter) {
    
    df1 <- as.data.frame(
      cbind(
        x,
        y
      )
    )
    
    p1 <- ggplot(
      df1,
      aes(
        x,
        y
      )
    ) +
      xlab(
        paste0(
          gene,
          " expression"
        )
      ) +
      ylab(
        paste0(
          drug,
          " drug sensitivity (IC50)"
        )
      ) +
      geom_point() +
      geom_smooth(
        method = "lm",
        formula = y ~ x
      ) +
      theme_bw() +
      stat_cor(
        method = cor_method,
        aes(
          x = x,
          y = y
        )
      )
    
    drugname <- gsub(
      "/",
      "_",
      drug
    )
    
    pdf(
      file = paste0(
        output_folder,
        "/drugSensitivity_Correlation_",
        drugname,
        ".pdf"
      ),
      width = plot_width,
      height = plot_height
    )
    
    print(p1)
    
    dev.off()
  }
}

write.csv(
  result_df,
  paste0(
    output_folder,
    "/drug_sensitivity_all_drugs_results.csv"
  ),
  row.names = FALSE
)


# 4. 模块二：指定药物和指定基因作图

input_file <- "after_group_TCGA.csv"

drug_fixed <- "Axitinib"

target_gene <- "RALA"

row_mean_filter <- 0.5

cor_method <- "pearson"
# 可选：
# cor_method <- "spearman"

plot_width <- 5

plot_height <- 5

output_folder <- "药物敏感性分析结果"

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

dataKK145 <- read.csv(
  input_file,
  header = TRUE,
  sep = ",",
  check.names = FALSE
)

rt <- dataKK145

rt <- as.matrix(rt)

rownames(rt) <- rt[, 1]

exp <- rt[, 2:ncol(rt)]

dimnames <- list(
  rownames(exp),
  colnames(exp)
)

data <- matrix(
  as.numeric(as.matrix(exp)),
  nrow = nrow(exp),
  dimnames = dimnames
)

data <- avereps(data)

data <- data[
  rowMeans(data) > row_mean_filter,
]

data <- t(data)

rownames(data) <- substr(
  rownames(data),
  1,
  12
)

data <- t(avereps(data))

possibleError <- tryCatch({
  sensitivity <- pRRopheticPredict(
    data,
    drug_fixed,
    selection = 1,
    dataset = "cgp2016"
  )
}, error = function(e) e)

if (inherits(possibleError, "error")) {
  stop("药物敏感性预测失败！")
}

sensitivity <- sensitivity[!is.nan(sensitivity)]

sensitivity[
  sensitivity > quantile(sensitivity, 0.99)
] <- quantile(sensitivity, 0.99)

geneExp <- as.data.frame(
  t(data[target_gene, , drop = FALSE])
)

geneExp$Type <- ifelse(
  geneExp[, target_gene] > median(geneExp[, target_gene]),
  "High",
  "Low"
)

sameSample <- intersect(
  rownames(geneExp),
  names(sensitivity)
)

geneExp1 <- geneExp[
  sameSample,
  "Type",
  drop = FALSE
]

geneExp2 <- geneExp[
  sameSample,
  target_gene,
  drop = FALSE
]

sens <- sensitivity[sameSample]

rt <- cbind(
  geneExp1,
  sens
)

rt$Type <- factor(
  rt$Type,
  levels = c("Low", "High")
)

type <- levels(rt$Type)

comp <- combn(type, 2)

my_comparisons <- list()

for (i in 1:ncol(comp)) {
  my_comparisons[[i]] <- comp[, i]
}

test <- wilcox.test(
  sens ~ Type,
  data = rt
)

diffPvalue <- test$p.value

x <- as.numeric(
  geneExp2[, target_gene]
)

y <- as.numeric(
  rt$sens
)

corT <- cor.test(
  x,
  y,
  method = cor_method
)

cor <- corT$estimate

pvalue <- corT$p.value

boxplot <- ggboxplot(
  rt,
  x = "Type",
  y = "sens",
  fill = "Type",
  xlab = target_gene,
  ylab = paste0(
    drug_fixed,
    " sensitivity (IC50)"
  ),
  palette = c("#0066FF", "#FF0000")
) +
  labs(
    fill = target_gene
  ) +
  stat_compare_means(
    comparisons = my_comparisons
  )

geneName <- gsub(
  "/",
  "_",
  target_gene
)

pdf(
  file = paste0(
    output_folder,
    "/drugSensitivity_Boxplot_",
    geneName,
    ".pdf"
  ),
  width = plot_width,
  height = plot_height
)

print(boxplot)

dev.off()

df1 <- as.data.frame(
  cbind(
    x,
    y
  )
)

p1 <- ggplot(
  df1,
  aes(
    x,
    y
  )
) +
  xlab(
    paste0(
      target_gene,
      " expression"
    )
  ) +
  ylab(
    paste0(
      drug_fixed,
      " drug sensitivity (IC50)"
    )
  ) +
  geom_point() +
  geom_smooth(
    method = "lm",
    formula = y ~ x
  ) +
  theme_bw() +
  stat_cor(
    method = cor_method,
    aes(
      x = x,
      y = y
    )
  )

geneName <- gsub(
  "/",
  "_",
  target_gene
)

pdf(
  file = paste0(
    output_folder,
    "/drugSensitivity_Correlation_",
    geneName,
    ".pdf"
  ),
  width = plot_width,
  height = plot_height
)

print(p1)

dev.off()


# 5. 模块三：指定药物和全部基因批量相关性结果

input_file <- "after_group_TCGA.csv"

drug_fixed <- "Axitinib"

row_mean_filter <- 0.5

cor_method <- "pearson"
# 可选：
# cor_method <- "spearman"

output_folder <- "药物敏感性分析结果"

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

dataKK145 <- read.csv(
  input_file,
  header = TRUE,
  sep = ",",
  check.names = FALSE
)

rt <- dataKK145

rt <- as.matrix(rt)

rownames(rt) <- rt[, 1]

exp <- rt[, 2:ncol(rt)]

dimnames <- list(
  rownames(exp),
  colnames(exp)
)

data <- matrix(
  as.numeric(as.matrix(exp)),
  nrow = nrow(exp),
  dimnames = dimnames
)

data <- avereps(data)

data <- data[
  rowMeans(data) > row_mean_filter,
]

data <- t(data)

rownames(data) <- substr(
  rownames(data),
  1,
  12
)

data <- t(avereps(data))

possibleError <- tryCatch({
  sensitivity <- pRRopheticPredict(
    data,
    drug_fixed,
    selection = 1,
    dataset = "cgp2016"
  )
}, error = function(e) e)

if (inherits(possibleError, "error")) {
  stop("药物敏感性预测失败！")
}

sensitivity <- sensitivity[!is.nan(sensitivity)]

sensitivity[
  sensitivity > quantile(sensitivity, 0.99)
] <- quantile(sensitivity, 0.99)

sensdata <- as.data.frame(
  sensitivity
)

write.csv(
  sensdata,
  file = paste0(
    output_folder,
    "/sensitivity_data.csv"
  ),
  row.names = TRUE
)

result_df <- data.frame(
  Gene = character(),
  Pvalue = numeric(),
  Correlation = numeric(),
  stringsAsFactors = FALSE
)

for (gene in rownames(data)) {
  
  geneExp <- as.data.frame(
    t(data[gene, , drop = FALSE])
  )
  
  geneExp$Type <- ifelse(
    geneExp[, gene] > median(geneExp[, gene]),
    "High",
    "Low"
  )
  
  sameSample <- intersect(
    rownames(geneExp),
    names(sensitivity)
  )
  
  geneExp1 <- geneExp[
    sameSample,
    "Type",
    drop = FALSE
  ]
  
  geneExp2 <- geneExp[
    sameSample,
    gene,
    drop = FALSE
  ]
  
  sens <- sensitivity[sameSample]
  
  rt <- cbind(
    geneExp1,
    sens
  )
  
  rt$Type <- factor(
    rt$Type,
    levels = c("Low", "High")
  )
  
  type <- levels(rt$Type)
  
  comp <- combn(type, 2)
  
  my_comparisons <- list()
  
  for (i in 1:ncol(comp)) {
    my_comparisons[[i]] <- comp[, i]
  }
  
  test <- wilcox.test(
    sens ~ Type,
    data = rt
  )
  
  diffPvalue <- test$p.value
  
  x <- as.numeric(
    geneExp2[, gene]
  )
  
  y <- as.numeric(
    rt$sens
  )
  
  corT <- cor.test(
    x,
    y,
    method = cor_method
  )
  
  cor <- corT$estimate
  
  pvalue <- corT$p.value
  
  result_df <- rbind(
    result_df,
    data.frame(
      Gene = gene,
      Pvalue = diffPvalue,
      Correlation = cor
    )
  )
}

write.csv(
  result_df,
  paste0(
    output_folder,
    "/drug_sensitivity_all_genes_results.csv"
  ),
  row.names = FALSE
)
