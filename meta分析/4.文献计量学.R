suppressPackageStartupMessages({
  library(bibliometrix)
})

# 1、导入文献数据

# 文献数据库导出的文件路径，可填写一个文件或多个文件
input_files <- c("savedrecs.txt")

# 数据库来源，可设置为wos、scopus、pubmed、cochrane、dimensions或generic
dbsource <- "wos"

# 文件格式，可设置为plaintext、bibtex或csv
file_format <- "plaintext"

if (length(input_files) == 0 || any(!file.exists(input_files))) {
  stop("没有找到以下输入文件：", paste(input_files[!file.exists(input_files)], collapse = "、"))
}

M <- bibliometrix::convert2df(file = input_files, dbsource = dbsource, format = file_format)

if (is.null(M) || nrow(M) == 0) {
  stop("导入后数据为空，请检查文件路径、数据库来源和文件格式。")
}

message("文献数据导入完成：", nrow(M), "篇文献，", ncol(M), "个字段。")
print(head(M[, seq_len(min(10, ncol(M))), drop = FALSE], 10))

# 2、创建结果保存文件夹

# 结果保存文件夹
out_dir <- "文献计量学分析"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# 3、进行基础文献计量分析

# summary()和plot()显示前K项
summary_k <- 5

if (summary_k < 1) {
  stop("summary_k必须大于等于1。")
}

results <- bibliometrix::biblioAnalysis(M, sep = ";")
S <- summary(results, k = summary_k, pause = FALSE, verbose = FALSE)
summary_text <- paste(capture.output(print(S)), collapse = "\n")

cat(summary_text, "\n")
writeLines(summary_text, file.path(out_dir, "基础文献计量分析结果.txt"))

# plot(results)会主动绘图，先用临时PDF接收绘图，再保留返回的ggplot列表
temp_overview_pdf <- tempfile(fileext = ".pdf")
grDevices::pdf(temp_overview_pdf, width = 10, height = 7)
overview_plots <- try(plot(results, k = summary_k, pause = FALSE), silent = TRUE)
grDevices::dev.off()
unlink(temp_overview_pdf)

if (inherits(overview_plots, "try-error")) {
  message("基础概览图生成失败：", as.character(overview_plots))
  overview_plots <- NULL
} else {
  overview_plots <- overview_plots[!vapply(overview_plots, is.null, logical(1))]
}

# 4、导出基础信息并生成结果解读

basic_info_df <- data.frame(n_records = nrow(M), n_fields = ncol(M))
write.csv(M, file.path(out_dir, "导入后的文献数据.csv"), row.names = FALSE)
write.csv(basic_info_df, file.path(out_dir, "基础信息.csv"), row.names = FALSE)
writeLines(paste(names(M), collapse = ", "), file.path(out_dir, "字段名称.txt"))

interpretation_text <- c(
  paste0("1. 当前数据共纳入 ", nrow(M), " 篇文献，包含 ", ncol(M), " 个字段。")
)

if ("PY" %in% names(M)) {
  publication_year <- suppressWarnings(as.numeric(as.character(M$PY)))
  publication_year <- publication_year[!is.na(publication_year)]
  
  if (length(publication_year) > 0) {
    interpretation_text <- c(
      interpretation_text,
      paste0(
        "2. 文献发表时间范围约为 ", min(publication_year), " 至 ", max(publication_year),
        "，可结合年发文趋势图判断该领域是否处于增长阶段。"
      )
    )
  }
}

if ("AU" %in% names(M) && any(nchar(trimws(as.character(M$AU))) > 0, na.rm = TRUE)) {
  interpretation_text <- c(
    interpretation_text,
    "3. 作者合作网络可用于观察核心作者群体及合作紧密程度。若网络较分散，说明该领域作者合作可能相对松散。"
  )
}

if ("DE" %in% names(M) && any(nchar(trimws(as.character(M$DE))) > 0, na.rm = TRUE)) {
  interpretation_text <- c(
    interpretation_text,
    "4. 关键词共现网络可用于识别当前研究热点。节点较大、连接较密的关键词通常代表热点主题或核心概念。"
  )
} else if ("ID" %in% names(M) && any(nchar(trimws(as.character(M$ID))) > 0, na.rm = TRUE)) {
  interpretation_text <- c(
    interpretation_text,
    "4. Keywords Plus或ID字段可辅助识别研究热点与主题结构。"
  )
}

if ("CR" %in% names(M) && any(nchar(trimws(as.character(M$CR))) > 0, na.rm = TRUE)) {
  interpretation_text <- c(
    interpretation_text,
    "5. 共被引网络和历史引文图谱有助于识别该领域的知识基础、经典文献和发展脉络。"
  )
}

if ("TC" %in% names(M)) {
  interpretation_text <- c(
    interpretation_text,
    "6. 高被引文献通常代表领域内影响力较高的研究，但应结合发表时间、研究类型和主题背景综合判断。"
  )
}

interpretation_text <- paste(interpretation_text, collapse = "\n")
cat(interpretation_text, "\n")
writeLines(interpretation_text, file.path(out_dir, "结果解读.txt"))

# 5、保存基础概览图

# 保存第几张基础概览图
overview_index <- 1

# PDF保存参数
overview_width <- 10
overview_height <- 7
overview_current_file <- "0.基础概览图_当前.pdf"
overview_all_file <- "0.基础概览图_全部.pdf"

if (!is.null(overview_plots) && length(overview_plots) > 0) {
  overview_index <- max(1, min(as.integer(overview_index), length(overview_plots)))
  
  grDevices::pdf(
    file.path(out_dir, overview_current_file),
    width = overview_width,
    height = overview_height
  )
  print(overview_plots[[overview_index]])
  grDevices::dev.off()
  
  grDevices::pdf(
    file.path(out_dir, overview_all_file),
    width = overview_width,
    height = overview_height,
    onefile = TRUE
  )
  
  for (i in seq_along(overview_plots)) {
    print(overview_plots[[i]])
  }
  
  grDevices::dev.off()
} else {
  message("没有可保存的基础概览图。")
}

# 6、绘制年发文趋势图

# 年发文趋势图参数
annual_type <- "b"
annual_pch <- 19
annual_color <- "#E64B35FF"
annual_xlab <- "Year"
annual_ylab <- "Number of Publications"
annual_main <- "Annual Scientific Production"

# PDF保存参数
annual_width <- 10
annual_height <- 7
annual_file <- "1.年发文趋势图.pdf"

if (!"PY" %in% names(M)) {
  message("跳过年发文趋势图：当前数据缺少PY字段。")
} else {
  annual_year <- suppressWarnings(as.numeric(as.character(M$PY)))
  annual_year <- annual_year[!is.na(annual_year)]
  
  if (length(annual_year) == 0) {
    message("跳过年发文趋势图：PY字段没有有效年份。")
  } else {
    annual_prod <- table(annual_year)
    annual_x <- as.numeric(names(annual_prod))
    annual_y <- as.numeric(annual_prod)
    annual_order <- order(annual_x)
    
    grDevices::pdf(file.path(out_dir, annual_file), width = annual_width, height = annual_height)
    
    plot(
      annual_x[annual_order],
      annual_y[annual_order],
      type = annual_type,
      pch = annual_pch,
      xlab = annual_xlab,
      ylab = annual_ylab,
      main = annual_main,
      col = annual_color
    )
    
    grDevices::dev.off()
  }
}

# 7、绘制国家合作网络

# 国家合作网络参数
country_type <- "circle"
country_labelsize <- 1
country_node_size <- 3
country_size_by_frequency <- TRUE

# PDF保存参数
country_width <- 10
country_height <- 8
country_file <- "2.国家合作网络.pdf"

if (!("C1" %in% names(M) || "RP" %in% names(M))) {
  message("跳过国家合作网络：当前数据缺少C1和RP字段。")
} else {
  M_country <- try(
    bibliometrix::metaTagExtraction(M, Field = "AU_CO", sep = ";"),
    silent = TRUE
  )
  
  if (
    inherits(M_country, "try-error") ||
    !"AU_CO" %in% names(M_country) ||
    !any(nchar(trimws(as.character(M_country$AU_CO))) > 0, na.rm = TRUE)
  ) {
    message("跳过国家合作网络：未成功生成有效的AU_CO字段。")
  } else {
    net_country <- try(
      bibliometrix::biblioNetwork(
        M_country,
        analysis = "collaboration",
        network = "countries",
        sep = ";"
      ),
      silent = TRUE
    )
    
    if (inherits(net_country, "try-error")) {
      message("国家合作网络构建失败：", as.character(net_country))
    } else {
      grDevices::pdf(
        file.path(out_dir, country_file),
        width = country_width,
        height = country_height
      )
      
      country_plot <- try(
        bibliometrix::networkPlot(
          net_country,
          type = country_type,
          Title = "国家合作网络",
          size = country_node_size,
          size.cex = country_size_by_frequency,
          labelsize = country_labelsize,
          verbose = TRUE
        ),
        silent = TRUE
      )
      
      grDevices::dev.off()
      
      if (inherits(country_plot, "try-error")) {
        message("国家合作网络绘制失败：", as.character(country_plot))
      }
    }
  }
}

# 8、绘制作者合作网络

# 作者合作网络参数
author_n <- 30
author_type <- "fruchterman"
author_labelsize <- 0.8
author_node_size <- 3
author_size_by_frequency <- TRUE
author_remove_isolates <- TRUE

# PDF保存参数
author_width <- 10
author_height <- 8
author_file <- "3.作者合作网络.pdf"

if (!"AU" %in% names(M) || !any(nchar(trimws(as.character(M$AU))) > 0, na.rm = TRUE)) {
  message("跳过作者合作网络：当前数据缺少有效的AU字段。")
} else {
  net_author <- try(
    bibliometrix::biblioNetwork(
      M,
      analysis = "collaboration",
      network = "authors",
      sep = ";"
    ),
    silent = TRUE
  )
  
  if (inherits(net_author, "try-error")) {
    message("作者合作网络构建失败：", as.character(net_author))
  } else {
    grDevices::pdf(
      file.path(out_dir, author_file),
      width = author_width,
      height = author_height
    )
    
    author_plot <- try(
      bibliometrix::networkPlot(
        net_author,
        n = author_n,
        type = author_type,
        Title = "作者合作网络",
        size = author_node_size,
        size.cex = author_size_by_frequency,
        remove.isolates = author_remove_isolates,
        labelsize = author_labelsize,
        verbose = TRUE
      ),
      silent = TRUE
    )
    
    grDevices::dev.off()
    
    if (inherits(author_plot, "try-error")) {
      message("作者合作网络绘制失败：", as.character(author_plot))
    }
  }
}

# 9、绘制机构合作网络

# 机构合作网络参数
inst_n <- 30
inst_type <- "fruchterman"
inst_labelsize <- 0.8
inst_node_size <- 3
inst_size_by_frequency <- TRUE
inst_remove_isolates <- TRUE

# PDF保存参数
inst_width <- 10
inst_height <- 8
inst_file <- "4.机构合作网络.pdf"

if (!"C1" %in% names(M) || !any(nchar(trimws(as.character(M$C1))) > 0, na.rm = TRUE)) {
  message("跳过机构合作网络：当前数据缺少有效的C1字段。")
} else {
  net_inst <- try(
    bibliometrix::biblioNetwork(
      M,
      analysis = "collaboration",
      network = "universities",
      sep = ";"
    ),
    silent = TRUE
  )
  
  if (inherits(net_inst, "try-error")) {
    message("机构合作网络构建失败：", as.character(net_inst))
  } else {
    grDevices::pdf(
      file.path(out_dir, inst_file),
      width = inst_width,
      height = inst_height
    )
    
    inst_plot <- try(
      bibliometrix::networkPlot(
        net_inst,
        n = inst_n,
        type = inst_type,
        Title = "机构合作网络",
        size = inst_node_size,
        size.cex = inst_size_by_frequency,
        remove.isolates = inst_remove_isolates,
        labelsize = inst_labelsize,
        verbose = TRUE
      ),
      silent = TRUE
    )
    
    grDevices::dev.off()
    
    if (inherits(inst_plot, "try-error")) {
      message("机构合作网络绘制失败：", as.character(inst_plot))
    }
  }
}

# 10、绘制关键词共现网络

# 关键词共现网络参数
keywords_n <- 15
keywords_type <- "fruchterman"
keywords_labelsize <- 0.7
keywords_edgesize <- 2
keywords_remove_isolates <- TRUE
keywords_cluster <- "walktrap"

# PDF保存参数
keywords_width <- 10
keywords_height <- 8
keywords_file <- "5.关键词共现网络.pdf"

keyword_available <-
  ("DE" %in% names(M) && any(nchar(trimws(as.character(M$DE))) > 0, na.rm = TRUE)) ||
  ("ID" %in% names(M) && any(nchar(trimws(as.character(M$ID))) > 0, na.rm = TRUE))

if (!keyword_available) {
  message("跳过关键词共现网络：当前数据缺少有效的DE和ID字段。")
} else {
  net_keywords <- try(
    bibliometrix::biblioNetwork(
      M,
      analysis = "co-occurrences",
      network = "keywords",
      sep = ";"
    ),
    silent = TRUE
  )
  
  if (inherits(net_keywords, "try-error")) {
    message("关键词共现网络构建失败：", as.character(net_keywords))
  } else {
    grDevices::pdf(
      file.path(out_dir, keywords_file),
      width = keywords_width,
      height = keywords_height
    )
    
    keywords_plot <- try(
      bibliometrix::networkPlot(
        net_keywords,
        type = keywords_type,
        Title = "关键词共现网络",
        n = keywords_n,
        remove.isolates = keywords_remove_isolates,
        edgesize = keywords_edgesize,
        labelsize = keywords_labelsize,
        cluster = keywords_cluster,
        verbose = TRUE
      ),
      silent = TRUE
    )
    
    grDevices::dev.off()
    
    if (inherits(keywords_plot, "try-error")) {
      message("关键词共现网络绘制失败：", as.character(keywords_plot))
    }
  }
}

# 11、绘制文献共被引网络

# 文献共被引网络参数
ref_n <- 30
ref_type <- "kamada"
ref_labelsize <- 0.7
ref_node_size <- 3
ref_size_by_frequency <- TRUE
ref_remove_isolates <- TRUE

# PDF保存参数
ref_width <- 10
ref_height <- 8
ref_file <- "6.文献共被引网络.pdf"

if (!"CR" %in% names(M) || !any(nchar(trimws(as.character(M$CR))) > 0, na.rm = TRUE)) {
  message("跳过文献共被引网络：当前数据缺少有效的CR字段。")
} else {
  net_ref <- try(
    bibliometrix::biblioNetwork(
      M,
      analysis = "co-citation",
      network = "references",
      sep = ";"
    ),
    silent = TRUE
  )
  
  if (inherits(net_ref, "try-error")) {
    message("文献共被引网络构建失败：", as.character(net_ref))
  } else {
    grDevices::pdf(file.path(out_dir, ref_file), width = ref_width, height = ref_height)
    
    ref_plot <- try(
      bibliometrix::networkPlot(
        net_ref,
        n = ref_n,
        type = ref_type,
        Title = "文献共被引网络",
        size = ref_node_size,
        size.cex = ref_size_by_frequency,
        remove.isolates = ref_remove_isolates,
        labelsize = ref_labelsize,
        verbose = TRUE
      ),
      silent = TRUE
    )
    
    grDevices::dev.off()
    
    if (inherits(ref_plot, "try-error")) {
      message("文献共被引网络绘制失败：", as.character(ref_plot))
    }
  }
}

# 12、绘制期刊共被引网络

# 期刊共被引网络参数
source_n <- 30
source_type <- "fruchterman"
source_labelsize <- 0.8
source_node_size <- 3
source_size_by_frequency <- TRUE
source_remove_isolates <- TRUE

# PDF保存参数
source_width <- 10
source_height <- 8
source_file <- "7.期刊共被引网络.pdf"

if (!"CR" %in% names(M) || !any(nchar(trimws(as.character(M$CR))) > 0, na.rm = TRUE)) {
  message("跳过期刊共被引网络：当前数据缺少有效的CR字段。")
} else {
  M_crso <- try(
    bibliometrix::metaTagExtraction(M, Field = "CR_SO", sep = ";"),
    silent = TRUE
  )
  
  if (
    inherits(M_crso, "try-error") ||
    !"CR_SO" %in% names(M_crso) ||
    !any(nchar(trimws(as.character(M_crso$CR_SO))) > 0, na.rm = TRUE)
  ) {
    message("跳过期刊共被引网络：未成功生成有效的CR_SO字段。")
  } else {
    net_source <- try(
      bibliometrix::biblioNetwork(
        M_crso,
        analysis = "co-citation",
        network = "sources",
        sep = ";"
      ),
      silent = TRUE
    )
    
    if (inherits(net_source, "try-error")) {
      message("期刊共被引网络构建失败：", as.character(net_source))
    } else {
      grDevices::pdf(
        file.path(out_dir, source_file),
        width = source_width,
        height = source_height
      )
      
      source_plot <- try(
        bibliometrix::networkPlot(
          net_source,
          n = source_n,
          type = source_type,
          Title = "期刊共被引网络",
          size = source_node_size,
          size.cex = source_size_by_frequency,
          remove.isolates = source_remove_isolates,
          labelsize = source_labelsize,
          verbose = TRUE
        ),
        silent = TRUE
      )
      
      grDevices::dev.off()
      
      if (inherits(source_plot, "try-error")) {
        message("期刊共被引网络绘制失败：", as.character(source_plot))
      }
    }
  }
}

# 13、绘制作者共被引网络

# 作者共被引网络参数
author_cocite_n <- 30
author_cocite_type <- "fruchterman"
author_cocite_labelsize <- 0.8
author_cocite_node_size <- 3
author_cocite_size_by_frequency <- TRUE
author_cocite_remove_isolates <- TRUE

# PDF保存参数
author_cocite_width <- 10
author_cocite_height <- 8
author_cocite_file <- "8.作者共被引网络.pdf"

if (!"CR" %in% names(M) || !any(nchar(trimws(as.character(M$CR))) > 0, na.rm = TRUE)) {
  message("跳过作者共被引网络：当前数据缺少有效的CR字段。")
} else {
  M_crau <- try(
    bibliometrix::metaTagExtraction(M, Field = "CR_AU", sep = ";"),
    silent = TRUE
  )
  
  if (
    inherits(M_crau, "try-error") ||
    !"CR_AU" %in% names(M_crau) ||
    !any(nchar(trimws(as.character(M_crau$CR_AU))) > 0, na.rm = TRUE)
  ) {
    message("跳过作者共被引网络：未成功生成有效的CR_AU字段。")
  } else {
    net_author_cocite <- try(
      bibliometrix::biblioNetwork(
        M_crau,
        analysis = "co-citation",
        network = "authors",
        sep = ";"
      ),
      silent = TRUE
    )
    
    if (inherits(net_author_cocite, "try-error")) {
      message("作者共被引网络构建失败：", as.character(net_author_cocite))
    } else {
      grDevices::pdf(
        file.path(out_dir, author_cocite_file),
        width = author_cocite_width,
        height = author_cocite_height
      )
      
      author_cocite_plot <- try(
        bibliometrix::networkPlot(
          net_author_cocite,
          n = author_cocite_n,
          type = author_cocite_type,
          Title = "作者共被引网络",
          size = author_cocite_node_size,
          size.cex = author_cocite_size_by_frequency,
          remove.isolates = author_cocite_remove_isolates,
          labelsize = author_cocite_labelsize,
          verbose = TRUE
        ),
        silent = TRUE
      )
      
      grDevices::dev.off()
      
      if (inherits(author_cocite_plot, "try-error")) {
        message("作者共被引网络绘制失败：", as.character(author_cocite_plot))
      }
    }
  }
}

# 14、进行概念结构分析

# 概念结构分析字段，可设置为自动优先ID否则DE、ID或DE
conceptual_field_priority <- "自动优先ID否则DE"

# 概念结构分析参数
conceptual_method <- "MCA"
conceptual_minDegree <- 4
conceptual_clust <- 4
conceptual_labelsize <- 10
conceptual_documents <- 5
conceptual_stemming <- FALSE
conceptual_index <- 1

# PDF保存参数
conceptual_width <- 10
conceptual_height <- 8
conceptual_current_file <- "9.概念结构分析_当前.pdf"
conceptual_all_file <- "9.概念结构分析_全部.pdf"

conceptual_field <- NULL

if (conceptual_field_priority == "自动优先ID否则DE") {
  if ("ID" %in% names(M) && any(nchar(trimws(as.character(M$ID))) > 0, na.rm = TRUE)) {
    conceptual_field <- "ID"
  } else if ("DE" %in% names(M) && any(nchar(trimws(as.character(M$DE))) > 0, na.rm = TRUE)) {
    conceptual_field <- "DE"
  }
} else {
  conceptual_field <- conceptual_field_priority
}

if (
  is.null(conceptual_field) ||
  !conceptual_field %in% names(M) ||
  !any(nchar(trimws(as.character(M[[conceptual_field]]))) > 0, na.rm = TRUE)
) {
  message("跳过概念结构分析：当前数据缺少有效的ID和DE字段。")
  conceptual_result <- NULL
  conceptual_plots <- NULL
} else {
  conceptual_result <- try(
    bibliometrix::conceptualStructure(
      M,
      field = conceptual_field,
      method = conceptual_method,
      minDegree = conceptual_minDegree,
      clust = conceptual_clust,
      stemming = conceptual_stemming,
      labelsize = conceptual_labelsize,
      documents = conceptual_documents,
      graph = FALSE
    ),
    silent = TRUE
  )
  
  if (inherits(conceptual_result, "try-error")) {
    message("概念结构分析失败：", as.character(conceptual_result))
    conceptual_plots <- NULL
  } else {
    conceptual_plots <- list()
    
    if (
      !is.null(conceptual_result$graph_terms) &&
      inherits(conceptual_result$graph_terms, "ggplot")
    ) {
      conceptual_plots[[length(conceptual_plots) + 1]] <- conceptual_result$graph_terms
    }
    
    if (
      !is.null(conceptual_result$graph_documents_Contrib) &&
      inherits(conceptual_result$graph_documents_Contrib, "ggplot")
    ) {
      conceptual_plots[[length(conceptual_plots) + 1]] <-
        conceptual_result$graph_documents_Contrib
    }
    
    if (
      !is.null(conceptual_result$graph_documents_TC) &&
      inherits(conceptual_result$graph_documents_TC, "ggplot")
    ) {
      conceptual_plots[[length(conceptual_plots) + 1]] <-
        conceptual_result$graph_documents_TC
    }
  }
  
  if (!is.null(conceptual_plots) && length(conceptual_plots) > 0) {
    conceptual_index <- max(
      1,
      min(as.integer(conceptual_index), length(conceptual_plots))
    )
    
    grDevices::pdf(
      file.path(out_dir, conceptual_current_file),
      width = conceptual_width,
      height = conceptual_height
    )
    
    print(conceptual_plots[[conceptual_index]])
    grDevices::dev.off()
    
    grDevices::pdf(
      file.path(out_dir, conceptual_all_file),
      width = conceptual_width,
      height = conceptual_height,
      onefile = TRUE
    )
    
    for (i in seq_along(conceptual_plots)) {
      print(conceptual_plots[[i]])
    }
    
    grDevices::dev.off()
  } else {
    message("概念结构分析未返回可保存的图形。")
  }
}

# 15、绘制历史引文图谱

# 历史引文图谱参数
hist_n <- 20
hist_size <- 10
hist_labelsize <- 5
hist_remove_isolates <- TRUE
hist_label <- "short"

# PDF保存参数
hist_width <- 10
hist_height <- 8
hist_file <- "10.历史引文图谱.pdf"

if (!"CR" %in% names(M) || !any(nchar(trimws(as.character(M$CR))) > 0, na.rm = TRUE)) {
  message("跳过历史引文图谱：当前数据缺少有效的CR字段。")
  hist_results <- NULL
} else {
  hist_results <- try(
    bibliometrix::histNetwork(
      M,
      sep = ";",
      network = TRUE,
      verbose = FALSE
    ),
    silent = TRUE
  )
  
  if (inherits(hist_results, "try-error")) {
    message("历史引文网络构建失败：", as.character(hist_results))
    hist_results <- NULL
  } else {
    grDevices::pdf(file.path(out_dir, hist_file), width = hist_width, height = hist_height)
    
    hist_plot <- try(
      bibliometrix::histPlot(
        hist_results,
        n = hist_n,
        size = hist_size,
        labelsize = hist_labelsize,
        remove.isolates = hist_remove_isolates,
        label = hist_label,
        verbose = TRUE
      ),
      silent = TRUE
    )
    
    grDevices::dev.off()
    
    if (inherits(hist_plot, "try-error")) {
      message("历史引文图谱绘制失败：", as.character(hist_plot))
    }
  }
}

# 16、绘制并保存三域图

# 中间字段，可设置为自动优先DE否则ID、DE或ID
three_field_mid_priority <- "自动优先DE否则ID"

# 三个字段分别显示的项目数量
three_n_left <- 10
three_n_mid <- 10
three_n_right <- 10

# HTML保存参数
threefields_file <- "11.三域图.html"
threefields_selfcontained <- TRUE

three_mid_field <- NULL

if (three_field_mid_priority == "自动优先DE否则ID") {
  if ("DE" %in% names(M) && any(nchar(trimws(as.character(M$DE))) > 0, na.rm = TRUE)) {
    three_mid_field <- "DE"
  } else if ("ID" %in% names(M) && any(nchar(trimws(as.character(M$ID))) > 0, na.rm = TRUE)) {
    three_mid_field <- "ID"
  }
} else {
  three_mid_field <- three_field_mid_priority
}

threefields_available <-
  "AU" %in% names(M) &&
  "SO" %in% names(M) &&
  any(nchar(trimws(as.character(M$AU))) > 0, na.rm = TRUE) &&
  any(nchar(trimws(as.character(M$SO))) > 0, na.rm = TRUE) &&
  !is.null(three_mid_field) &&
  three_mid_field %in% names(M) &&
  any(nchar(trimws(as.character(M[[three_mid_field]]))) > 0, na.rm = TRUE)

if (!threefields_available) {
  message("跳过三域图：当前数据缺少有效的AU、SO或DE/ID字段。")
  threefields_plot <- NULL
} else {
  threefields_plot <- try(
    bibliometrix::threeFieldsPlot(
      M,
      fields = c("AU", three_mid_field, "SO"),
      n = c(three_n_left, three_n_mid, three_n_right)
    ),
    silent = TRUE
  )
  
  if (inherits(threefields_plot, "try-error")) {
    message("三域图绘制失败：", as.character(threefields_plot))
    threefields_plot <- NULL
  } else {
    save_threefields_result <- try(
      htmlwidgets::saveWidget(
        threefields_plot,
        file = file.path(out_dir, threefields_file),
        selfcontained = threefields_selfcontained
      ),
      silent = TRUE
    )
    
    if (inherits(save_threefields_result, "try-error")) {
      message("三域图HTML保存失败：", as.character(save_threefields_result))
    }
  }
}

# 17、提取高被引文献

# 提取高被引文献的数量
highly_cited_n <- 10

if (!"TC" %in% names(M)) {
  message("跳过高被引文献提取：当前数据缺少TC字段。")
  top_cited <- NULL
} else {
  TC_clean <- trimws(as.character(M$TC))
  TC_clean[TC_clean == ""] <- NA_character_
  TC_clean <- gsub(",", "", TC_clean)
  TC_num <- suppressWarnings(as.numeric(TC_clean))
  
  if (all(is.na(TC_num))) {
    message("跳过高被引文献提取：TC字段没有可用的数值型被引次数。")
    top_cited <- NULL
  } else {
    keep_cols <- intersect(c("TI", "AU", "PY", "SO", "TC"), names(M))
    top_cited <- M[!is.na(TC_num), keep_cols, drop = FALSE]
    top_cited$TC_num <- TC_num[!is.na(TC_num)]
    top_cited <- top_cited[order(-top_cited$TC_num), , drop = FALSE]
    top_cited <- head(top_cited, highly_cited_n)
    
    write.csv(
      top_cited,
      file.path(out_dir, paste0("高被引文献前", highly_cited_n, "条.csv")),
      row.names = FALSE
    )
  }
}

# 18、记录并保存本次分析参数

parameter_text <- paste0(
  "运行流程说明：\n",
  "1. 使用普通R脚本读取文献数据库导出的文件。\n",
  "2. 使用bibliometrix包的convert2df()将原始文献文件转换为标准化数据框。\n",
  "3. 使用biblioAnalysis()和summary()进行基础文献计量分析。\n",
  "4. 使用plot(results)生成并保存基础概览图。\n",
  "5. 根据PY字段绘制年发文趋势图。\n",
  "6. 使用metaTagExtraction()提取AU_CO、CR_SO和CR_AU字段。\n",
  "7. 使用biblioNetwork()构建合作、共现和共被引网络。\n",
  "8. 使用networkPlot()绘制国家、作者、机构、关键词及共被引网络。\n",
  "9. 使用conceptualStructure()进行概念结构分析。\n",
  "10. 使用histNetwork()和histPlot()绘制历史引文图谱。\n",
  "11. 使用threeFieldsPlot()绘制作者、关键词和期刊三域图，并导出HTML。\n",
  "12. 根据TC字段提取高被引文献。\n\n",
  "本次参数总结：\n",
  "- dbsource：", dbsource, "\n",
  "- format：", file_format, "\n",
  "- summary()/plot()显示前K项：", summary_k, "\n",
  "- conceptualStructure字段：",
  if (is.null(conceptual_field)) "未使用" else conceptual_field,
  "\n",
  "- conceptualStructure method：", conceptual_method, "\n",
  "- conceptualStructure minDegree：", conceptual_minDegree, "\n",
  "- conceptualStructure clust：", conceptual_clust, "\n",
  "- conceptualStructure labelsize：", conceptual_labelsize, "\n",
  "- conceptualStructure documents：", conceptual_documents, "\n",
  "- conceptualStructure stemming：", conceptual_stemming, "\n",
  "- 三域图中间字段：",
  if (is.null(three_mid_field)) "未使用" else three_mid_field,
  "\n\n",
  "写作提示词：\n",
  "1. 基于bibliometrix包对目标研究领域文献进行文献计量学分析。\n",
  "2. 从发文时间分布、作者合作、机构合作、国家合作、关键词共现和共被引关系等方面描述研究领域的发展特征。\n",
  "3. 结合概念结构分析与三域图，识别研究热点、知识结构与主题特征。\n",
  "4. 结合高被引文献与历史引文图谱，分析关键文献及知识基础。\n",
  "5. 网络图中的节点大小通常反映频次或影响力，连线反映合作、共现或共被引关系。"
)

writeLines(parameter_text, file.path(out_dir, "bibliometrix_parameters.txt"))

message("文献计量学分析完成。")
message("结果保存位置：", normalizePath(out_dir, winslash = "/", mustWork = FALSE))