# 加载必要的包
library(vegan)
library(lme4)
library(car)
library(multcomp)
library(reshape2)

# 函数：读取数据
load_data <- function() {
  # 这里读取您的数据
  # 例如：gene_copies <- read.csv("data/processed/gene_copies.csv")
  # 返回数据列表
  return(list(
    gene_data = gene_copies,
    soil_data = soil_properties,
    plant_data = plant_measurements
  ))
}

# 函数：差异显著性检验
perform_statistical_tests <- function(data) {
  # 例如：ANOVA分析处理间差异
  aov_nifH <- aov(nifH ~ treatment + year + treatment:year, data = data$gene_data)
  summary(aov_nifH)
  
  # Tukey's HSD事后检验
  tukey_nifH <- TukeyHSD(aov_nifH, "treatment")
  
  # 相关性分析（如RDA/dbRDA）
  rda_result <- rda(nifH_community ~ Ks + PAD + EC, data = data$soil_data)
  
  # Mantel检验
  mantel_test <- mantel(gene_dist, soil_dist, method = "pearson")
  
  return(list(
    anova_results = summary(aov_nifH),
    tukey_results = tukey_nifH,
    rda_results = rda_result,
    mantel_results = mantel_test
  ))
}

# 函数：保存统计结果
save_statistical_results <- function(results, output_dir = "results/stats/") {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 保存为文本文件
  sink(paste0(output_dir, "statistical_results.txt"))
  print("ANOVA Results:")
  print(results$anova_results)
  print("\nTukey HSD Results:")
  print(results$tukey_results)
  sink()
  
  # 保存R对象
  saveRDS(results, paste0(output_dir, "statistical_results.rds"))
  
  cat("统计结果已保存至:", output_dir, "\n")
}

# 主执行流程
main <- function() {
  cat("开始统计分析...\n")
  
  # 1. 加载数据
  data <- load_data()
  
  # 2. 执行统计检验
  results <- perform_statistical_tests(data)
  
  # 3. 保存结果
  save_statistical_results(results)
  
  cat("统计分析完成!\n")
}

# 运行主函数
if (sys.nframe() == 0) {
  main()
}
3. 网络分析的R脚本示例 (03_network_analysis.R)
r
# 网络分析脚本
library(igraph)
library(SpiecEasi)  # 用于微生物网络推断
library(qgraph)
library(WGCNA)

# 函数：构建共现网络
build_cooccurrence_network <- function(otu_table, method = "sparcc", threshold = 0.3) {
  # 方法1: SparCC相关性网络
  if (method == "sparcc") {
    # 使用SpiecEasi包
    se.mb <- spiec.easi(otu_table, method='mb', lambda.min.ratio=1e-2,
                        nlambda=20, pulsar.params=list(rep.num=50))
    adj_matrix <- getRefit(se.mb)
  }
  
  # 方法2: Spearman相关网络
  else if (method == "spearman") {
    cor_matrix <- cor(t(otu_table), method = "spearman")
    adj_matrix <- ifelse(abs(cor_matrix) > threshold, cor_matrix, 0)
  }
  
  # 转换为igraph对象
  network <- graph_from_adjacency_matrix(adj_matrix, 
                                         mode = "undirected", 
                                         weighted = TRUE,
                                         diag = FALSE)
  
  return(list(
    network = network,
    adjacency_matrix = adj_matrix
  ))
}

# 函数：计算网络稳定性指标
calculate_network_stability <- function(network, n_iterations = 100) {
  # 1. 节点移除分析（随机）
  robustness_random <- numeric(n_iterations)
  for(i in 1:n_iterations) {
    # 随机移除50%节点
    nodes_to_remove <- sample(V(network), round(vcount(network) * 0.5))
    subgraph <- delete_vertices(network, nodes_to_remove)
    
    # 计算剩余网络的连通性
    robustness_random[i] <- max(components(subgraph)$csize) / vcount(network)
  }
  
  # 2. 针对性移除（按度中心性）
  degree_order <- order(degree(network), decreasing = TRUE)
  robustness_targeted <- numeric(length(degree_order))
  for(j in seq_along(degree_order)) {
    subgraph <- delete_vertices(network, degree_order[1:j])
    robustness_targeted[j] <- ifelse(vcount(subgraph) > 0, 
                                     max(components(subgraph)$csize) / vcount(network),
                                     0)
  }
  
  # 3. 计算脆弱性
  vulnerability <- 1 - mean(robustness_random)
  
  return(list(
    robustness_random = mean(robustness_random),
    robustness_targeted = mean(robustness_targeted),
    vulnerability = vulnerability
  ))
}

# 函数：计算网络拓扑属性
calculate_network_properties <- function(network) {
  properties <- list(
    # 基本属性
    n_nodes = vcount(network),
    n_edges = ecount(network),
    density = graph.density(network),
    
    # 中心性指标
    average_degree = mean(degree(network)),
    average_clustering_coefficient = transitivity(network, type = "average"),
    
    # 模块性
    modularity = modularity(cluster_fast_greedy(network)),
    
    # 连通性
    connectedness = ecount(network) / choose(vcount(network), 2),
    
    # 路径长度
    average_path_length = average.path.length(network)
  )
  
  return(properties)
}

# 主函数
main <- function() {
  cat("开始网络分析...\n")
  
  # 1. 加载OTU数据
  # otu_data <- read.csv("data/processed/otu_table_normalized.csv", row.names = 1)
  
  # 2. 为每个处理构建网络
  treatments <- c("CK", "FMF", "MOF")
  networks <- list()
  stability_results <- list()
  
  for(trt in treatments) {
    cat("处理", trt, "...\n")
    
    # 提取该处理的OTU数据
    # trt_data <- subset(otu_data, treatment == trt)
    
    # 构建网络
    network_result <- build_cooccurrence_network(trt_data, method = "spearman")
    networks[[trt]] <- network_result$network
    
    # 计算网络属性
    properties <- calculate_network_properties(network_result$network)
    
    # 计算稳定性
    stability <- calculate_network_stability(network_result$network)
    
    # 保存结果
    stability_results[[trt]] <- list(
      properties = properties,
      stability = stability
    )
  }
  
  # 3. 保存所有结果
  output_dir <- "results/networks/"
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  saveRDS(networks, paste0(output_dir, "cooccurrence_networks.rds"))
  saveRDS(stability_results, paste0(output_dir, "network_stability_results.rds"))
  
  # 生成汇总表格（类似您论文中的Table 4）
  summary_table <- do.call(rbind, lapply(stability_results, function(x) {
    data.frame(
      Nodes = x$properties$n_nodes,
      Edges = x$properties$n_edges,
      AvgDegree = x$properties$average_degree,
      ClusteringCoeff = x$properties$average_clustering_coefficient,
      Modularity = x$properties$modularity,
      RobustnessRandom = x$stability$robustness_random,
      Vulnerability = x$stability$vulnerability
    )
  }))
  
  write.csv(summary_table, paste0(output_dir, "network_summary_table.csv"))
  
  cat("网络分析完成! 结果保存在:", output_dir, "\n")
}

if (sys.nframe() == 0) {
  main()
}
4. 图形生成的R脚本示例 (04_figure_generation.R)
r
# 图形生成脚本
library(ggplot2)
library(ggpubr)
library(vegan)
library(RColorBrewer)
library(patchwork)

# 主题设置
theme_custom <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
    legend.position = "right",
    legend.box.background = element_rect(colour = "black"),
    axis.text = element_text(color = "black"),
    strip.background = element_rect(fill = "gray95")
  )

# 函数：生成图1 (NH3挥发和N吸收)
generate_figure1 <- function(data, output_path = "results/figures/figure1.png") {
  # 准备数据
  plot_data <- data.frame(
    Treatment = rep(c("CK", "FMF", "MOF"), each = 2),
    Year = rep(c("2019", "2020"), 3),
    NH3_emission = c(100, 95, 121, 115, 90, 85),  # 示例数据
    N_uptake = c(50, 55, 71.5, 105, 90, 99)        # 示例数据
  )
  
  # A: NH3挥发
  p1 <- ggplot(plot_data, aes(x = Treatment, y = NH3_emission, fill = Year)) +
    geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
    scale_fill_brewer(palette = "Set2") +
    labs(x = "Treatment", y = expression("NH"[3]*" emission (kg/ha)"),
         title = "A) NH₃ volatilization") +
    theme_custom +
    theme(legend.position = c(0.1, 0.85))
  
  # B: N吸收
  p2 <- ggplot(plot_data, aes(x = Treatment, y = N_uptake, fill = Year)) +
    geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
    scale_fill_brewer(palette = "Set2") +
    labs(x = "Treatment", y = "Nitrogen uptake (kg/ha)",
         title = "B) Nitrogen uptake") +
    theme_custom +
    theme(legend.position = "none")
  
  # C: 相关性图（示例）
  cor_data <- data.frame(
    Ks = runif(20, 0.5, 2),
    PAD = runif(20, 20, 50),
    NH3 = runif(20, 80, 120),
    Nup = runif(20, 50, 100)
  )
  
  p3 <- ggplot(cor_data, aes(x = Ks, y = NH3)) +
    geom_point(size = 3, alpha = 0.7, color = "steelblue") +
    geom_smooth(method = "lm", se = FALSE, color = "red") +
    labs(x = "Saturated hydraulic conductivity (Ks)",
         y = expression("NH"[3]*" emission"),
         title = "C) Relationship with soil properties") +
    theme_custom
  
  # 组合图形
  figure1 <- (p1 | p2) / p3 +
    plot_layout(heights = c(1, 1.2)) +
    plot_annotation(tag_levels = 'A')
  
  # 保存图形
  ggsave(output_path, figure1, width = 10, height = 8, dpi = 300)
  cat("图1已保存至:", output_path, "\n")
  
  return(figure1)
}

# 函数：生成网络图（图5-6）
generate_network_figures <- function(networks, output_dir = "results/figures/") {
  # 为每个处理生成网络图
  for(trt in names(networks)) {
    network <- networks[[trt]]
    
    # 设置节点颜色（按模块）
    communities <- cluster_louvain(network)
    V(network)$color <- communities$membership
    
    # 设置节点大小（按度中心性）
    V(network)$size <- sqrt(degree(network)) * 2
    
    # 设置边宽度（按权重）
    E(network)$width <- abs(E(network)$weight) * 3
    
    # 绘图
    png(paste0(output_dir, "network_", trt, ".png"), 
        width = 800, height = 800, res = 150)
    
    plot(communities, network,
         vertex.label = NA,
         vertex.frame.color = "white",
         layout = layout_with_fr,
         main = paste("Co-occurrence network -", trt))
    
    dev.off()
  }
  
  cat("网络图已保存至:", output_dir, "\n")
}

# 函数：生成概念模型图（图8）
generate_conceptual_model <- function(output_path = "results/figures/figure8_conceptual_model.png") {
  # 使用ggplot2创建简单概念图
  # 这里可以替换为更专业的绘图工具如Inkscape或DiagrammeR
  
  library(ggplot2)
  
  # 创建基础画布
  p <- ggplot() +
    xlim(0, 10) + ylim(0, 10) +
    theme_void()
  
  # 添加文字和箭头（简化版）
  p <- p +
    annotate("text", x = 5, y = 9, label = "MOF Application", size = 6, fontface = "bold") +
    annotate("segment", x = 5, xend = 5, y = 8.5, yend = 7, 
             arrow = arrow(type = "closed", length = unit(0.3, "cm"))) +
    annotate("text", x = 3, y = 6, label = "Improved Soil\nStructure", size = 4) +
    annotate("text", x = 7, y = 6, label = "Enriched nifH\nCommunity", size = 4) +
    annotate("segment", x = 3, xend = 2, y = 5.5, yend = 4, 
             arrow = arrow(type = "closed", length = unit(0.2, "cm"))) +
    annotate("segment", x = 7, xend = 8, y = 5.5, yend = 4, 
             arrow = arrow(type = "closed", length = unit(0.2, "cm"))) +
    annotate("text", x = 1.5, y = 3.5, label = "Reduced NH₃\nVolatilization", size = 4) +
    annotate("text", x = 8.5, y = 3.5, label = "Increased N\nUptake", size = 4) +
    annotate("text", x = 5, y = 1, label = "Enhanced Nitrogen Use Efficiency", 
             size = 5, fontface = "bold", color = "darkgreen")
  
  ggsave(output_path, p, width = 8, height = 6, dpi = 300)
  cat("概念模型图已保存至:", output_path, "\n")
}

# 主函数
main <- function() {
  cat("开始生成图形...\n")
  
  # 创建输出目录
  dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)
  
  # 加载数据
  # data <- load_data()
  # networks <- readRDS("results/networks/cooccurrence_networks.rds")
  
  # 生成所有图形
  # fig1 <- generate_figure1(data)
  # generate_network_figures(networks)
  # generate_conceptual_model()
  
  cat("图形生成完成!\n")
}

if (sys.nframe() == 0) {
  main()
}