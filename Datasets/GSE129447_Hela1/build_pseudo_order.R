library(BAMBI)

# 读取过滤后的cpm数据 -----------------
Hela1_cpm <- read.csv("./Hela1/cpm/Hela1_cpm.csv", row.names = 1)
cc_gene_test <- read.csv("./table/cc_gene_test.txt")[,1]

load("cc_genes.Rdata")

# 根据各个基因集取子集 -----------------
Seurat_cc_cpm <- Hela1_cpm[intersect(rownames(Hela1_cpm),seurat_cc_genes),]
GO0007049_cc_cpm <- as.matrix(Hela1_cpm[intersect(rownames(Hela1_cpm),GO0007049_genes),])
# 选前500个高变基因：
GO0007049_cc_cpm <- GO0007049_cc_cpm[order(apply(GO0007049_cc_cpm, 1, var))[1:500],]
cyclebase3.0_cc_cpm <- Hela1_cpm[intersect(rownames(Hela1_cpm),cyclebase3.0_genes),]
cc_gene_test_cpm <- Hela1_cpm[cc_gene_test,]

# Seurat基因子集 ==============
source("functions.R")
p_x <- pca_plot(Seurat_cc_cpm)
p_x$plot

# sort(p_x$x[,2], decreasing = T)[1:2]
# index_pc2 <- which(p_x$x[,2] < 8)

Seurat_cc_cpm_new <- Seurat_cc_cpm#[,index_pc2]

# 构建伪时序：
res_order_Seurat <- build_pseudo_order(Seurat_cc_cpm_new)

# pca作图：
data <- as.data.frame(res_order_Seurat$score)

p_seurat <- pca_scatter_plot(data, res_order_Seurat$sdev,
                             title = "PCA plot of Seurat cc genes in Hela1",
                             sub_title = "")
p_seurat

# GO0007049基因子集 ==============
p_x <- pca_plot(GO0007049_cc_cpm)
p_x$plot

# index_pc1 <- which(p_x$x[,1] > -10 & p_x$x[,1] < 10)
# index_pc2 <- which(p_x$x[,2] < 10)
# index <- intersect(index_pc1, index_pc2)

GO0007049_cc_cpm_new <- GO0007049_cc_cpm#[,index]

res_order_GO0007049 <- build_pseudo_order(GO0007049_cc_cpm_new)

# pca作图：
data <- as.data.frame(res_order_GO0007049$score)

p_GO0007049 <- pca_scatter_plot(data, res_order_GO0007049$sdev,
                             title = "PCA plot of GO0007049 cc genes in Hela1",
                             sub_title = "")
p_GO0007049

# cyclebase3.0基因子集 ==============
p_x <- pca_plot(cyclebase3.0_cc_cpm)
p_x$plot

# sort(p_x[["x"]][,1])[1:6]
# index_pc1 <- which(p_x$x[,1] > -7)
# index_pc2 <- which(p_x$x[,2] < 10)
# index <- intersect(index_pc1, index_pc2)

cyclebase3.0_cc_cpm_new <- cyclebase3.0_cc_cpm#[,index]

res_order_cyclebase3.0 <- build_pseudo_order(cyclebase3.0_cc_cpm_new)

# pca作图：
data <- as.data.frame(res_order_cyclebase3.0$score)

p_cyclebase3.0 <- pca_scatter_plot(data, res_order_cyclebase3.0$sdev,
                                title = "PCA plot of cyclebase3.0 cc genes in Hela1",
                                sub_title = "")
p_cyclebase3.0


# cc_gene_test基因子集 ===========
source("functions.R")
p_x <- pca_plot(cc_gene_test_cpm)
p_x$plot

# sort(p_x$x[,2], decreasing = T)[1:2]
# index_pc2 <- which(p_x$x[,2] < 9)

cc_gene_test_cpm_new <- cc_gene_test_cpm#[,index_pc2]

# 构建伪时序：
res_order_cc_gene_test <- build_pseudo_order(cc_gene_test_cpm_new)

# pca作图：
data <- as.data.frame(res_order_Seurat$score)

p_cc_gene_test <- pca_scatter_plot(data, res_order_cc_gene_test$sdev,
                             title = "PCA plot of finded cc genes in Hela1",
                             sub_title = "")
p_cc_gene_test

# 拼图：
library(cowplot)

p <- plot_grid(p_seurat, p_GO0007049, 
               p_cyclebase3.0, p_cc_gene_test, ncol = 2)

ggsave(plot = p, filename = "./Hela1/plots/pca_plot_cc_genes.pdf", height = 9, width = 9)

# 三组比较+可视化 -------------------
source("functions.R")

Theta_data <- data.frame(
  THETA_of_Seurat = res_order_Seurat$THETA,
  THETA_of_cc_gene_test = res_order_cc_gene_test$THETA,
  THETA_of_GO0007049 = res_order_GO0007049$THETA,
  THETA_of_cyclebase3.0 = res_order_cyclebase3.0$THETA
)

# 调用cor_scatter_plot绘制最大化相关性散点图：
p_list <- circ_cor_scatter_plot(Theta_data)

library(cowplot)

p <- plot_grid(plotlist = p_list, ncol = 3)
p
ggsave("./Hela1/plots/cor_of_THETA.pdf", plot = p, height = 8, width = 12)

