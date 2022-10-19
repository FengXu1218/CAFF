## Data preprocessing ------------------------
library(edgeR)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(pracma)
source("./R/functions.R")
library(data.table)

dir <- "./DataSets/GSE129447_Hela1/"
Hela1 <- read.table(paste0(dir, "rawdata/GSM3713084_HeLa_1.txt"), header = T, row.names = 1)

usethis::use_data(Hela1)
# Check the sequencing depth distribution:
quantile(colSums(Hela1), probs = c(0,0.1, 0.25,0.5,0.75,0.9, 1))

# Cells with sequencing depth less than 10000 were removed:
colids <- which(colSums(Hela1) > 10000)

# Genes with 0 expression in more than 70% of all cells were removed:
rowids <- which(apply(Hela1, 1, function(x) length(which(x == 0))/length(x) <= 0.7))

exp <- Hela1[rowids, colids]

# cpm transformation：
exp_cpm <- cpm(exp, log = T) %>% round(2)
exp_cpm_nolog <- cpm(exp, log = F) %>% round(2)

usethis::use_data(exp, exp_cpm, exp_cpm_nolog)

# extracte the expression matrix of cell cycle related genes:
load("cc_genes/cc_genes.Rdata")

cc_genes <- intersect(seurat_cc_genes, rownames(exp_cpm))
data_cc <- exp_cpm[cc_genes,]

# the pca_plot function to see how the data is distributed on pc1 and pc2:
p_x <- pca_plot(data_cc, time.iqr = 1)
p_x$plot1
p_x$plot2

# remove_pc_outliers
final_index <- p_x$final_index
data_cc_new <- data_cc[, final_index]

usethis::use_data(data_cc, data_cc_new)

## build_pseudo_order --------------------------------
data(data_cc_new)
## Our method ==============
pseudo_order_list <- build_pseudo_order(data_cc_new, method = "Default")
pseudo_order <- pseudo_order_list$pseudo_order
pseudo_order_rank <- pseudo_order_list$pseudo_order_rank

## tricycle method ==============
tricycle_order_rank <- build_pseudo_order(data_cc_new, method = "tricycle")


## reCAT method ==============
data_cc_cyclebase <- exp_cpm[intersect(cyclebase3.0_genes, rownames(exp_cpm)),
                             intersect(colnames(data_cc_new), colnames(exp_cpm))]
dim(data_cc_cyclebase)

reCAT_order <- build_pseudo_order(data_cc_cyclebase, method = "reCAT")
reCAT_order_rank <- reCAT_order[[2]]

## Comparison above the three groups ==============
rank_data <- data.frame(pseudo_order_rank = pseudo_order_rank,
                         tricycle_order_rank = tricycle_order_rank,
                         reCAT_order_rank = reCAT_order_rank)
rank_data <- rank_data[order(rank_data$pseudo_order_rank),]
data_cc_ordered <- data_cc_new[, pseudo_order]

p_list <- cor_scatter_plot(rank_data)
library(cowplot)

p <- plot_grid(plotlist = p_list, ncol = 3)
p

usethis::use_data(rank_data, pseudo_order, data_cc_ordered, overwrite = T)

## 傅里叶滤波 ---------------------
# 傅里叶变换
source("fftf.R")

# 单组:
# X_fftf <- data_cc_new[, pseudo_order]
#
# for (i in 1:nrow(data_cc_new)) {
#   res <- fftf(1:281, exp(as.numeric(data_cc_new[i,pseudo_order])), 1/(281/3))
#   X_fftf[i,] <- res[[1]]
#   print(i)
# }

# 担心在单一周期数据外曲线拟合差，人工复制制造三周期数据并分析
data_cc_bind <- cbind(data_cc_new[,pseudo_order], data_cc_new[,pseudo_order], data_cc_new[,pseudo_order])
X_fftf2 <- data_cc_bind

for (i in 1:nrow(data_cc_bind)) {
  res <- fftf(1:ncol(data_cc_bind), exp(as.numeric(data_cc_bind[i,])), 1/(ncol(data_cc_bind)/6), plot = F)
  X_fftf2[i,] <- res[[1]]
  print(i)
}

write.csv(X_fftf2, "./Hela1/X_FFTF_data/Hela1_X_fftf2.csv")

################# 对傅里叶滤波后的结果可视化 #########################
# 读取滤波后的数据：
X_fftf2 <- read.csv("./Hela1/X_FFTF_data/Hela1_X_fftf2.csv", header = T, row.names = 1)
data_cc_new <- read.csv("./Hela1/filtered_data/Hela1_cc_cpm_filter.csv", row.names = 1)

# 对傅里叶滤波后的结果重新排序，并绘制热图：
p_ind <- heatmap_plot(X_fftf2, name = "Hela1")
p1 <- p_ind$plot

pdf(paste0("./Hela1/plots/", "Hela1", "_heatmap.pdf"), height = 7, width = 5)
p1
dev.off()


