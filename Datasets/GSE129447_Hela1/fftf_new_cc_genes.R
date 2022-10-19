################ 数据预处理 ################
library(dplyr)
library(edgeR)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(pracma)
source("functions.R")
library(data.table)

load("./Hela1/cpm/exp_cpm_nolog.Rdata")
Hela1_cc_gene_test <- read.csv("./Hela1/cc_gene/cc_gene_test.txt")[,1]
Hela2_cc_gene_test <- read.csv("./Hela2/cc_gene/cc_gene_test.txt")[,1]
# GSE184542_TNBC_cc_gene_test <- read.csv("./GSE184542_TNBC/cc_gene/cc_gene_test.txt")[,1]

cc_gene_test <- Reduce(intersect, list(Hela1_cc_gene_test, Hela2_cc_gene_test#,
                                       #GSE184542_TNBC_cc_gene_test
                                       ))

################# 第一轮pca #########################
# 调用pca_plot函数查看数据在pc1和pc2上的分布情况
data_cc <- exp_cpm[cc_gene_test, ]
p_x <- pca_plot(data_cc)
p_x$plot

# plot(p_x$x[,1], p_x$x[,2])
sort(p_x$x[,2], decreasing = T)[1:2]
index_pc2 <- which(p_x$x[,2] < 8)

data_cc_new <- data_cc[,index_pc2]

################# 伪时序构建 ########################
res_order_cc_gene_test <- build_pesudo_order(data_cc_new)
pseudo_order <- res_order_cc_gene_test$pseudo_order

################# 傅里叶滤波 ########################
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

write.csv(X_fftf2, "./Hela1/X_FFTF_data/Hela1_X_fftf2_new_cc_genes.csv")

################# 对傅里叶滤波后的结果可视化 #########################
# 对傅里叶滤波后的结果重新排序，并绘制热图：
p_ind <- heatmap_plot(X_fftf2, name = "Hela1")
p1 <- p_ind$plot

pdf(paste0("./Hela1/plots/", "Hela1_new_cc_genes", "_heatmap.pdf"), height = 12, width = 5)
p1
dev.off()


