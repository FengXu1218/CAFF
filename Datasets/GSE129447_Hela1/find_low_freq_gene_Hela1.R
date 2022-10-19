# 探索是否存在有些基因是以更低的频率在震荡表达 
source("functions.R")
library(tidyverse)

load("./Hela1/cpm/exp_cpm_nolog.Rdata")

Hela1_cpm_nolog <- exp_cpm_nolog[, intersect(colnames(data_cc_new), colnames(exp_cpm_nolog))]
##################### 第一步 ###################
# 先通过自相关检验找出p值接近1的基因（排除掉低表达对结果的）
# 只看一阶：
# 先排除低表达基因：
# p_value_filter <- p_value[!(rownames(p_value) %in% low_exp_gene), ]
# 
# # 选择p值大于0.9的基因，作为备选低频率基因
# length(which(p_value_filter[,1] > 0.9))
# 
# low_freq_gene_tmp <- rownames(p_value_filter)[which(p_value_filter[,1] > 0.9)]

#################### 第二步: 寻找低频率基因 ###############
source("functions.R")
############ 方法一：归一化方法采用除总平均数 ############
# 调用find_low_freq_gene函数，计算low_freq_gene_tmp中基因的方差矩阵
res_var1 <- find_low_freq_gene2(gene_list = rownames(Hela1_cpm_nolog), 
                                cpm_data = Hela1_cpm_nolog,
                                pseudo_order = pseudo_order, 
                                n = 10,
                                normalize_method = "all_avg")

# 提取方差矩阵和自回归检验p值矩阵
var_mat1 <- res_var1$var_mat
p_value_var1 <- res_var1$p_value_var

# 随意挑一个显著的基因预览：
# plot(1:28, var_mat1["GINS2",])

# 对p值进行fdr矫正：
p.adj_var1 <- apply(p_value_var1, 2, function(x) p.adjust(x, method = "fdr"))

# 查看显著的个数
length(which(p_value_var1[,3]<0.05))
length(which(p.adj_var1[,3]<0.05))

# # 预览：
# group_number <- ncol(Hela1_cpm_nolog) %/% 10
# 
# # 对矫正后依然显著的基因方差进行预览：
# plot(1:group_number, var_mat[which(p.adj_var[,3]<0.05)[1],])
# plot(1:group_number, var_mat[which(p.adj_var[,3]<0.05)[2],])
# plot(1:group_number, var_mat[which(p.adj_var[,3]<0.05)[3],])

# 保存鉴定出了低频率基因：
low_freq_genes1 <- rownames(Hela1_cpm_nolog)[which(p.adj_var1[,3]<0.05)]

write.csv(low_freq_genes1, row.names = F,
          "./Hela1/periodic_analysis_of_variance/Hela1_low_freq_genes1.csv")

## 可视化优化
# 主要对矫正p值<0.05的基因方差进行可视化
library(ggplot2)
library(tidyverse)
source("functions.R")

res_list1 <- scatter_plot(var_mat1, low_freq_genes1, n_col = 4)

p <- res_list1$plot

ggsave("./Hela1/plots/hela1_variance_of_gene1.pdf", plot = p, height = 10, width = 8)

save(p_value_var1, p.adj_var1, var_mat1, low_freq_genes1, res_list1,
     file = "./Hela1/periodic_analysis_of_variance/var_Hela1_method1.Rdata")


############ 方法二：归一化方法采用除总分组平均数 ############
res_var2 <- find_low_freq_gene2(gene_list = rownames(Hela1_cpm_nolog), 
                                cpm_data = Hela1_cpm_nolog,
                                pseudo_order = pseudo_order, n = 10,
                                normalize_method = "group_avg")

# 提取方差矩阵和自回归检验p值矩阵
var_mat2 <- res_var2$var_mat
p_value_var2 <- res_var2$p_value_var

# 对p值进行fdr矫正：
p.adj_var2 <- apply(p_value_var2, 2, function(x) p.adjust(x, method = "fdr"))

# 查看显著的个数
length(which(p_value_var2[,3]<0.05))
length(which(p.adj_var2[,3]<0.05))

# 保存鉴定出了低频率基因：
low_freq_genes2 <- rownames(Hela1_cpm_nolog)[which(p.adj_var2[,3]<0.05)]

write.csv(low_freq_genes2, row.names = F,
          "./Hela1/periodic_analysis_of_variance/Hela1_low_freq_genes2.csv")

#################### 第三步：可视化优化 ###############
# 主要对矫正p值<0.05的基因方差进行可视化
library(ggplot2)
library(tidyverse)
source("functions.R")

res_list2 <- scatter_plot(var_mat2, low_freq_genes2, n_col = 4)

p <- res_list2$plot

ggsave("./Hela1/plots/hela1_variance_of_gene2.pdf", plot = p, height = 4, width = 8)

save(p_value_var2, p.adj_var2, var_mat2, low_freq_genes2, res_list2,
     file = "./Hela1/periodic_analysis_of_variance/var_Hela1_method2.Rdata")




