################# 细胞周期相关基因识别和比较 ########################
load("./Hela1/cpm/exp_cpm_nolog.Rdata")

# 基于重构的伪时序进行Ljung-Box检验
Hela1_cpm <- exp_cpm[,intersect(colnames(data_cc_new), colnames(exp_cpm))]
res <- list()
p_value <- matrix(NA, nrow = nrow(Hela1_cpm), ncol = 3)
lags <- c(1,2,5)

# 分别计算1、2、5阶滞后回归结果：
for (j in 1:length(lags)) {
  print(j)
  for (i in 1:nrow(Hela1_cpm)) {
    res[[i]] <- Box.test(2^as.numeric(Hela1_cpm[i,pseudo_order]), type = "Ljung", lag = lags[j])
    p_value[i, j] <- res[[i]][["p.value"]]
  }
}

length(which(p_value[,1] < 0.05))
hist(p_value[,3])
# names(p_value) <- rownames(Hela1_cpm)
# write.csv(p_value, "p_value_Hela1.csv")

p.adj <- apply(p_value, 2, function(x) p.adjust(x, method = "fdr"))
length(which(p.adj[,3] < 0.05))

cc_gene_test <- rownames(Hela1_cpm)[which(p.adj[,3] < 0.05)]

write.table(cc_gene_test, "./Hela1/cc_gene/cc_gene_test.txt", row.names = F, quote = F)

# 观察基因表达量对检验结果的影响:
# hist(rowSums(Hela1_cpm[which(p_value < 0.05),]))
# hist(rowSums(Hela1_cpm))

# 与三组细胞周期基因集比较 ----------------
load("cc_genes.Rdata")

# 与Seurat包中的G2M和S期比较：
length(intersect(seurat_cc_genes, cc_gene_test))

# 与GO0007049通路比较：
length(intersect(GO0007049_genes, cc_gene_test))

# 与cyclebase 3.0比较：
length(intersect(cyclebase3.0_genes, cc_gene_test))

# 四者韦恩图：
library(VennDiagram)

venn.diagram(
  x = list(cc_gene_test, seurat_cc_genes,
           GO0007049_genes, cyclebase3.0_genes),
  category.names = c("find in Hela1" , "Seurat" , "GO0007049", "cyclebase3.0"),
  filename = "./Hela1/plots/cc_genes_venn.png",
  
  # 设置输出：
  imagetype="png" ,
  height = 1000 , 
  width = 1000 , 
  resolution = 300,
  compression = "lzw",
  
  # 圆圈属性：
  col = "white",  # 描边颜色
  lty = 1, # 虚线形式:1,2,3,4,5可选
  lwd = 1,  # 粗细
  fill = c("#ffd7d8", "#d8f2e7", "#d9e7f2", "#eadff0"),  # 填充颜色；
  alpha = 0.90, # 透明度
  
  # 标签属性：
  label.col = "black",
  cex = .5, # 字体大小
  fontfamily = "serif",
  fontface = "bold",
  
  # 集合名称属性：
  cat.col = c("#cb6274", "#7ba498", "#687d94", "#81668b"),
  cat.cex = .6,
  cat.fontfamily = "serif"
)






