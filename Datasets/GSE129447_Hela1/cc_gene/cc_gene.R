library(xlsx)

data1 <- read.xlsx("NIHMS687993-supplement-supp_data_2.xlsx", sheetIndex = 1)

data2 <- read.xlsx("NIHMS687993-supplement-supp_data_2.xlsx", sheetIndex = 2)

h_g <- data1$human.gene
h_g <- na.omit(gsub(" ", "", h_g))
cc_gene <- as.character(unlist(data2))

cc_gene <- na.omit(gsub(" ", "", cc_gene))

length(intersect(h_g, cc_gene))

cc_gene_test <- read.csv("./table/cc_gene_test.txt")[,1]

length(intersect(cc_gene, cc_gene_test))

length(intersect(cc_gene_test, h_g))

library(VennDiagram)
venn.diagram(
  x = list(cc_gene, cc_gene_test, h_g),
  category.names = c("from 2002", "cc gene test" , "from Cell"),
  filename = "./plots/cc_genes_venn2.png",
  
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
  fill = c("#ffd7d8", "#d8f2e7", "#d9e7f2"),  # 填充颜色；
  alpha = 0.90, # 透明度
  
  # 标签属性：
  label.col = "black",
  cex = .5, # 字体大小
  fontfamily = "serif",
  fontface = "bold",
  
  # 集合名称属性：
  cat.col = c("#cb6274", "#7ba498", "#687d94"),
  cat.cex = .6,
  cat.fontfamily = "serif",
  #cat.pos = c(0,-30,-130),
  cat.dist = 0.18  # 外部多少距离
)
