library(Seurat)
dir_path <- "./cc_genes/"

# cc.genes in Seurat ==============
seurat_cc_genes <- c(Seurat::cc.genes.updated.2019$s.genes,
                     Seurat::cc.genes.updated.2019$g2m.genes)

# cc.genes in GO0007049 ==============
library(org.Hs.eg.db)
db.o <- "org.Hs.eg.db"
id.type <- "SYMBOL"
GO0007049_genes <- unique(AnnotationDbi::select(get(db.o), keytype="GOALL", keys="GO:0007049", 
                                                columns=id.type)[, id.type])

# cc.genes in cyclebase3.0 ==============
cyclebase3.0 <- read.csv(paste0(dir_path, "cyclebase3.0_genes_phase.csv"), header = T)
cyclebase3.0_genes <- unique(cyclebase3.0$Genename)

save(seurat_cc_genes, GO0007049_genes, cyclebase3.0_genes, file = paste0(dir_path, "cc_genes.Rdata"))

# Venn-plot
library(VennDiagram)
library(RColorBrewer)

if (!dir.exists(paste0(dir_path, "plots"))) {
  dir.create(paste0(dir_path, "plots"))
}

if (T) {
  venn.diagram(
    x = list(seurat_cc_genes, GO0007049_genes, cyclebase3.0_genes),
    category.names = c("seurat" , "GO0007049 " , "cyclebase3.0"),
    filename = paste0(dir_path, "plots/Venn_of_cc_gene.png"),
    output=TRUE,
    
    # output：
    imagetype="png" ,
    height = 1000 , 
    width = 1000 , 
    resolution = 300,
    compression = "lzw",
    
    # cycle：
    lwd = 2, 
    lty = "blank",  
    fill = brewer.pal(3, "Pastel2"), 
    
    # font：
    cex = .5, 
    fontface = "bold",
    fontfamily = "sans",
    
    # name of sets:
    cat.cex = 0.6,
    cat.fontface = "bold",
    cat.default.pos = "outer",  
    cat.pos = c(-27, 27, 180),  
    cat.dist = c(0.055, 0.055, 0.055), 
    cat.fontfamily = "sans",
    rotation = 1  
  )
}

# Upset-Plot
library(UpSetR)        

upset_list <- list(seurat_cc_genes, GO0007049_genes, cyclebase3.0_genes)  
names(upset_list) <- c("seurat", "GO0007049", "cyclebase3.0")  

# Plot
pdf(paste0(dir_path, "plots/upset.pdf"), height = 7, width = 10)
upset(fromList(upset_list), order.by = "freq")
dev.off()
