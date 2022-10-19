library(Seurat)

# Read in the expression matrix The first row is a header row, the first column is rownames
exp.mat <- read.csv("./Hela1/cpm/Hela1_cpm.csv", row.names = 1)

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Create our Seurat object and complete the initalization steps
marrow <- CreateSeuratObject(counts = exp.mat)
# marrow <- NormalizeData(marrow)
# marrow <- FindVariableFeatures(marrow, selection.method = "vst")
# marrow <- ScaleData(marrow, features = rownames(marrow))

# marrow <- RunPCA(marrow, features = VariableFeatures(marrow), ndims.print = 6:10, nfeatures.print = 10)
# DimHeatmap(marrow, dims = c(8, 10))

marrow <- CellCycleScoring(marrow, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)

# view cell cycle scores and phase assignments
head(marrow[[]])


library(ggplot2)
library(tidyverse)

data <- data.frame(G1S_score = marrow$S.Score,
           G2M_score = marrow$G2M.Score,
           Phase = marrow$Phase)

ggplot(data)+
  geom_point(aes(G1S_score, G2M_score, color = Phase))+
  theme_bw()


