cc_gene_test <- read.table("Hela1/cc_gene/cc_gene_test.txt")

intersect(low_freq_genes1, cc_gene_test[,1])

intersect(low_freq_genes2, cc_gene_test[,1])

p1 <- scatter_plot(Hela1_cpm_nolog[, pseudo_order],
             gene_list = low_freq_genes2,
             y_lab = "Expression",n_col = 5)


p2 <- scatter_plot(var_mat2,
             gene_list = low_freq_genes2,
             y_lab = "Variance",n_col = 5)

library(patchwork)
p <- p1$plot/p2$plot
ggsave("p2.pdf", plot = p, height = 8, width = 18)
