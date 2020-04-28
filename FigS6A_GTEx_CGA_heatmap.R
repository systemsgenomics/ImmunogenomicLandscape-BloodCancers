
# Plot heatmaps of CGA (this study and CTdatabase) using GTEx data (Figure S6A)

library(data.table)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

# load data
data <- fread(file = "GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.txt", data.table = F)

# CGA list
genelist <- as.character(read.table("cga.txt")[,1])

# import list of CTdatabase antigens and merge to GTEX
ctdatabase <- as.character(read.table(file = "CTdatabase.txt", header = FALSE)$V1)


# heatmap of Hemap CGAs
mat <- data[data$Description%in%genelist,-c(1:2)]
rownames(mat) <- data$Description[data$Description%in%genelist]
mat <- log2(mat+0.25)
colnames(mat) <- gsub("_", " ", gsub(")", "", colnames(data[data$Description%in%genelist,-c(1:2)])))
mat <- t(mat) # transpose

ht1 <- Heatmap(mat,
              column_title = "Hemap cancer-germline antigens",
              column_title_gp = gpar(fontsize = 9),
              name = "Expression (log2 RPKM)", 
              col = colorRamp2(seq(log2(0.25), 9, length.out = 9), (brewer.pal(9, "Reds"))),
              rect_gp = gpar(col = "white", lty = 1, lwd = 1),
              show_row_names = FALSE,
              column_names_gp = gpar(fontsize = 6),
              cluster_rows = FALSE,
              cluster_columns = TRUE,
              heatmap_legend_param = list(title_gp = gpar(fontsize = 7, fontface = "bold"),
                                          labels_gp = gpar(fontsize = 7),
                                          grid_height = unit(0.2, "cm"),
                                          grid_width = unit(2, "mm")),
              width = unit(7, "cm")
)



# heatmap of CTdatabase genes

mat <- data[data$Description%in%ctdatabase,-c(1:2)]
rownames(mat) <- data$Description[data$Description%in%ctdatabase]
mat <- log2(mat+0.25)
colnames(mat) <- gsub("_", " ", gsub(")", "", colnames(data[data$Description%in%genelist,-c(1:2)])))
mat <- t(mat) # transpose

ht2 <- Heatmap(mat,
        column_title = "CTdatabase cancer-germline antigens",
        column_title_gp = gpar(fontsize = 9),
        name = "Expression (log2 RPKM)", 
        col = colorRamp2(seq(log2(0.25), 9, length.out = 9), (brewer.pal(9, "Reds"))),
        row_names_side = "right", 
        row_names_gp = gpar(fontsize = 6),
        show_column_names = FALSE,
        cluster_rows = FALSE,
        cluster_columns = TRUE,
        heatmap_legend_param = list(title_gp = gpar(fontsize = 7, fontface = "bold"),
                                    labels_gp = gpar(fontsize = 7),
                                    grid_height = unit(0.2, "cm"),
                                    grid_width = unit(2, "mm")),
        width = unit(7, "cm")
)

pdf("FigureS6A_GTEx_CGA_heatmap.pdf", height = 6, width = 10)
draw(ht1 + ht2,
     gap = unit(0.5, "cm"),
     padding = unit(c(1, 1, 0.2, 0.2), "cm"),
     heatmap_legend_side = "bottom"
)
dev.off()
