
# Plot Figure 3D oncoprint of genetic alterations correlated with cytolytic score in TCGA AML

library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(grid)

# load data
dufva=get(load("DUFVA_TCGA_AML_FM_meth.Rdata"))

# select samples with RNA-seq data
data <- dufva[,is.na(dufva["N:SAMP:CytolyticScore",])==FALSE]

# order samples by cytolytic score
data <- data[,order(data["N:SAMP:CytolyticScore",], decreasing = TRUE)]

# load correlations
cor <- read.table("TableS3_TCGA_AML_cytScore_correlations.tsv", sep = "\t", header = TRUE)

# create matrix with selected alterations

# top 4 mutations
gnab_top4 <- cor[order(cor$p),]
gnab_top4 <- gnab_top4[grepl("GNAB", gnab_top4$featureB),]
gnab_top4 <- as.character(gnab_top4[c(1:4), "featureB"])

# top 4 Gistic copy number alterations
cnvr_top4 <- cor[order(cor$p),]
cnvr_top4 <- cnvr_top4[grepl("Gistic", cnvr_top4$featureB),]
cnvr_top4 <- as.character(cnvr_top4[c(1:4), "featureB"])

feat <- c(gnab_top4, cnvr_top4)

mat <- data.matrix(data[feat,])


## change codes for alterations
mat[abs(mat)<0.3] = 0
mat[abs(mat)>0.3] = 1
mat[grepl("_amp", rownames(mat)),][mat[grepl("_amp", rownames(mat)),]==1] <- "Amplification"
mat[grepl("_del|Arm", rownames(mat)),][mat[grepl("_del|Arm", rownames(mat)),]==1] <- "Deletion"
mat[grepl("GNAB", rownames(mat)),][mat[grepl("GNAB", rownames(mat)),]==1] <- "Nonsynonymous mutation"
mat[is.na(mat)] <- "NA"

# change row names
rownames(mat) <- gsub(":.*", "", gsub("^.......", "", rownames(mat)))

## create heatmap annotations
annot <- data.frame(cluster = factor(as.character(data["C:SAMP:cancermap_cluster",])), FAB = as.factor(gsub("_Undifferentiated", "", gsub("Not_Classified", "Unclassified", as.character(data["C:CLIN:leukemia_french_american_british_morphology_code:::::",])))))

frag <- data["N:SAMP:GENOME_FRAGMENTATION_RATE",]
frag[is.na(frag)] <- 0
frag <- as.numeric(frag)
frag <- frag*100

ha1 <- HeatmapAnnotation(Cytolytic_score = anno_barplot(as.numeric(data["N:SAMP:CytolyticScore",]), 
                                                  bar_width = 0.75, 
                                                  border = FALSE, 
                                                  axis = TRUE,
                                                  axis_param = list(gp = gpar(fontsize = 5, lwd = 0.5)),
                                                  gp = gpar(col = NA, fill = "grey20"),
                                                  ylim = c(1.5,10.5)),
                         show_annotation_name = F,
                         height = unit(0.75, "cm"))


ha2 = HeatmapAnnotation(Mutations = anno_barplot(as.numeric(data["N:SAMP:MUTATION_RATE",]), 
                                           bar_width = 0.75, 
                                           border = FALSE, 
                                           axis = TRUE,
                                           axis_param = list(gp = gpar(fontsize = 3, lwd = 0.5)),
                                           gp = gpar(col = NA, fill = "black")
),
Fragmentation = anno_barplot(frag, 
                   bar_width = 0.75,
                   border = FALSE, 
                   axis = TRUE,
                   axis_param = list(gp = gpar(fontsize = 3, lwd = 0.5)),
                   gp = gpar(col = NA, fill = "black")
),
Blasts = anno_barplot(as.numeric(data["N:CLIN:Percentage_BM_Blast",]), 
                   bar_width = 0.75, 
                   border = FALSE, 
                   axis = TRUE,
                   axis_param = list(gp = gpar(fontsize = 3, lwd = 0.5)),
                   gp = gpar(col = NA, fill = "black")
),
df = annot,
col = list(cluster = structure(brewer.pal(length(unique(annot$cluster)), "Set1"), 
                                          names = as.character(unique(annot$cluster))[c(4,6,2,1,5,3,7)]),
                      FAB = structure(c("grey70", rev(brewer.pal(length(unique(annot$FAB))-1, "Spectral"))), 
                                      names = as.character(unique(annot$FAB))[c(9,1,4,2,8,3,6,7,5)])
),
annotation_height = c(3, 3, 3, 3, 3),
annotation_legend_param = list(cluster = list(title = "t-SNE cluster", title_gp = gpar(fontsize = 5, fontface = "bold"), 
                                                labels_gp = gpar(fontsize = 5), grid_height = unit(0.2, "cm"), grid_width = unit(2, "mm")),
                               FAB = list(title = "FAB subtype", title_gp = gpar(fontsize = 5, fontface = "bold"), 
                                         labels_gp = gpar(fontsize = 5), grid_height = unit(0.2, "cm"), grid_width = unit(2, "mm"))
),
gap = unit(0.75, "mm"),
show_annotation_name = F,
height = unit(1.5, "cm")
)

cor_adj_p <- cor$adj.p[match(feat, cor$featureB)]

zero_col_mat = matrix(nrow = nrow(mat), ncol = 0)
rownames(zero_col_mat) <- cor_adj_p


# oncoprint 

ht <- Heatmap(mat,
        name = "TCGA AML",
        col = c("0" = "#EDEDED", "Amplification" = "#DC0000FF", "Deletion" = "#3C5488FF", "Nonsynonymous mutation" = "black", "NA" = "grey70"), 
        rect_gp = gpar(col= "white", lwd = unit(0.5, "mm")),
        top_annotation = ha1, 
        bottom_annotation = ha2,
        row_names_side = "left",
        row_names_gp = gpar(fontsize = 7),
        show_column_names = FALSE, 
        show_column_dend = FALSE, 
        cluster_columns = FALSE,
        show_heatmap_legend = TRUE,
        heatmap_legend_param = list(title = "Alteration",
                                    title_gp = gpar(fontsize = 5, fontface = "bold"),
                                    labels_gp = gpar(fontsize = 5),
                                    grid_height = unit(0.2, "cm"),
                                    grid_width = unit(2, "mm"))
) +
  Heatmap(zero_col_mat,
          row_names_gp = gpar(fontsize = 5))


# print pdf
pdf("FIgure3D_TCGA_AML_cytscore_oncoprint.pdf", height = 2, width = 6)
draw(ht, padding = unit(c(2, 15, 2, 2), "mm"))
decorate_annotation("Cytolytic_score", {
  grid.text("Cytolytic score", unit(0, "npc") - unit(6, "mm"), 0.5, 
            default.units = "npc", just = "right", gp = gpar(fontsize = 7))
})
decorate_annotation("Mutations", {
  grid.text("Mutations", unit(0, "npc") - unit(6, "mm"), 0.5, 
            default.units = "npc", just = "right", gp = gpar(fontsize = 5))
})
decorate_annotation("Fragmentation", {
  grid.text("Genome fragmentation %", unit(0, "npc") - unit(6, "mm"), 0.5, 
            default.units = "npc", just = "right", gp = gpar(fontsize = 5))
})
decorate_annotation("Blasts", {
  grid.text("BM blast %", unit(0, "npc") - unit(6, "mm"), 0.5, 
            default.units = "npc", just = "right", gp = gpar(fontsize = 5))
})
decorate_annotation("cluster", {
  grid.text("t-SNE cluster", unit(0, "npc") - unit(6, "mm"), 0.5, 
            default.units = "npc", just = "right", gp = gpar(fontsize = 5))
})
decorate_annotation("FAB", {
  grid.text("FAB", unit(0, "npc") - unit(6, "mm"), 0.5, 
            default.units = "npc", just = "right", gp = gpar(fontsize = 5))
})
dev.off()

