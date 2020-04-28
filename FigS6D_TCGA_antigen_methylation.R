
# Plot Figure S6D heatmaps of antigen methylation for individual probes ordered by expression (TCGA DLBCL and AML)

library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(dplyr)
library(readxl)
library(edgeR)
library(limma)

# load data
# gene expression
aml_gexp <- get(load("LAML.rnaseqv2.counts.Rdata"))
dlbcl_gexp <- get(load("DLBC.rnaseqv2.counts.Rdata"))

# normalize and transform using voom
aml_gexp <- voom(aml_gexp+0.01, normalize.method = "quantile")$E
dlbcl_gexp <- voom(dlbcl_gexp+0.01, normalize.method = "quantile")$E

# methylation
aml_meth <- get(load("TCGA_AML_meth_probes_genelist.Rdata"))
dlbcl_meth <- get(load("TCGA_DLBCL_meth_probes_genelist.Rdata"))

# select matching cases
aml_gexp <- aml_gexp[,colnames(aml_gexp) %in% colnames(aml_meth)]
aml_meth <- aml_meth[,colnames(aml_meth) %in% colnames(aml_gexp)]
aml_meth <- aml_meth[,colnames(aml_gexp)]

dlbcl_gexp <- dlbcl_gexp[,colnames(dlbcl_gexp) %in% colnames(dlbcl_meth)]
dlbcl_meth <- dlbcl_meth[,colnames(dlbcl_meth) %in% colnames(dlbcl_gexp)]
dlbcl_meth <- dlbcl_meth[,colnames(dlbcl_gexp)]


# function for heatmaps
plot_heatmap <- function(gene){


# select gene 
meth_annot <- meth_annot_extra
meth_annot_gene <- meth_annot_extra[meth_annot_extra$nearestTSS %in% c(gene) & meth_annot_extra$distanceToTSS < 100,]

aml_meth_gene <- merge(aml_meth, meth_annot_gene, by.x = "row.names", by.y = "methProbeIDs")
dlbcl_meth_gene <- merge(dlbcl_meth, meth_annot_gene, by.x = "row.names", by.y = "methProbeIDs")


## -------------------------------------------------------------------------------

## plot heatmaps

# DLBCL

mat_dlbcl_meth <- data.matrix(dlbcl_meth_gene[2:49])

rownames(mat_dlbcl_meth) <- dlbcl_meth_gene$Row.names 
mat_dlbcl_gexp <- data.matrix(dlbcl_gexp)

mat_dlbcl_gexp <- mat_dlbcl_gexp[,order(mat_dlbcl_gexp[gene,])]

mat_dlbcl_meth <- mat_dlbcl_meth[,colnames(mat_dlbcl_gexp)]

# make 0 the smallest expression value
mat_dlbcl_gexp[mat_dlbcl_gexp<0] = 0


ha1 = HeatmapAnnotation(exprs = anno_barplot(mat_dlbcl_gexp[gene,],
                                             ylim = c(0, 6),
                                             bar_width = 0.75, 
                                             border = FALSE, 
                                             axis = TRUE,
                                             axis_param = list(gp = gpar(fontsize = 5, lwd = 0.5)),
                                             gp = gpar(col = NA, fill = "black")),
                        height = unit(0.75, "cm"),
                        show_annotation_name = F)

ha2 = rowAnnotation(DistanceToTSS = anno_barplot(meth_annot_gene$distanceToTSS, 
                                       which = "row",
                                       bar_width = 0.75, 
                                       border = FALSE, 
                                       axis = TRUE,
                                       axis_param = list(gp = gpar(fontsize = 5, lwd = 0.5)),
                                       gp = gpar(col = NA, fill = "grey30")))

hm_dlbcl <- Heatmap(mat_dlbcl_meth,
                    column_title = "TCGA DLBCL",
                    column_title_gp = gpar(fontsize = 9),
                    col = colorRamp2(seq(0.5, 1, 0.25), rev(c("#fea719", "#d4d3d1", "#a52a2a"))),
                    show_column_names = FALSE,
                    show_row_names = FALSE,
                    cluster_rows = TRUE,
                    cluster_columns = FALSE,
                    clustering_distance_rows = "spearman",
                    clustering_method_rows = "ward.D",
                    top_annotation = ha1, 
                    show_heatmap_legend = FALSE
)




## ----------------------------------------------------------------------------------

# AML

mat_aml_meth <- data.matrix(aml_meth_gene[2:171])

rownames(mat_aml_meth) <- aml_meth_gene$Row.names 
mat_aml_gexp <- data.matrix(aml_gexp)

mat_aml_gexp <- mat_aml_gexp[,order(mat_aml_gexp[gene,])]

mat_aml_meth <- mat_aml_meth[,colnames(mat_aml_gexp)]

# make 0 the smallest expression value
mat_aml_gexp[mat_aml_gexp<0] = 0


ha1 = HeatmapAnnotation(exprs = anno_barplot(mat_aml_gexp[gene,], 
                                             ylim = c(0, 6),
                                             bar_width = 0.75, 
                                             border = FALSE, 
                                             axis = TRUE,
                                             axis_param = list(gp = gpar(fontsize = 5, lwd = 0.5)),
                                             gp = gpar(col = NA, fill = "black")),
                        height = unit(0.75, "cm"),
                        show_annotation_name = F)

ha2 = rowAnnotation(DistanceToTSS = anno_barplot(meth_annot_gene$distanceToTSS, 
                                       which = "row",
                                       bar_width = 0.75, 
                                       border = FALSE, 
                                       axis = TRUE,
                                       axis_param = list(gp = gpar(fontsize = 5, lwd = 0.5)),
                                       gp = gpar(col = NA, fill = "grey30")),
                    width = unit(1, "cm"))

hm_aml <- Heatmap(mat_aml_meth,
                  column_title = "TCGA AML",
                  column_title_gp = gpar(fontsize = 9),
                  name = "Methylation beta value", 
                  col = colorRamp2(seq(0.5, 1, 0.25), rev(c("#fea719", "#d4d3d1", "#a52a2a"))),
                  row_names_side = "right", 
                  row_names_gp = gpar(fontsize = 7),
                  show_column_names = FALSE,
                  cluster_rows = FALSE,
                  cluster_columns = FALSE,
                  clustering_distance_rows = "spearman",
                  clustering_method_rows = "ward.D",
                  top_annotation = ha1, 
                  width = unit(2.25, "cm") # MAGEB2
) + ha2

## ----------------------------------------------------------------------------------

pdf(paste0("FigureS6D_TCGA_DLBCL_AML_methylation_", gene, "_heatmap.pdf"), height = 2, width = 6)
draw(hm_dlbcl + hm_aml,
     gap = unit(1, "cm"),
     padding = unit(c(1, 1, 0.2, 0.2), "cm"))
decorate_annotation("exprs", {
  grid.text(paste(gene), unit(0, "npc") - unit(6, "mm"), 0.5, 
            default.units = "npc", just = "right", gp = gpar(fontsize = 7, fontface = "italic"))
})
dev.off()
}

# plot heatmaps for MAGEB1 and MAGEB2
lapply(c("MAGEB1", "MAGEB2"), plot_heatmap)
