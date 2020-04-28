
# Plot oncoprint of CoMMpass data CGA results (Figure 6G)

library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(grid)
library(dplyr)
library(viridis)


# load data
load("MM_COMPASS_FM.Rdata")
load("GSVA_MM_COMPASS_scores.Rdata")

# read correlation results
cor <- read.table("TableS6_CoMMpass_MM_nCGA_correlations.tsv", sep = "\t", header = TRUE)

# only samples with gexp data
fm <- fm[, !is.na(fm["N:GEXP:KRAS",])]

# add subtypes
subtypes <- rep("CCND1", length(fm))
subtypes[fm["B:SAMP:cancermap_subtypes_WHSC1_FGFR3_Ig",]==1] <- "FGFR3"
subtypes[fm["B:SAMP:cancermap_subtypes_Hyperdiploid",]==1] <- "Hyperdiploid(gain(11q)"
subtypes[fm["B:SAMP:cancermap_subtypes_Hyperdiploid_amp1q",]==1] <- "Hyperdiploid/gain(1q)"
subtypes[fm["B:SAMP:cancermap_subtypes_MAF_Ig",]==1] <- "MAF"
subtypes[fm["B:SAMP:cancermap_subtypes_TRAF3_Aberrated",]==1] <- "TRAF3"
subtypes[fm["B:SAMP:cancermap_cluster_CGA_Prolif",]==1] <- "CGA/Proliferative"

# order fm and viz by CGA number
cga_order <- order(fm["N:SAMP:nCGA",])
fm <- fm[,cga_order]
viz_scores <- viz_scores[,cga_order]

# create matrix with selected GSVA scores
feats <- c("MYC_TARGETS_V1-MSIGDB_HALLMARKS",
           "DNA_REPLICATION-KEGG_MSIGDB_C2",
           "E2F_TARGETS-MSIGDB_HALLMARKS",
           "G2M_CHECKPOINT-MSIGDB_HALLMARKS",
           "MITOTIC_M_M_G1_PHASES-REACTOME_MSIGDB_C2",
           "INFLAMMATORY_RESPONSE-MSIGDB_HALLMARKS",
           "TNFA_SIGNALING_VIA_NFKB-MSIGDB_HALLMARKS",
           "BCR_SIGNALING_PATHWAY-SIG_MSIGDB_C2")
mat <- viz_scores[feats,]
mat <- mat[!is.na(mat[,1]),] # remove CGA genes not found

mat_scaled <- t(apply(mat, 1, scale))
colnames(mat_scaled) <- colnames(mat)
rownames(mat_scaled) <- gsub("_", " ",gsub("\\-.*", "", rownames(mat_scaled)))

# create matrix with selected genetic alteration features
alt <- data.frame(NRAS = as.character(fm["B:GNAB:NRAS",]),
                          `CCND1-Ig` = as.character(fm["B:CNVR:SeqWGS_CCND1_Ig_translocation",]),
                          `MAF-Ig` = as.character(fm["B:CNVR:SeqWGS_MAF_Ig_translocation",]))
colnames(alt) <- c("NRAS", "CCND1-Ig", "MAF-Ig")

# create heatmap annotations
annot <- data.frame(subtype = subtypes[cga_order])

gexp <- data.frame(`HLA II score` = as.numeric(fm["N:SAMP:HLAIIScore",]))
colnames(gexp) <- "HLA II score"

gexp_scaled <- as.data.frame(apply(gexp, 2, scale))
mat_hla2 <- t(gexp_scaled)
colnames(mat_hla2) <- colnames(mat)

ramp <- colorRamp2(seq(-2, 2, length.out = 11), rev(brewer.pal(11, "RdBu")))
green <- structure(c("darkgreen", "#EDEDED", "grey50"), names = c("1", "0", "NA"))
blackgrey <- structure(c("black", "#EDEDED", "grey50"), names = c("1", "0", "NA"))

mut <- as.numeric(fm["N:CLIN:NS_Non_IG_Variants",])
mut[mut > 150] <- 150 # mutation load cutoff due to hypermutated samples

ha1 <- HeatmapAnnotation(`# CGA` = anno_barplot(as.numeric(fm["N:SAMP:nCGA",]), 
                                            bar_width = 1, 
                                            border = FALSE, 
                                            axis = TRUE,
                                            axis_param = list(gp = gpar(fontsize = 5, lwd = 0.5)),
                                            gp = gpar(col = NA, fill = "grey30"),
                                            ylim = c(0, 13)),
                         annotation_name_gp = gpar(fontsize = 10),
                         height = unit(1, "cm"))


ha2 <- HeatmapAnnotation(df = alt, col = list(NRAS = blackgrey,
                                              `CCND1-Ig` = green,
                                              `MAF-Ig` = green),
                         annotation_name_gp = gpar(fontsize = 10),
                         show_legend = F,
                         height = unit(1, "cm"))

ha3 <- HeatmapAnnotation(`# Hyperdiploid chr` = anno_barplot(as.numeric(fm["N:CNVR:Hyperdiploid_Chr_Count",]), 
                                                      bar_width = 0.75, 
                                                      border = FALSE, 
                                                      axis = TRUE,
                                                      axis_param = list(gp = gpar(fontsize = 5, lwd = 0.5)),
                                                      gp = gpar(col = NA, fill = "grey20")),
                         
                         Proliferation = anno_barplot(as.numeric(fm["N:CLIN:Prolif_Index",]), 
                                                      bar_width = 0.75, 
                                                      border = FALSE, 
                                                      axis = TRUE,
                                                      axis_param = list(gp = gpar(fontsize = 5, lwd = 0.5)),
                                                      gp = gpar(col = NA, fill = "grey20")),
                         
                         Mutations = anno_barplot(mut, 
                                                   bar_width = 0.75, 
                                                   border = FALSE, 
                                                   axis = TRUE,
                                                   axis_param = list(gp = gpar(fontsize = 5, lwd = 0.5)),
                                                   gp = gpar(col = NA, fill = "grey20"),
                                                   ylim = c(0, 150)),
                         annotation_name_gp = gpar(fontsize = 10),
                         show_legend = F,
                         height = unit(1, "cm")
)
                                


# make heatmaps
ht <- Heatmap(mat_scaled,
              name = "oncoprint",
              col = colorRamp2(seq(-1, 1, length.out = 11), viridis(11, begin = 0.1, end = 0.9)),
              row_names_side = "right",
              row_names_gp = gpar(fontsize = 6),
              show_column_names = FALSE, 
              show_column_dend = FALSE,
              cluster_columns = FALSE,
              cluster_rows = FALSE,
              show_row_dend = FALSE,
              row_title_gp = gpar(fontsize = 8),
              show_heatmap_legend = T,
              heatmap_legend_param = list(title = "GSVA\nZ-score",
                                          title_gp = gpar(fontsize = 10),
                                          labels_gp = gpar(fontsize = 10),
                                          grid_height = unit(0.2, "cm"),
                                          grid_width = unit(2, "mm"),
                                          legend_direction = "horizontal",
                                          title_position = "topcenter")
) 

ht_hla <- Heatmap(mat_hla2,
                  name = "HLA II score",
                  col = ramp,
                  row_names_side = "right",
                  row_names_gp = gpar(fontsize = 8),
                  show_column_names = FALSE, 
                  show_column_dend = FALSE,
                  cluster_columns = FALSE,
                  cluster_rows = FALSE,
                  show_row_dend = FALSE,
                  show_heatmap_legend = T,
                  heatmap_legend_param = list(title = "HLA II\nZ-score",
                                      title_gp = gpar(fontsize = 10),
                                      labels_gp = gpar(fontsize = 10),
                                      grid_height = unit(0.2, "cm"),
                                      grid_width = unit(2, "mm"),
                                      legend_direction = "horizontal",
                                      title_position = "topcenter")
  )

# combine heatmaps and annotations
ht_list <- ha1 %v% ha2 %v% ha3 %v% ht_hla %v% ht

# print oncoprint
pdf("Figure6G_CoMMpass_CGA_oncoprint.pdf", height = 3.5, width = 5)
draw(ht_list, padding = unit(c(2, 5, 2, 2), "mm"), heatmap_legend_side = "bottom")
dev.off()

