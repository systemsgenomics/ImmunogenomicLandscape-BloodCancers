
# Plot heatmap of CGA expression in CoMMpass data (Figure S6E)

library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(grid)
library(dplyr)

# load data
load("MM_COMPASS_FM.Rdata")
load("GSVA_MM_COMPASS_scores.Rdata")

# read correlation results
cor <- read.table("TableS6_CoMMpass_MM_nCGA_correlations.tsv", sep = "\t", header = TRUE)

# only samples with gexp data
fm <- fm[, !is.na(fm["N:GEXP:KRAS",])]

# order fm by subtype
subtypes <- rep("CCND1", length(fm))
subtypes[fm["B:SAMP:cancermap_subtypes_WHSC1_FGFR3_Ig",]==1] <- "FGFR3"
subtypes[fm["B:SAMP:cancermap_subtypes_Hyperdiploid",]==1] <- "Hyperdiploid(gain(11q)"
subtypes[fm["B:SAMP:cancermap_subtypes_Hyperdiploid_amp1q",]==1] <- "Hyperdiploid/gain(1q)"
subtypes[fm["B:SAMP:cancermap_subtypes_MAF_Ig",]==1] <- "MAF"
subtypes[fm["B:SAMP:cancermap_subtypes_TRAF3_Aberrated",]==1] <- "TRAF3"
subtypes[fm["B:SAMP:cancermap_cluster_CGA_Prolif",]==1] <- "CGA/Proliferative"

subtypes_order <- order(subtypes, fm["N:SAMP:nCGA",])
fm <- fm[,subtypes_order]

# create matrix with selected alterations
cga <- as.character(read.table("cga.txt")[,1])

feats <- paste0("N:GEXP:", cga)
mat <- fm[feats,]
mat <- mat[!is.na(mat[,1]),] # remove CGA genes not found

mat_scaled <- t(apply(mat, 1, scale))
colnames(mat_scaled) <- colnames(mat)
rownames(mat_scaled) <- gsub("N:GEXP:", "", rownames(mat_scaled))


# create heatmap annotations
annot <- data.frame(subtype = subtypes[subtypes_order])

ha1 <- HeatmapAnnotation(df = annot, col = list(subtype = structure(brewer.pal(length(unique(annot$subtype)), "Set1"), 
                                                                                             names = as.character(unique(annot$subtype)))),
                                                  
                          CGA = anno_barplot(as.numeric(fm["N:SAMP:nCGA",]), 
                                                  bar_width = 0.75, 
                                                  border = FALSE, 
                                                  axis = TRUE,
                                                  axis_param = list(gp = gpar(fontsize = 5, lwd = 0.5)),
                                                  gp = gpar(col = NA, fill = "#a52a2a"),
                                                  ylim = c(0, 13)),
                         
                         Proliferation = anno_barplot(as.numeric(fm["N:CLIN:Prolif_Index",]), 
                                               bar_width = 0.75, 
                                               border = FALSE, 
                                               axis = TRUE,
                                               axis_param = list(gp = gpar(fontsize = 5, lwd = 0.5)),
                                               gp = gpar(col = NA, fill = "grey20")),
                         
                         Mutations = anno_barplot(as.numeric(fm["N:CLIN:NS_Non_IG_Variants",]), 
                                              bar_width = 0.75, 
                                              border = FALSE, 
                                              axis = TRUE,
                                              axis_param = list(gp = gpar(fontsize = 5, lwd = 0.5)),
                                              gp = gpar(col = NA, fill = "grey20"),
                                              ylim = c(0, 900)),
                         
                                           gap = unit(0.75, "mm"),
annotation_legend_param = list(subtype = list(title = "Subtype", title_gp = gpar(fontsize = 5), 
                                                labels_gp = gpar(fontsize = 5), grid_height = unit(0.2, "cm"), grid_width = unit(2, "mm"))),
height = unit(1.5, "cm"),
simple_anno_size_adjust = T,
show_annotation_name = F
)

alt <- data.frame(ccnd1 = as.character(fm["B:CNVR:RNASeq_CCND1_Ig_translocation",]),
                  fgfr3 = as.character(fm["B:CNVR:RNASeq_WHSC1_Ig_translocation",]),
                  hyperdiploid = as.character(fm["B:CNVR:Hyperdiploid_Call",]),
                  maf = as.character(fm["B:CNVR:RNASeq_MAF_Ig_translocation",]),
                  traf3 = as.character(fm["B:GNAB:TRAF3",]),
                  nras = as.character(fm[c("B:GNAB:NRAS"),]),
                  traf3_cnv = as.numeric(fm["N:CNVR:TRAF3",]),
                  gain11q = as.numeric(fm["N:CNVR:11q23",]),
                  gain1q =  as.numeric(fm["N:CNVR:1q21",])
                  )

green <- structure(c("darkgreen", "#EDEDED", "grey50"), names = c("1", "0", "NA"))
blackgrey <- structure(c("black", "#EDEDED", "grey50"), names = c("1", "0", "NA"))
ramp <- colorRamp2(seq(-1, 1, length.out = 11), rev(brewer.pal(11, "PuOr")))

ha2 = HeatmapAnnotation(df = alt, col = list(ccnd1 = green,
                                             fgfr3 = green,
                                             hyperdiploid = green,
                                             maf = green,
                                             traf3 = blackgrey,
                                             nras = blackgrey,
                                             traf3_cnv = ramp,
                                             gain11q = ramp,
                                             gain1q = ramp),
                        height = unit(1.75, "cm"),
                        simple_anno_size_adjust = T,
                        show_legend = F,
                        show_annotation_name = F,
                        gap = unit(0, "mm"))

# make heatmap
ht <- Heatmap(mat_scaled,
              name = "heatmap",
              col = colorRamp2(seq(-3, 3, length.out = 11), rev(brewer.pal(11, "RdBu"))),
              use_raster = T,
              top_annotation = ha1, 
              bottom_annotation = ha2,
              row_names_side = "right",
              row_names_gp = gpar(fontsize = 5),
              show_column_names = FALSE, 
              show_column_dend = FALSE,
              cluster_columns = FALSE,
              cluster_rows = TRUE,
              show_row_dend = FALSE,
              row_title_gp = gpar(fontsize = 5),
              show_heatmap_legend = T,
              heatmap_legend_param = list(title = "Expression (log2)\nZ-score",
                                          title_gp = gpar(fontsize = 5),
                                          labels_gp = gpar(fontsize = 5),
                                          grid_height = unit(0.2, "cm"),
                                          grid_width = unit(2, "mm"))
)


# print heatmap
pdf("FigureS6E_CoMMpass_CGA_heatmap.pdf", height = 3.25, width = 5)
draw(ht)
decorate_annotation("CGA", {
  grid.text("CGAs expressed", unit(0, "npc") + unit(0.15, "cm"), 0.5, 
            default.units = "npc", just = c("left", "bottom"), gp = gpar(fontsize = 5))
})
decorate_annotation("Proliferation", {
  grid.text("Prolif. index", unit(0, "npc") + unit(0.15, "cm"), 0.5, 
            default.units = "npc", just = c("left", "bottom"), gp = gpar(fontsize = 5))
})
decorate_annotation("Mutations", {
  grid.text("Mutations", unit(0, "npc") + unit(0.15, "cm"), 0.5, 
            default.units = "npc", just = c("left", "bottom"), gp = gpar(fontsize = 5))
})
decorate_annotation("ccnd1", {
  grid.text("CCND1", unit(0, "npc") + unit(7.1, "cm"), 0.1, 
            default.units = "npc", just = c("left", "bottom"), gp = gpar(fontsize = 5))
})
decorate_annotation("fgfr3", {
  grid.text("FGFR3/MMSET", unit(0, "npc") + unit(7.1, "cm"), 0.1, 
            default.units = "npc", just = c("left", "bottom"), gp = gpar(fontsize = 5))
})
decorate_annotation("hyperdiploid", {
  grid.text("Hyperdiploid", unit(0, "npc") + unit(7.1, "cm"), 0.1, 
            default.units = "npc", just = c("left", "bottom"), gp = gpar(fontsize = 5))
})
decorate_annotation("maf", {
  grid.text("MAF", unit(0, "npc") + unit(7.1, "cm"), 0.1, 
            default.units = "npc", just = c("left", "bottom"), gp = gpar(fontsize = 5))
})
decorate_annotation("traf3", {
  grid.text("TRAF3", unit(0, "npc") + unit(7.1, "cm"), 0.1, 
            default.units = "npc", just = c("left", "bottom"), gp = gpar(fontsize = 5))
})
decorate_annotation("nras", {
  grid.text("NRAS", unit(0, "npc") + unit(7.1, "cm"), 0.1, 
            default.units = "npc", just = c("left", "bottom"), gp = gpar(fontsize = 5))
})
decorate_annotation("traf3_cnv", {
  grid.text("TRAF3", unit(0, "npc") + unit(7.1, "cm"), 0.1, 
            default.units = "npc", just = c("left", "bottom"), gp = gpar(fontsize = 5))
})
decorate_annotation("gain11q", {
  grid.text("11q23", unit(0, "npc") + unit(7.1, "cm"), 0.1, 
            default.units = "npc", just = c("left", "bottom"), gp = gpar(fontsize = 5))
})
decorate_annotation("gain1q", {
  grid.text("1q21", unit(0, "npc") + unit(7.1, "cm"), 0.1, 
            default.units = "npc", just = c("left", "bottom"), gp = gpar(fontsize = 5))
})

dev.off()
