
# Plot oncoprint of DLBCL (GSE98588 Chapuy et al.) data CGA results (Figure 6H)

library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(grid)
library(dplyr)
library(readxl)


# load data
load("GSE98588_fm.Rdata")

# read correlation results
cor <- read.table("TableS6_GSE98588_DLBCL_antigen_correlations.tsv", sep = "\t", header = TRUE)
cor_gcb <- read.table("TableS6_GSE98588_DLBCL_antigen_correlations_ABC.tsv", sep = "\t", header = TRUE)
cor_abc <- read.table("TableS6_GSE98588_DLBCL_antigen_correlations_GCB.tsv", sep = "\t", header = TRUE)

# remove case with testicular involvement
fm <- fm[,!fm["B:CLIN:hodz_Testicular_invovlement",]%in%c(1)]

# order samples by CGA number
fm_cga_order <- fm[,order(as.numeric(fm["N:SAMP:nCGA",]), decreasing = FALSE)]

# create matrix with selected alterations
cor_noduplicates <- cor[grepl("nCGA", cor$featureA)&!grepl("DEL_LOSS|AMP_GAIN", cor$featureB),]

gnab_top <- cor_noduplicates[order(cor_noduplicates$p),]
gnab_top <- gnab_top[grepl("GNAB", gnab_top$featureB),]
gnab_top <- gnab_top[gnab_top$FDR<0.1,]

cnvr_top <- cor_noduplicates[order(cor_noduplicates$p),]
cnvr_top <- cnvr_top[grepl("CNVR", cnvr_top$featureB),]
cnvr_top <- cnvr_top[cnvr_top$FDR<0.05,]

# combine cnv and mutations and edit feat to match fm
feat_original <- rbind(cnvr_top, gnab_top) 

feat <- feat_original %>%
  mutate(featureB = gsub("AMP:GAIN", "AMP_GAIN", 
                         gsub("_GAIN", ":GAIN",
                              gsub("_LOSS", ":LOSS",
                                   gsub("_AMP", ":AMP",
                                        gsub("_DEL", "_DEL",
                                             gsub("SV_BCL2", "BCL2:SV", featureB))))))) %>%
  as.data.frame()

mat <- fm_cga_order[as.character(feat[,"featureB"]),]

# change codes for CNVs 
# mat[grepl(":AMP|:GAIN", rownames(mat)),][mat[grepl(":AMP|:GAIN", rownames(mat)),]==1] <- "Amplification/gain"
mat[grepl(":LOSS|:DEL", rownames(mat)),][mat[grepl(":LOSS|:DEL", rownames(mat)),]==1] <- "Deletion/loss"
mat[grepl("GNAB", rownames(mat)),][mat[grepl("GNAB", rownames(mat)),]==1] <- "Nonsynonymous mutation" 

# change row names
rownames(mat) <- gsub("Q$", "q",
                      gsub("([0-9])P$", "\\1p", 
                           gsub("([0-9])_([0-9])", "\\1.\\2",
                                gsub("_nonsynonymous|:AMP|:GAIN|:DEL|:LOSS|:LOSS_DEL|_GAIN|_DEL|_low_grade", "",
                                     gsub(":::::", "",
                                          gsub(".?:....:", "", 
                                               gsub("SV:BCL2", "BCL2 (SV)", 
                                                    gsub("([0-9])Q([0-9])", "\\1q\\2",
                                                         gsub("([0-9])P([0-9])", "\\1p\\2", 
                                                              rownames(mat))
                                                    )
                                               )
                                          )
                                     )))))

mat <- as.matrix(mat)

# create data frames for heatmap annotations
coo <- as.character(fm_cga_order["B:SAMP:COO_byGEP_ABC",])
coo[coo=="1"] <- "ABC"
coo[fm_cga_order["B:SAMP:COO_byGEP_GCB",]==1] <- "GCB"
coo[fm_cga_order["B:SAMP:COO_byGEP_Unclassified",]==1] <- "Unclassified"

# data frame of sample annotations
annot <- droplevels(data.frame(COO_byGEP = coo,
                               IPI = as.factor(fm_cga_order["N:CLIN:IPI",])))
levels(annot$IPI) <- c("NA", "0", "1", "2", "3", "4", "5")
annot[is.na(annot)] <- "NA"

# change mutation load of hypermutated sample
fm_cga_order["N:SAMP:numberOfMutations",][fm_cga_order["N:SAMP:numberOfMutations",]==5956] <- 505

## create heatmap annotations

ha1 <- HeatmapAnnotation(CGA = anno_barplot(as.numeric(fm_cga_order["N:SAMP:nCGA",]), 
                                                  bar_width = 0.75, 
                                                  border = FALSE, 
                                                  axis = TRUE,
                                                  axis_param = list(gp = gpar(fontsize = 5, lwd = 0.5)),
                                                  gp = gpar(col = NA, fill = "grey50"),
                                                  ylim = c(0,8)),
                         show_annotation_name = F,
                         height = unit(0.5, "cm"))


ha2 = HeatmapAnnotation(Mutations = anno_barplot(as.numeric(fm_cga_order["N:SAMP:numberOfMutations",]), 
                                           bar_width = 0.75, 
                                           border = FALSE, 
                                           axis = TRUE,
                                           axis_param = list(gp = gpar(fontsize = 5, lwd = 0.5)),
                                           gp = gpar(col = NA, fill = "black")
),
CNVs = anno_barplot(as.numeric(fm_cga_order["N:SAMP:numberOfCNAs",]), 
                  bar_width = 0.75, 
                   border = FALSE, 
                   axis = TRUE,
                  axis_param = list(gp = gpar(fontsize = 5, lwd = 0.5)),
                  gp = gpar(col = NA, fill = "black")
),
Purity = anno_barplot(as.numeric(fm_cga_order["N:SAMP:purity_absolute_reviewed",]), 
                   bar_width = 0.75, 
                   border = FALSE, 
                   axis = TRUE,
                   axis_param = list(gp = gpar(fontsize = 5, lwd = 0.5)),
                   gp = gpar(col = NA, fill = "black")
),
df = annot,
col = list(COO_byGEP = structure(brewer.pal(length(unique(annot$COO_byGEP)), "Set1"), 
                                 names = as.character(unique(annot$COO_byGEP))),
           IPI = structure(gsub("#D53E4F", "grey", rev(brewer.pal(length(unique(annot$IPI)), "Spectral"))), 
                           names = sort(as.character(unique(annot$IPI))))
),
annotation_height = c(3, 3, 3, 3, 3),
annotation_legend_param = list(COO_byGEP = list(title = "Subtype", title_gp = gpar(fontsize = 5), 
                                                labels_gp = gpar(fontsize = 5), grid_height = unit(0.2, "cm"), grid_width = unit(2, "mm")),
                               IPI = list(title = "IPI", title_gp = gpar(fontsize = 5), 
                                          labels_gp = gpar(fontsize = 5), grid_height = unit(0.2, "cm"), grid_width = unit(2, "mm"))
),
gap = unit(0.75, "mm"),
show_annotation_name = F,
height = unit(1.5, "cm")
)

cor_adj_p <- cor_noduplicates$FDR[match(feat_original$featureB, cor_noduplicates$featureB)]

zero_col_mat = matrix(nrow = nrow(mat), ncol = 0)
rownames(zero_col_mat) <- cor_adj_p

# make heatmap
ht <- Heatmap(mat,
              name = "oncoprint",
              col = c("0" = "#EDEDED", "Nonsynonymous mutation" = "black", "Amplification/gain" = "#DC0000FF", "Deletion/loss" = "#3C5488FF", "Structural variation" = "#00A087FF"), 
              rect_gp = gpar(col= "white", lwd = unit(0.4, "mm")),
              top_annotation = ha1, 
              bottom_annotation = ha2,
              row_names_side = "left",
              row_names_gp = gpar(fontsize = 7),
              show_column_names = FALSE, 
              show_column_dend = FALSE,
              cluster_columns = FALSE,
              cluster_rows = FALSE,
              split = c(rep("CNV/SV",1), rep("Mutation",3)),
              row_title_gp = gpar(fontsize = 5),
              show_heatmap_legend = TRUE,
              heatmap_legend_param = list(title = "Alteration",
                                          title_gp = gpar(fontsize = 5),
                                          labels_gp = gpar(fontsize = 5),
                                          grid_height = unit(0.2, "cm"),
                                          grid_width = unit(2, "mm"))
) +
  Heatmap(zero_col_mat,
          row_names_gp = gpar(fontsize = 5))


# oncoprint
pdf("Figure6H_DLBCL_GSE98588_CGA_oncoprint.pdf", height = 1.5, width = 4.5)
draw(ht, padding = unit(c(2, 5, 2, 2), "mm"))
decorate_annotation("CGA", {
  grid.text("CGA number", unit(0, "npc") - unit(6, "mm"), 0.5, 
            default.units = "npc", just = "right", gp = gpar(fontsize = 7))
})
decorate_annotation("Mutations", {
  grid.text("Mutations", unit(0, "npc") - unit(6, "mm"), 0.5, 
            default.units = "npc", just = "right", gp = gpar(fontsize = 5))
})
decorate_annotation("CNVs", {
  grid.text("CNVs", unit(0, "npc") - unit(6, "mm"), 0.5, 
            default.units = "npc", just = "right", gp = gpar(fontsize = 5))
})
decorate_annotation("Purity", {
  grid.text("Purity", unit(0, "npc") - unit(6, "mm"), 0.5, 
            default.units = "npc", just = "right", gp = gpar(fontsize = 5))
})
decorate_annotation("COO_byGEP", {
  grid.text("Subtype", unit(0, "npc") - unit(6, "mm"), 0.5, 
            default.units = "npc", just = "right", gp = gpar(fontsize = 5))
})
decorate_annotation("IPI", {
  grid.text("IPI", unit(0, "npc") - unit(6, "mm"), 0.5, 
            default.units = "npc", just = "right", gp = gpar(fontsize = 5))
})
dev.off()


