
# Plot Figure 3A and Figure S3B-C oncoprints of genetic alterations correlated with cytolytic score in DLBCL (GSE98588 Chapuy et al.)

library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(grid)
library(dplyr)


# load data
load("GSE98588_fm.Rdata")

# read correlation results
cor <- read.table("TableS3_GSE98588_DLBCL_CytScore_correlations.tsv", sep = "\t", header = TRUE, stringsAsFactors = F)
cor_gcb <- read.table("TableS3_GSE98588_DLBCL_CytScore_correlations_GCB_all.tsv", sep = "\t", header = TRUE, stringsAsFactors = F)
cor_abc <- read.table("TableS3_GSE98588_DLBCL_CytScore_correlations_ABC_all.tsv", sep = "\t", header = TRUE, stringsAsFactors = F)

# order samples by Cytolytic score
fm_cyt_order <- fm[,order(as.numeric(fm["N:SAMP:CytolyticScore",]), decreasing = TRUE)]

# create matrix with selected alterations
cor_noduplicates <- cor[!grepl("6P_GAIN|7Q_GAIN|7P_AMP$|11Q_GAIN", cor$featureB),]

gnab_top4 <- cor_noduplicates[order(cor_noduplicates$p),]
gnab_top4 <- gnab_top4[grepl("GNAB", gnab_top4$featureB),]
gnab_top4 <- gnab_top4[c(1:6),]

cnvr_top4 <- cor_noduplicates[order(cor_noduplicates$p),]
cnvr_top4 <- cnvr_top4[grepl("CNVR", cnvr_top4$featureB),]
cnvr_top4 <- cnvr_top4[c(1:6),]

# combine cnv and mutations and edit feat to match fm
feat_original <- rbind(cnvr_top4, gnab_top4) 

mat <- fm_cyt_order[as.character(c(cnvr_top4$featureB, gnab_top4$featureB)),]

# change codes for CNVs 
mat[grepl(":AMP|:GAIN", rownames(mat)),][mat[grepl(":AMP|:GAIN", rownames(mat)),]==1] <- "Amplification"
mat[grepl(":LOSS|:DEL", rownames(mat)),][mat[grepl(":LOSS|:DEL", rownames(mat)),]==1] <- "Deletion"
mat[grepl(":SV", rownames(mat)),][mat[grepl(":SV", rownames(mat)),]==1] <- "Structural variation"
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
coo <- as.character(fm_cyt_order["B:SAMP:COO_byGEP_ABC",])
coo[coo=="1"] <- "ABC"
coo[fm_cyt_order["B:SAMP:COO_byGEP_GCB",]==1] <- "GCB"
coo[fm_cyt_order["B:SAMP:COO_byGEP_Unclassified",]==1] <- "Unclassified"

# data frame of sample annotations
annot <- droplevels(data.frame(COO_byGEP = coo,
                               IPI = as.factor(fm_cyt_order["N:CLIN:IPI",])))
levels(annot$IPI) <- c("NA", "0", "1", "2", "3", "4", "5")
annot[is.na(annot)] <- "NA"


# change mutation load of hypermutated sample
fm_cyt_order["N:SAMP:numberOfMutations",][fm_cyt_order["N:SAMP:numberOfMutations",]==5956] <- 505

## create heatmap annotations

ha1 <- HeatmapAnnotation(Cytolytic_score = anno_barplot(as.numeric(fm_cyt_order["N:SAMP:CytolyticScore",]), 
                                                  bar_width = 0.75, 
                                                  border = FALSE, 
                                                  axis = TRUE,
                                                  axis_param = list(gp = gpar(fontsize = 5, lwd = 0.5)),
                                                  gp = gpar(col = NA, fill = "grey50"),
                                                  ylim = c(3.75,11)),
                         height = unit(0.5, "cm"),
                         show_annotation_name = F
)


ha2 = HeatmapAnnotation(Mutations = anno_barplot(as.numeric(fm_cyt_order["N:SAMP:numberOfMutations",]), 
                                           bar_width = 0.75, 
                                           border = FALSE, 
                                           axis = TRUE,
                                           axis_param = list(gp = gpar(fontsize = 3, lwd = 0.5)),
                                           gp = gpar(col = NA, fill = "black")
),
CNVs = anno_barplot(as.numeric(fm_cyt_order["N:SAMP:numberOfCNAs",]), 
                   bar_width = 0.75, 
                   border = FALSE, 
                   axis = TRUE,
                   axis_param = list(gp = gpar(fontsize = 3, lwd = 0.5)),
                   gp = gpar(col = NA, fill = "black")
),
Purity = anno_barplot(as.numeric(fm_cyt_order["N:SAMP:purity_absolute_reviewed",]), 
                   bar_width = 0.75, 
                   border = FALSE, 
                   axis = TRUE,
                   axis_param = list(gp = gpar(fontsize = 3, lwd = 0.5)),
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

mat=gsub(" 0", "0", mat)

# make heatmap
ht <- Heatmap(mat,
              name = "oncoprint",
              col = c("0" = "#EDEDED", "Nonsynonymous mutation" = "black", "Amplification" = "#DC0000FF", "Deletion" = "#3C5488FF", "Structural variation" = "#00A087FF"), 
              rect_gp = gpar(col= "white", lwd = unit(0.5, "mm")),
              top_annotation = ha1, 
              #top_annotation_height = unit(0.75, "cm"), 
              bottom_annotation = ha2,
              #bottom_annotation_height = unit(1.5, "cm"),
              row_names_side = "left",
              row_names_gp = gpar(fontsize = 7),
              show_column_names = FALSE, 
              show_column_dend = FALSE,
              cluster_columns = FALSE,
              cluster_rows = FALSE,
              split = c(rep("CNV/SV",6), rep("Mutation",6)),
              row_title_gp = gpar(fontsize = 5),
              show_heatmap_legend = TRUE,
              heatmap_legend_param = list(title = "Alteration",
                                          title_gp = gpar(fontsize = 5),
                                          labels_gp = gpar(fontsize = 5),
                                          grid_height = unit(0.2, "cm"),
                                          grid_width = unit(2, "mm"))
) +
  #ha3 +
  Heatmap(zero_col_mat,
          row_names_gp = gpar(fontsize = 5))


# oncoprint
pdf("Figure3A_DLBCL_cytscore_oncoprint.pdf", height = 2.25, width = 5)
draw(ht, padding = unit(c(2, 5, 2, 2), "mm"))
decorate_annotation("Cytolytic_score", {
  grid.text("Cytolytic score", unit(0, "npc") - unit(6, "mm"), 0.5, 
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



## ----------------------



## ABC

cor_abc_noduplicates <- cor_abc[!grepl("6P_GAIN|7Q_GAIN|7P_AMP$|11Q_GAIN", cor_abc$featureB),] # the top AMP_GAIN features are the same as GAIN only, remove these

gnab_top4_abc <- cor_abc_noduplicates[order(cor_abc_noduplicates$p),]
gnab_top4_abc <- gnab_top4_abc[grepl("GNAB", gnab_top4_abc$featureB),]
gnab_top4_abc <- gnab_top4_abc[c(1:6),]

cnvr_top4_abc <- cor_abc_noduplicates[order(cor_abc_noduplicates$p),]
cnvr_top4_abc <- cnvr_top4_abc[grepl("CNVR", cnvr_top4_abc$featureB),]
cnvr_top4_abc <- cnvr_top4_abc[c(1:6),]

feat_abc <- rbind(cnvr_top4_abc, gnab_top4_abc)

# combine cnv and mutations and edit feat to match fm
feat_all_abc_original <- rbind(feat_abc, feat_original[!feat_original$featureB%in%feat_abc$featureB,]) # combine top 4 features from all samples and ABC only and remove duplicates

feat_all_abc <- feat_all_abc_original %>%
  mutate(featureB = gsub("AMP:GAIN", "AMP_GAIN", 
                         gsub(":DEL:LOSS", ":DEL_LOSS",
                              gsub("_GAIN", ":GAIN",
                                   gsub("_LOSS", ":LOSS",
                                        gsub("_AMP", ":AMP",
                                             gsub("_DEL", ":DEL",
                                                  gsub("SV_BCL2", "BCL2:SV", featureB))))))))

mat <- fm_cyt_order[as.character(feat_all_abc[, "featureB"]), coo=="ABC"]

mat <- as.matrix(mat)

## change codes for CNVs 
mat[grepl(":AMP|:GAIN", rownames(mat)),][mat[grepl(":AMP|:GAIN", rownames(mat)),]==1] <- "Amplification"
mat[grepl(":LOSS|:DEL", rownames(mat)),][mat[grepl(":LOSS|:DEL", rownames(mat)),]==1] <- "Deletion"
mat[grepl(":SV", rownames(mat)),][mat[grepl(":SV", rownames(mat)),]==1] <- "Structural variation"
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


# data frame of sample annotations
annot <- droplevels(data.frame(COO_byGEP = coo[coo=="ABC"],
                               IPI = as.factor(fm_cyt_order["N:CLIN:IPI", coo == "ABC"])))
levels(annot$IPI) <- c("NA", "0", "1", "2", "3", "4")
annot[is.na(annot)] <- "NA"

# create heatmap annotations
ha1 <- HeatmapAnnotation(Cytolytic_score = anno_barplot(as.numeric(fm_cyt_order["N:SAMP:CytolyticScore", coo == "ABC"]), 
                                                  bar_width = 0.75, 
                                                  border = FALSE, 
                                                  axis = TRUE,
                                                  axis_param = list(gp = gpar(fontsize = 5, lwd = 0.5)),
                                                  gp = gpar(col = NA, fill = "grey50"),
                                                  ylim = c(3.75,11)),
                         height = unit(0.5, "cm"),
                         show_annotation_name = F
)


ha2 = HeatmapAnnotation(Mutations = anno_barplot(as.numeric(fm_cyt_order["N:SAMP:numberOfMutations", coo == "ABC"]), 
                                           bar_width = 0.75, 
                                           border = FALSE, 
                                           axis = TRUE,
                                           axis_param = list(gp = gpar(fontsize = 3, lwd = 0.5)),
                                           gp = gpar(col = NA, fill = "black")
),
CNVs = anno_barplot(as.numeric(fm_cyt_order["N:SAMP:numberOfCNAs", coo == "ABC"]), 
                   bar_width = 0.75, 
                   border = FALSE, 
                   axis = TRUE,
                   axis_param = list(gp = gpar(fontsize = 3, lwd = 0.5)),
                   gp = gpar(col = NA, fill = "black")
),
Purity = anno_barplot(as.numeric(fm_cyt_order["N:SAMP:purity_absolute_reviewed", coo == "ABC"]), 
                   bar_width = 0.75, 
                   border = FALSE, 
                   axis = TRUE,
                   axis_param = list(gp = gpar(fontsize = 3, lwd = 0.5)),
                   gp = gpar(col = NA, fill = "black")
),
df = annot,
col = list(COO_byGEP = structure(brewer.pal(length(unique(annot$COO_byGEP)), "Set1")[1], 
                                 names = as.character(unique(annot$COO_byGEP))),
           IPI = structure(gsub("#FC8D59", "grey", rev(brewer.pal(7, "Spectral"))[c(1:6)]), 
                           names = sort(as.character(unique(annot$IPI))))
),
annotation_height = c(2, 2, 3, 3, 3),
gap = unit(0.5, "mm"),
annotation_legend_param = list(COO_byGEP = list(title = "Subtype", title_gp = gpar(fontsize = 5), 
                                                labels_gp = gpar(fontsize = 5), grid_height = unit(0.2, "cm"), grid_width = unit(2, "mm")),
                               IPI = list(title = "IPI", title_gp = gpar(fontsize = 5), 
                                          labels_gp = gpar(fontsize = 5), grid_height = unit(0.2, "cm"), grid_width = unit(2, "mm"))
),
show_annotation_name = F,
height = unit(1.5, "cm")
)


cor_adj_p_abc <- cor_abc_noduplicates$FDR[match(feat_all_abc_original$featureB, cor_abc_noduplicates$featureB)]

zero_col_mat = matrix(nrow = nrow(mat), ncol = 0)
rownames(zero_col_mat) <- cor_adj_p_abc


cor_p_abc <- cor_abc$p[match(feat_all_abc_original$featureB, cor_abc_noduplicates$featureB)]
cor_p_abc <- substr(cor_p_abc, 1,5)

zero_col_mat_2 = matrix(nrow = nrow(mat), ncol = 0)
rownames(zero_col_mat_2) <- cor_p_abc


# make heatmap

sampleorder <- order(c(rep(1,6), rep(2,6), rep(3,nrow(mat)-12)),cor_p_abc) # order samples first by CNV/SV, Mutation, non-top ABC and then by p value

ht <- Heatmap(mat[sampleorder,],
              col = c("0" = "#EDEDED", "Nonsynonymous mutation" = "black", "Amplification" = "#DC0000FF", "Deletion" = "#3C5488FF", "Structural variation" = "#00A087FF"), 
              rect_gp = gpar(col= "white", lwd = unit(0.5, "mm")),
              top_annotation = ha1, 
              #top_annotation_height = unit(0.5, "cm"), 
              bottom_annotation = ha2,
              #bottom_annotation_height = unit(1.5, "cm"),
              row_names_side = "left",
              row_names_gp = gpar(fontsize = 7),#, fontface = "italic"),
              show_column_names = FALSE, 
              show_column_dend = FALSE, 
              cluster_columns = FALSE,
              cluster_rows = FALSE,
              split = c(rep("Top CNV/SV\n(ABC)",6), rep("Top mutation\n(ABC)",6), rep("Top mutation/\nCNV/SV\nAll DLBCL", nrow(mat)-12)),
              row_title_gp = gpar(fontsize = 5),
              show_heatmap_legend = TRUE,
              heatmap_legend_param = list(title = "Alteration",
                                          title_gp = gpar(fontsize = 5),
                                          labels_gp = gpar(fontsize = 5),
                                          grid_height = unit(0.2, "cm"),
                                          grid_width = unit(2, "mm"))
) +
  Heatmap(zero_col_mat_2[sampleorder,],
          row_names_gp = gpar(fontsize = 5)) +
  Heatmap(zero_col_mat[sampleorder,],
          row_names_gp = gpar(fontsize = 5)
  )


# oncoprint
pdf("FigureS3B_DLBCL_ABC_cytscore_oncoprint.pdf", height = 2.75, width = 4.75)
draw(ht, padding = unit(c(2, 15, 2, 2), "mm"))
decorate_annotation("Cytolytic_score", {
  grid.text("Cytolytic score", unit(0, "npc") - unit(6, "mm"), 0.5, 
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






## ----------------------



## GCB

cor_gcb_noduplicates <- cor_gcb[!grepl("6P_GAIN|7Q_GAIN|7P_AMP$|11Q_GAIN", cor_gcb$featureB),] # not useful 

gnab_top4_gcb <- cor_gcb_noduplicates[order(cor_gcb_noduplicates$p),]
gnab_top4_gcb <- gnab_top4_gcb[grepl("GNAB", gnab_top4_gcb$featureB),]
gnab_top4_gcb <- gnab_top4_gcb[c(1:6),]

cnvr_top4_gcb <- cor_noduplicates[order(cor_noduplicates$p),]
cnvr_top4_gcb <- cnvr_top4_gcb[grepl("CNVR", cnvr_top4_gcb$featureB),]
cnvr_top4_gcb <- cnvr_top4_gcb[c(1:6),]

feat_gcb <- rbind(cnvr_top4_gcb, gnab_top4_gcb)

feat_all_gcb_original <- rbind(feat_gcb, feat_original[!feat_original$featureB%in%feat_gcb$featureB,]) # combine top 4 features from all samples and gcb only and remove duplicates

feat_all_gcb <- feat_all_gcb_original %>%
  mutate(featureB = gsub("AMP:GAIN", "AMP_GAIN", 
                         gsub("_GAIN", ":GAIN",
                              gsub("_LOSS", ":LOSS",
                                   gsub("_AMP", ":AMP",
                                        gsub("_DEL", "_DEL",
                                             gsub("SV_BCL2", "BCL2:SV", featureB)))))))

mat <- fm_cyt_order[as.character(feat_all_gcb[, "featureB"]), coo=="GCB"]

mat <- as.matrix(mat)

## change codes for CNVs 
mat[grepl(":AMP|:GAIN", rownames(mat)),][mat[grepl(":AMP|:GAIN", rownames(mat)),]==1] <- "Amplification"
mat[grepl(":LOSS|:DEL", rownames(mat)),][mat[grepl(":LOSS|:DEL", rownames(mat)),]==1] <- "Deletion"
mat[grepl(":SV", rownames(mat)),][mat[grepl(":SV", rownames(mat)),]==1] <- "Structural variation"
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

## create data frames for heatmap annotations

# data frame of sample annotations
annot <- droplevels(data.frame(COO_byGEP = coo[coo == "GCB"],
                               IPI = as.factor(fm_cyt_order["N:CLIN:IPI", coo == "GCB"])))
levels(annot$IPI) <- c("NA", "0", "1", "2", "3", "4", "5")
annot[is.na(annot)] <- "NA"


## create heatmap annotations
ha1 <- HeatmapAnnotation(Cytolytic_score = anno_barplot(as.numeric(fm_cyt_order["N:SAMP:CytolyticScore", coo=="GCB"]), 
                                                  bar_width = 0.75, 
                                                  border = FALSE, 
                                                  axis = TRUE,
                                                  axis_param = list(gp = gpar(fontsize = 5, lwd = 0.5)),
                                                  gp = gpar(col = NA, fill = "grey50"),
                                                  ylim = c(3.75,11)),
                         show_annotation_name = F,
                         height = unit(0.5, "cm"))


ha2 = HeatmapAnnotation(Mutations = anno_barplot(as.numeric(fm_cyt_order["N:SAMP:numberOfMutations", coo == "GCB"]), 
                                           bar_width = 0.75, 
                                           border = FALSE, 
                                           axis = TRUE,
                                           axis_param = list(gp = gpar(fontsize = 3, lwd = 0.5)),
                                           gp = gpar(col = NA, fill = "black")
),
CNVs = anno_barplot(as.numeric(fm_cyt_order["N:SAMP:numberOfCNAs", coo == "GCB"]), 
                   bar_width = 0.75, 
                   border = FALSE, 
                   axis = TRUE,
                   axis_param = list(gp = gpar(fontsize = 3, lwd = 0.5)),
                   gp = gpar(col = NA, fill = "black")
),
Purity = anno_barplot(as.numeric(fm_cyt_order["N:SAMP:purity_absolute_reviewed", coo == "GCB"]), 
                   bar_width = 0.75, 
                   border = FALSE, 
                   axis = TRUE,
                   axis_param = list(gp = gpar(fontsize = 3, lwd = 0.5)),
                   gp = gpar(col = NA, fill = "black")
),
df = annot,
col = list(COO_byGEP = structure(brewer.pal(length(unique(annot$COO_byGEP)), "Set1")[3], 
                                 names = as.character(unique(annot$COO_byGEP))),
           IPI = structure(gsub("#D53E4F", "grey", rev(brewer.pal(length(unique(annot$IPI)), "Spectral"))), 
                           names = sort(as.character(unique(annot$IPI))))
),
annotation_height = c(2, 2, 3, 3, 3),
gap = unit(0.5, "mm"),
annotation_legend_param = list(COO_byGEP = list(title = "Subtype", title_gp = gpar(fontsize = 5), 
                                                labels_gp = gpar(fontsize = 5), grid_height = unit(0.2, "cm"), grid_width = unit(2, "mm")),
                               IPI = list(title = "IPI", title_gp = gpar(fontsize = 5), 
                                          labels_gp = gpar(fontsize = 5), grid_height = unit(0.2, "cm"), grid_width = unit(2, "mm"))
),
show_annotation_name = F,
height = unit(1.5, "cm")
)


cor_adj_p_gcb <- cor_gcb_noduplicates$FDR[match(feat_all_gcb_original$featureB, cor_gcb_noduplicates$featureB)]

zero_col_mat = matrix(nrow = nrow(mat), ncol = 0)
rownames(zero_col_mat) <- cor_adj_p_gcb


cor_p_gcb <- cor_gcb_noduplicates$p[match(feat_all_gcb_original$featureB, cor_gcb_noduplicates$featureB)]
cor_p_gcb <- substr(cor_p_gcb, 1,5)

zero_col_mat_2 = matrix(nrow = nrow(mat), ncol = 0)
rownames(zero_col_mat_2) <- cor_p_gcb


# make heatmap

sampleorder <- order(c(rep(1,6), rep(2,6), rep(3,nrow(mat)-12)),cor_p_gcb) # order samples first by CNV/SV, Mutation, non-top ABC and then by p value

ht <- Heatmap(mat[sampleorder,],
              col = c("0" = "#EDEDED", "Nonsynonymous mutation" = "black", "Amplification" = "#DC0000FF", "Deletion" = "#3C5488FF", "Structural variation" = "#00A087FF"), 
              rect_gp = gpar(col= "white", lwd = unit(0.5, "mm")),
              top_annotation = ha1, 
              #top_annotation_height = unit(0.5, "cm"), 
              bottom_annotation = ha2,
              #bottom_annotation_height = unit(1.5, "cm"),
              row_names_side = "left",
              row_names_gp = gpar(fontsize = 7),#, fontface = "italic"),
              show_column_names = FALSE, 
              show_column_dend = FALSE, 
              cluster_columns = FALSE,
              cluster_rows = FALSE,
              split = c(rep("Top CNV/SV\n(GCB)",6), rep("Top mutation\n(GCB)",6), rep("Top mutation/\nCNV/SV\nAll DLBCL", length(rownames(mat))-12)),
              row_title_gp = gpar(fontsize = 5),
              show_heatmap_legend = TRUE,
              heatmap_legend_param = list(title = "Alteration",
                                          title_gp = gpar(fontsize = 5),
                                          labels_gp = gpar(fontsize = 5),
                                          grid_height = unit(0.2, "cm"),
                                          grid_width = unit(2, "mm"))
) +
  Heatmap(zero_col_mat_2[sampleorder,],
          row_names_gp = gpar(fontsize = 5)) +
  Heatmap(zero_col_mat[sampleorder,],
          row_names_gp = gpar(fontsize = 5)
  )


# oncoprint
pdf("FigureS3C_DLBCL_GCB_cytscore_oncoprint.pdf", height = 2.75, width = 4.5)
draw(ht, padding = unit(c(2, 15, 2, 2), "mm"))
decorate_annotation("Cytolytic_score", {
  grid.text("Cytolytic score", unit(0, "npc") - unit(6, "mm"), 0.5, 
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




