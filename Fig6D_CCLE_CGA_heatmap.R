
# Plot heatmap of cancer-germline antuigen expression and methylation in CCLE hematologic cell line data (Figure 6D)

library(CePa)
library(biomaRt)
library(ggplot2)
library(cowplot)
library(reshape2)
library(matrixStats)
library(data.table)
library(ggrepel)
library(dplyr)
library(gridExtra)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

# load data
meth <- fread("CCLE_RRBS_TSS_1kb_20180614.txt", data.table = F)
meth_cpg <- fread("CCLE_RRBS_TSS_CpG_clusters_20180614.txt", data.table = F)
rpkm <- read.gct(file = "CCLE_DepMap_18q3_RNAseq_RPKM_20180718.gct") # read in CCLE RNA-seq rpkm table
annot <- fread(file = "DepMap-2018q3-celllines.csv", data.table = F)
annot2 <- fread(file = "CCLE_sample_info_file_2012-10-18_modified.csv", data.table = F, check.names = T) # read in CCLE annotation file (2012 publication) with lineage data added by Olli

# clean rpkm column names
colnames(rpkm) <- gsub("..ACH.*", "", colnames(rpkm))

# subset rpkm to methylation data
samples_rkpm_meth <- colnames(rpkm)[colnames(rpkm) %in% colnames(meth_cpg)]

# subset annotations to data
annot <- annot[annot$CCLE_Name%in%samples_rkpm_meth,]

# select hematopoietic cell lines from annot
annot_hem <- annot[grepl("HAEMATOPOIETIC", annot$CCLE_Name),]


colnames(meth) <- gsub("_name|TSS_|cluster_", "", colnames(meth_cpg))
colnames(meth_cpg) <- gsub("_name|TSS_|cluster_", "", colnames(meth_cpg))

# rpkm data frame with hematopoietic cell lines
rpkm_hem <- rpkm[,as.character(annot_hem$CCLE_Name)]

# combine TSS and CpG methylation data
colnames(meth[,c("id", "gene", as.character(annot_hem$CCLE_Name))])==colnames(meth_cpg[,c("id", "gene", as.character(annot_hem$CCLE_Name))]) # check that column names match before joining
meth_hem <- rbind(meth[,c("id", "gene", as.character(annot_hem$CCLE_Name))], meth_cpg[,c("id", "gene", as.character(annot_hem$CCLE_Name))])
meth_hem[,!colnames(meth_hem)%in%c("id", "gene")] <- sapply(meth_hem[,!colnames(meth_hem)%in%c("id", "gene")], as.numeric) # make numeric


# get gene symbols
ensembl_hs_mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
ensembl_df <- getBM(attributes=c("ensembl_gene_id", "ensembl_gene_id_version",
                                 "hgnc_symbol", "description", "chromosome_name", "start_position"),
                    mart=ensembl_hs_mart)
gene_annot <- ensembl_df[match(gsub("\\..*", "", rownames(rpkm)), ensembl_df$ensembl_gene_id),] # match biomart data using Ensembl gene symbols
rownames(gene_annot) <- rownames(rpkm)

# get cga genes:
genelist <- as.character(read.table("cga.txt")[,1])

gene_annot_cga <- gene_annot[gene_annot$hgnc_symbol%in%gsub("C11orf85", "MAJIN", genelist),]

rpkm_hem_cga <- rpkm_hem
rownames(rpkm_hem_cga) <- gsub("\\..*", "", rownames(rpkm))
rpkm_hem_cga <- rpkm_hem_cga[gene_annot_cga$ensembl_gene_id,]
rownames(rpkm_hem_cga) <- gsub("MAJIN", "C11orf85", gene_annot_cga$hgnc_symbol)
rpkm_hem_cga_log2 <- log2(rpkm_hem_cga + 0.25)

meth_hem_cga <- as.data.frame(meth_hem[meth_hem$gene%in%gene_annot_cga$hgnc_symbol,])
rownames(meth_hem_cga) <- gsub("MAJIN", "C11orf85", meth_hem_cga$id)
meth_hem_cga$id <- NULL
meth_hem_cga$gene <- NULL

hem_cga <- t(rbind(rpkm_hem_cga_log2, meth_hem_cga))
hem_cga <- merge(hem_cga, annot_hem, by.x = "row.names", by.y = "CCLE_Name")
hem_cga <- merge(hem_cga, annot2, by.x = "Row.names", by.y = "CCLE.name")
hem_cga <- hem_cga[grepl("lymphoma|leukaemia|myeloma", hem_cga$Hist.Subtype1),]
colnames(hem_cga)[1] <- "CCLE_Name"

# add shortened annotations
grep("lymphoma|leukaemia|myeloma", hem_cga$Hist.Subtype1, value=T)
hem_cga$Hist.Subtype_hem[hem_cga$Hist.Subtype1%in%"plasma_cell_myeloma"]="MM"
hem_cga$Hist.Subtype_hem[hem_cga$Hist.Subtype1%in%"mantle_cell_lymphoma"]="MCL"
hem_cga$Hist.Subtype_hem[hem_cga$Hist.Subtype1%in%"diffuse_large_B_cell_lymphoma"]="DLBCL"
hem_cga$Hist.Subtype_hem[hem_cga$Hist.Subtype1%in%"chronic_lymphocytic_leukaemia-small_lymphocytic_lymphoma"]="CLL"
hem_cga$Hist.Subtype_hem[hem_cga$Hist.Subtype1%in%"blast_phase_chronic_myeloid_leukaemia"]="CML"
hem_cga$Hist.Subtype_hem[hem_cga$Hist.Subtype1%in%"anaplastic_large_cell_lymphoma"]="ALCL"
hem_cga$Hist.Subtype_hem[hem_cga$Hist.Subtype1%in%"adult_T_cell_lymphoma-leukaemia"]="TCL other"
hem_cga$Hist.Subtype_hem[hem_cga$Hist.Subtype1%in%"acute_myeloid_leukaemia"]="AML"
hem_cga$Hist.Subtype_hem[hem_cga$Hist.Subtype1%in%"acute_lymphoblastic_T_cell_leukaemia"]="T-ALL"
hem_cga$Hist.Subtype_hem[hem_cga$Hist.Subtype1%in%"acute_lymphoblastic_B_cell_leukaemia"]="pre-B-ALL"
hem_cga$Hist.Subtype_hem[hem_cga$Hist.Subtype1%in%"Hodgkin_lymphoma"]="CHL"
hem_cga$Hist.Subtype_hem[hem_cga$Hist.Subtype1%in%"Burkitt_lymphoma"]="BL"
hem_cga$Hist.Subtype_hem[hem_cga$Hist.Subtype1%in%c("B_cell_lymphoma_unspecified")]="BCL other"
hem_cga$Hist.Subtype_hem[hem_cga$Hist.Subtype1%in%c("peripheral_T_cell_lymphoma_unspecified")]="TCL other"

# add lineage annotations
hem_cga$Lineage <- ifelse(grepl("ALCL|T-ALL|TCL", hem_cga$Hist.Subtype_hem), "T/NK cell", 
                          ifelse(grepl("AML|CML", hem_cga$Hist.Subtype_hem), "Myeloid", "B cell"))


# add methylation means to data frame for plotting
meth_mean <- function(gene){
  df <- data.frame(gene = rowMeans(data.frame(hem_cga[,colnames(hem_cga)[grepl(paste0(gene, "_"), colnames(hem_cga))]]), na.rm = T))
  return(df)
}

cga_meth_mean <- do.call(cbind, lapply(gsub("MAJIN", "C11orf85", genelist), meth_mean))
colnames(cga_meth_mean) <- paste0(gsub("MAJIN", "C11orf85", genelist), "_meth_mean")

hem_cga_meth_mean <- cbind(hem_cga, cga_meth_mean)
hem_cga_meth_mean <- hem_cga_meth_mean[order(hem_cga_meth_mean$Hist.Subtype_hem),]

# order cancers for plotting
hem_cga_meth_mean$Hist.Subtype_hem <- factor(hem_cga_meth_mean$Hist.Subtype_hem, levels = c("T-ALL", "pre-B-ALL", "AML", "CML", "CLL", "MM", "BL", "CHL", "DLBCL", "MCL", "BCL other", "ALCL", "TCL other"))

# data prepared

## --------------------------------------------------------------------------------------


## Test correlations for all methylation values averaged per gene

# function to correlate gexp to averaged same gene methylation features
corr_mean <- function(gene) {
  cor.result <- cor.test(hem_cga[,gene],
                         as.numeric( # numeric vector for cor.test
                           rowMeans( # average methylation values
                             data.frame( # data frame for rowMeans
                               hem_cga[,colnames(hem_cga)[grepl(paste0(gene, "_"), colnames(hem_cga))]]), na.rm = T)), method = "spearman") # grep columns with underscore to get methylation values only
  p <- cor.result$p.value
  r <- cor.result$estimate
  return(data.frame(gexp = gene, meth = gene, r = r, p = p))
}

# select genes that have methylation data
genes = genelist[genelist%in%gsub("_.*", "", colnames(hem_cga)[grepl("_", colnames(hem_cga))])]

# apply function on gene list
result_mean <- do.call(rbind, lapply(genes, corr_mean))
result_mean$q <- p.adjust(result_mean$p, method = "fdr")

# check results sorted by siginificance
result_mean %>%
  arrange(q)

# write sorted mean methylation correlations
result_mean_sorted <- result_mean %>%
  arrange(r)

# add significance codes to mean gexp to meth correlation for each gene
mean_corr <- result_mean %>%
  mutate(signifCode = ifelse(q<0.0001, "****",
                             ifelse(q<0.001, "***", 
                                    ifelse(q<0.01, "**",
                                           ifelse(q<0.05, "*", "")))))


## -------------------------------------------------------------------------------

# Heatmap 

mat <- data.matrix(rpkm_hem_cga_log2[,hem_cga_meth_mean$CCLE_Name])
mat <- t(apply(mat, 1, scale))
colnames(mat) <- gsub("_.*", "", colnames(rpkm_hem_cga_log2[,hem_cga_meth_mean$CCLE_Name]))

mat <- mat[gene_annot_cga$hgnc_symbol,]

mat 

# colors
cols <- read.table("colors_hemap_immunology.tsv", header = TRUE, sep = "\t", comment.char = " ")
cols <- cols %>% mutate(combinedcolor = ifelse(subtypecolor!="", as.character(subtypecolor), as.character(color)))
cols$sample <- gsub("^BCL", "BCL other", gsub("TCL", "TCL other", cols$sample))

# data frames for heatmap annotations
exprs_num <- colSums(rpkm_hem_cga[,hem_cga_meth_mean$CCLE_Name]>0.5)

meth_mat <- hem_cga_meth_mean[,grep("mean", colnames(hem_cga_meth_mean))]
colnames(meth_mat) <- gsub("_meth_mean", "", colnames(meth_mat))
meth_mat <- t(meth_mat)
colnames(meth_mat) <- gsub("_.*", "", hem_cga_meth_mean$CCLE_Name)
meth_mat <- meth_mat[,match(colnames(mat), colnames(meth_mat))]
meth_mat[meth_mat>0.5] = 1
meth_mat[meth_mat<0.5] = 0
meth_num <- colSums(meth_mat==0, na.rm = T)

disease_order <- data.frame(disease = unique(hem_cga_meth_mean$Hist.Subtype_hem), rank = c(13, 3, 11, 7, 8, 5, 4, 9, 10, 6, 2, 1, 12))
disease_order$disease <- as.character(disease_order$disease)
hem_cga_meth_mean <- merge(hem_cga_meth_mean, disease_order, by.x = "Hist.Subtype_hem", by.y = "disease")

sample_order <- order(hem_cga_meth_mean$rank, exprs_num)
exprs_num <- exprs_num[sample_order]
meth_num <- meth_num[sample_order]
annot_df <- data.frame(subtype = hem_cga_meth_mean$Hist.Subtype_hem[sample_order])

# heatmap annotations

ha1 = HeatmapAnnotation(exprs = anno_barplot(exprs_num, 
                                             bar_width = 0.75, 
                                             border = FALSE, 
                                             axis = TRUE,
                                             axis_param = list(gp = gpar(fontsize = 5, lwd = 0.5)),
                                             gp = gpar(col = NA, fill = "grey50")),
                        meth = anno_barplot(meth_num, 
                                            bar_width = 0.75, 
                                            border = FALSE, 
                                            axis = TRUE,
                                            axis_param = list(gp = gpar(fontsize = 5, lwd = 0.5)),
                                            gp = gpar(col = NA, fill = "#a52a2a")),
                        df = annot_df, 
                        col = list(subtype = structure(as.character(cols$combinedcolor[cols$sample %in% hem_cga_meth_mean$Hist.Subtype_hem]), 
                                                       names = as.character(cols$sample[cols$sample %in% hem_cga_meth_mean$Hist.Subtype_hem]))[as.character(disease_order$disease[order(disease_order$rank)])]
                        ),
                        annotation_legend_param = list(subtype = list(title = "Cancer type", title_gp = gpar(fontsize = 7, fontface = "bold"), 
                                                                      labels_gp = gpar(fontsize = 7), grid_height = unit(3, "mm"), grid_width = unit(3, "mm"))
                        ),
                        show_annotation_name = F,
                        height = unit(2, "cm"))

meth_corr_df <- data.frame(corr = mean_corr$r[match(rownames(mat), mean_corr$gexp)])
colnames(meth_corr_df) <- "Correlation to methylation"

mat <- mat[order(meth_corr_df$`Correlation to methylation`),sample_order]

meth_corr_df <- data.frame(corr = mean_corr$r[match(rownames(mat), mean_corr$gexp)])
colnames(meth_corr_df) <- "meth_corr"

ha2 = HeatmapAnnotation(df = meth_corr_df,
                        which = "row",
                        col = list(meth_corr = colorRamp2(c(-1, -0.5, -0.375, -0.25, -0.125, 0, 0.125, 0.25, 0.375, 0.5, 1), rev(brewer.pal(11, "PuOr")))),
                        annotation_legend_param = list(meth_corr = list(title = "Correlation to methylation", title_gp = gpar(fontsize = 7, fontface = "bold"), 
                                                                        labels_gp = gpar(fontsize = 7)), width = unit(0.1, "cm"), grid_height = unit(2, "mm"), grid_width = unit(2, "mm"), legend_direction = "horizontal", legend_width = unit(2, "cm"), title_position = "topcenter"))


# print heatmap
pdf(file = "Figure6D_CCLE_CGA_heatmap.pdf", height = 4, width = 5)
draw(
       ha2 +
       Heatmap(mat,
               name = "Expression Z-score", 
               col = colorRamp2(seq(-3, 3, length.out = 11), rev(brewer.pal(11, "RdBu"))),
               row_names_side = "right", 
               row_names_gp = gpar(fontsize = 6),
               show_column_names = FALSE,
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               top_annotation = ha1, 
               heatmap_legend_param = list(title_gp = gpar(fontsize = 7, fontface = "bold"),
                                           labels_gp = gpar(fontsize = 7),
                                           legend_direction = "horizontal",
                                           legend_width = unit(2, "cm"), title_position = "topcenter",
                                           grid_height = unit(2, "mm"),
                                           grid_width = unit(2, "mm")
               )
       ),
     auto_adjust = F,
     padding = unit(c(2, 15, 2, 2), "mm"), heatmap_legend_side = "bottom")
decorate_annotation("exprs", {
  grid.text("# CGAs expressed", unit(0, "npc") + unit(1, "mm"), 0.5, 
            default.units = "npc", just = c("left", "bottom"), gp = gpar(fontsize = 7))
})
decorate_annotation("meth", {
  grid.text("# CGAs hypomethylated", unit(0, "npc") + unit(1, "mm"), 0.5, 
            default.units = "npc", just = c("left", "bottom"), gp = gpar(fontsize = 7))
})
dev.off()

