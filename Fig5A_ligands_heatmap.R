
# Plot Figure 5A immunomodulatory molecule heatmap

library(ComplexHeatmap)
library(circlize)
library(dendextend)
library(readxl)
library(RColorBrewer)
library(matrixStats)

# load Hemap gene expression data
data = get(load("data9544_with_gene_symbols.RData"))

# load Hemap annotations
annot = get(load("Hemap_immunology_Annotations.Rdata"))

annot = annot[!(annot$Category.specifying.lineage.tumor.origin=="MM"&annot$CELLS_SORTED==0),] # remove non-malignant myeloma samples

data = data[annot[,1],]

# modify annotations
annot$Category.specifying.subtype.Heatmap <- annot$plotNormals

annot$Category.specifying.subtype.Heatmap[annot$Sample.type%in%c("Cancer", "Prolif")] <- annot$Category.specifying.subtype[annot$Sample.type%in%c("Cancer", "Prolif")]

annot$Category.specifying.subtype.Heatmap[annot$Category.specifying.lineage.tumor.origin=="MDS"] <- "MDS"
annot$Category.specifying.subtype.Heatmap[annot$Category.specifying.lineage.tumor.origin=="AML"] <- "AML"
annot$Category.specifying.subtype.Heatmap[annot$Category.specifying.subtype=="AILT"] <- "AITL"
annot$Category.specifying.subtype.Heatmap[annot$Category.specifying.subtype=="LC"] <- "LCH"


# -----------------------


# list of immunomodulatory genes
ligands = read_excel("ligands.xlsx")
genelist = ligands$Gene

# Select cancer samples and summarize by cancer type medians
ligands_cancer <- data[annot$Sample.type%in%c("Cancer", "Prolif")&annot$Category.specifying.subtype.Heatmap!="na", colnames(data) %in% genelist]
ligands_cancer_median <- aggregate(ligands_cancer, by = list(annot$Category.specifying.subtype.Heatmap[annot$Sample.type%in%c("Cancer", "Prolif")&annot$Category.specifying.subtype.Heatmap!="na"]), FUN = "median")
rownames(ligands_cancer_median) <- ligands_cancer_median$Group.1
ligands_cancer_median$Group.1 <- NULL

# -----------------------

# select genes from expression matrix
mat_unscaled = t(ligands_cancer_median)
mat_unscaled = mat_unscaled[genelist,]

# scaale
mat = t(apply(mat_unscaled, 1, scale))
colnames(mat) <- colnames(mat_unscaled)

# remove genes with max log2 expression <5 (not expressed)
mat = mat[!rowMaxs(mat_unscaled)<5,]

# make another data matrix for receptor names for later use
mat_receptor <- mat

# add more common names of ligands in parentheses
ligands <- ligands[match(rownames(mat), ligands$Gene),]

rownames(mat) <- paste(ligands$Gene, paste0("(", ligands$`Common name`, ")"))
rownames(mat) <- gsub("\\ \\(NA\\)", "", rownames(mat))
rownames(mat) <- gsub("PD-1H", "", rownames(mat))

# make lineage annotation to heatmap
lineage <- colnames(mat)
lineage[grepl("AML|CML|MDS|JMML|Monocyte|Neutrophil|Myeloid|Dendritic|Macrophage|^LC|Langerhans|Erythroid", lineage)] <- "Myeloid"
lineage[grepl("CLL|HCL|FL|MM|MALT|MCL|MZL|CHL|NLPHL|BL|B-ALL|DLBCL|B cell|Germinal|Plasma", lineage)] <- "B cell"
lineage[grepl("AITL|ALCL|ATL|HSTCL|ENKTL|CTCL|PTCLNOS|T-ALL|T/NKcell", lineage)] <- "T/NK cell"

lineage <- as.data.frame(lineage)

mainclass <- colnames(mat)
mainclass[grepl("AML|CML|JMML|T-ALL|CLL|pre-B-ALL|HCL", mainclass)] <- "Leukemia"
mainclass[grepl("DLBCL|FL|MALT|MCL|MZL|CHL|NLPHL|BL|AI|ALCL|ATL|HSTCL|ENKTL|CTCL|PTCLNOS", mainclass)] <- "Lymphoma"
mainclass[grepl("MM|LC|MDS", mainclass)] <- "Other"
mainclass <- as.data.frame(mainclass)

heatmap_annot <- cbind(lineage, mainclass)

# heatmap annotations
ha1 = HeatmapAnnotation(df = heatmap_annot,
                        col = list(lineage = structure(brewer.pal(length(unique(heatmap_annot$lineage)), "Set2")[c(1,3,2)],
                                                       names = as.character(unique(heatmap_annot$lineage))[c(2,1,3)]),
                                   mainclass = structure(brewer.pal(10, "Set3")[c(10,7,9)],
                                                         names = as.character(unique(heatmap_annot$mainclass))[c(2,1,3)])
                        ),
                        annotation_legend_param = list(lineage = list(title = "Lineage", title_gp = gpar(fontsize = 10, fontface = "plain"),
                                                                      labels_gp = gpar(fontsize = 10)),
                                                       mainclass = list(title = "Cancer type", title_gp = gpar(fontsize = 10, fontface = "plain"),
                                                                        labels_gp = gpar(fontsize = 10))
                        ))


# ligand cell type information annotation
ha2 = rowAnnotation(df = as.data.frame(ligands[5]),
                    col = list(`Immune checkpoint function` = structure(brewer.pal(length(unique(ligands$`Immune checkpoint function`)), "Paired")[c(3,2,1,4)],
                                                                        names = as.character(unique(ligands$`Immune checkpoint function`))[c(3,2,1,4)])),
                    annotation_legend_param = list(`Immune checkpoint function` = list(title = "Immune checkpoint function", title_gp = gpar(fontsize = 10, fontface = "plain"),
                                                                                       labels_gp = gpar(fontsize = 10)), width = unit(0.25, "cm")))


# receptor names
zero_col_mat = matrix(nrow = nrow(mat), ncol = 0)
rownames(zero_col_mat) <- ligands$`Receptor common name`
rownames(zero_col_mat) <- gsub("\\, TNFSF14", "", rownames(zero_col_mat))

# dendrogram colors
dend = hclust(as.dist(1-cor(mat, method="spearman")), method = "ward.D")
dend = color_branches(dend, k = 3, col = c("#66C2A5", "#B3DE69", "#FC8D62"))


# print heatmap
pdf(file = "Figure5A_Hemap_ligands_heatmap.pdf", height = 8.5, width = 10)
draw(Heatmap(mat,
        split = ligands$Split,
        gap = unit(2, "mm"),
        name = "Expression (log2)\nZ-score",
        col = colorRamp2(c(seq(-2, 2, length.out = 9), 4), rev(brewer.pal(9, "RdBu"))[c(1:9,9)]),
        rect_gp = gpar(col = "white", lty = 1, lwd = 1),
        row_names_side = "right",
        row_names_gp = gpar(fontsize = 9),
        column_names_gp = gpar(fontsize = 10),
        cluster_rows = FALSE,
        clustering_distance_rows = "spearman",
        clustering_method_rows = "ward.D",
        cluster_columns = dend,
        clustering_distance_columns = "spearman",
        clustering_method_columns = "ward.D",
        top_annotation = ha1,
        heatmap_legend_param = list(title_gp = gpar(fontsize = 10, fontface = "plain"))
) +
  Heatmap(zero_col_mat,
          right_annotation = ha2,
          row_names_gp = gpar(fontsize = 9)),
auto_adjust = F
)
dev.off()



