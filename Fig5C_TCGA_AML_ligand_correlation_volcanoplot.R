# Plot Figure 5C-D volcano plots of co-stimulatory gene expression vs genetic alterations

library(ggplot2)
library(cowplot)
library(ggrepel)
library(dplyr)
library(ggsci)

# load correlation data
dlbcl <- read.table("TableS5_GSE98588_DLBCL_costim_correlations.tsv", header = TRUE, sep = "\t")
aml_tcga <- read.table("TableS5_TCGA_AML_costim_LR_correlations.tsv", header = TRUE, sep = "\t")
dlbcl_tcga <- read.table("TCGA_DLBCL_costim_LR_correlations.tsv", header = TRUE, sep = "\t")

# load TCGA AML and DLBCL feature matrices
aml=get(load("DUFVA_TCGA_AML_FM_meth.Rdata"))
matrix=get(load("TCGA_DLBCL_FM_DUFVA.Rdata"))

# ------------------------------------------------------------------

# modify DLBCL data frame
dlbcl_gexp_gnab_cnvr <- dlbcl %>%
  mutate(log10.FDR = -log10(FDR),
         featureA_short = gsub("^.......", "",
                               gsub(":::::", "", featureA)),
         featureB_short = gsub("GNAB", "MUT",
                               gsub("_", ".",
                                    gsub("SV_", "SV:",
                                         gsub(".*@|CNVR:|_nonsynonymous", "",
                                              gsub("_DEL", " DEL",
                                                   gsub("_LOSS", " LOSS",
                                                        gsub("_DEL_LOSS", " DEL/LOSS",
                                                             gsub("_AMP", " AMP",
                                                                  gsub("_GAIN", " GAIN",
                                                                       gsub("_AMP_GAIN", " AMP/GAIN",
                                                                            gsub("^..", "",
                                                                                 gsub(":::::", "", featureB)))))))))))),
         alteration_type = ifelse(grepl("METH", featureB), "Methylation (METH)",
                                  ifelse(grepl("CNVR:SV", featureB), "Structural variation (SV)",
                                         ifelse(grepl("AMP|GAIN", featureB), "Amplification",
                                                ifelse(grepl("DEL|LOSS", featureB), "Deletion",
                                                       ifelse(grepl("GNAB", featureB), "Nonsynonymous mutation (MUT)", "")))))
  ) %>%
  filter(grepl("GEXP", featureA)&grepl("GNAB|CNVR", featureB), # select only genetic alterations
         #gsub("^.......", "", gsub(":::::", "", featureA))%in%ligands, # select only ligand gene expression (no receptors)
         !grepl("low_grade", featureB), # remove duplicate CNVR
         !grepl("_SV", featureB) # remove duplicate SV
  )

# ------------------------------------------------------------------

# get amplification/deletion info for TCGA AML

## select samples and CNV data
data <- aml[grepl("N\\:CNVR", rownames(aml)), is.na(aml["C:SAMP:cancermap_cluster",])==FALSE]

mat <- data.matrix(data)

## change codes for alterations
mat[abs(mat<0.3)] = 0
mat[mat>0.3] = 1
mat[mat<(-0.3)] = -1

cnv_sums <- rowSums(mat, na.rm = TRUE)


# modify TCGA AML data frame
aml_tcga_gexp_gnab_cnvr_meth <- aml_tcga %>%
  mutate(log10.FDR = -log10(FDR),
         featureA_short = gsub("^.......", "",
                               gsub(":chr.*", "", featureA)),
         featureB_short = gsub("GNAB", "MUT",
                               gsub("^..", "",
                                    gsub(":chr.*|:mean.*", "", featureB))),
         alteration_type = ifelse(grepl("METH", featureB), "Methylation (METH)",
                                  ifelse(grepl("CNVR", featureB)&cnv_sums[featureB]>0, "Amplification",
                                         ifelse(grepl("CNVR", featureB)&cnv_sums[featureB]<0, "Deletion",
                                                ifelse(grepl("GNAB", featureB), "Nonsynonymous mutation (MUT)", ""))))
  ) %>%
  filter(grepl("GEXP", featureA)&grepl("GNAB|CNVR|METH", featureB), # select only genetic alterations
         !is.na(alteration_type) # remove CNVs not called (< 0.3 absolute log2 ratio)
         #gsub("^.......", "", gsub(":chr.*", "", featureA))%in%ligands # select only ligand gene expression (no receptors)
  ) %>%
  mutate(featureAB_short = paste(featureA_short, featureB_short, sep = "_"))


aml_tcga_gexp_gnab_cnvr_meth[aml_tcga_gexp_gnab_cnvr_meth$log10.FDR=="Inf","log10.FDR"] = 20

# ------------------------------------------------------------------

# get amplification/deletion info for TCGA DLBCL

## select samples and CNV data
data <- matrix[grepl("N\\:CNVR", rownames(matrix)),]

mat <- data.matrix(data)

## change codes for alterations
mat[abs(mat<0.3)] = 0
mat[mat>0.3] = 1
mat[mat<(-0.3)] = -1

cnv_sums <- rowSums(mat, na.rm = TRUE)

# modify TCGA DLBCL data frame
dlbcl_tcga_gexp_gnab_cnvr_meth <- dlbcl_tcga %>%
  filter(grepl("mean|CNVR|GNAB", featureB)) %>% # filter out non-mean methylation features
  mutate(log10.FDR = -log10(FDR),
         featureA_short = gsub("^.......", "",
                               gsub(":chr.*", "", featureA)),
         featureB_short = gsub("GNAB", "MUT",
                               gsub("^..", "",
                                    gsub(":chr.*|:mean.*", "", featureB))),
         alteration_type = ifelse(grepl("METH", featureB), "Methylation (METH)",
                                  ifelse(grepl("CNVR", featureB)&cnv_sums[featureB]>0, "Amplification",
                                         ifelse(grepl("CNVR", featureB)&cnv_sums[featureB]<0, "Deletion",
                                                ifelse(grepl("GNAB", featureB), "Nonsynonymous mutation (MUT)", ""))))
  ) %>%
  filter(grepl("GEXP", featureA)&grepl("GNAB|CNVR|METH", featureB), # select only genetic alterations
         !is.na(alteration_type) # remove CNVs not called (< 0.3 absolute log2 ratio)
         #gsub("^.......", "", gsub(":chr.*", "", featureA))%in%ligands # select only ligand gene expression (no receptors)
  )

# ------------------------------------------------------------------

## plots

p_dlbcl <- ggplot(dlbcl_gexp_gnab_cnvr, aes(x = cor, y = log10.FDR)) +
  geom_point(aes(fill = alteration_type, size = log10.FDR), pch = 21, color = "grey20") +
  geom_label_repel(aes(label=ifelse((cor>0.2&FDR<=0.035)|(cor<(-0.35)&FDR<=0.01), paste(featureA_short, featureB_short, sep = "\n"), '')),
                   point.padding = 0.5,
                   label.padding = 0.15,
                   segment.size = 0.25,
                   #box.padding = 1,
                   size = 3) +
  ylab("-log10(adjusted P value)") +
  xlab("Correlation between expression and alteration") +
  theme(legend.title=element_blank()) +
  guides(size = FALSE,
         color = guide_legend(override.aes = list(size=3))) +
  xlim(c(-0.5, 0.5)) +
  scale_fill_manual(values = c("Amplification" = "#DC0000FF", "Deletion" = "#3C5488FF", "Nonsynonymous mutation (MUT)" = "black", "Structural variation (SV)" = "#00A087FF"))


p_aml_tcga <- ggplot(aml_tcga_gexp_gnab_cnvr_meth, aes(x = cor, y = log10.FDR)) +
  geom_point(aes(fill = alteration_type, size = log10.FDR), pch = 21, color = "grey20") +#, size=2.5) +
  geom_label_repel(aes(label=ifelse((cor>0.3&FDR<=0.005)|(cor<(-0.45)|FDR<=0.0000005), 
                                    ifelse(grepl("METH", featureB_short), featureB_short, paste(featureA_short, featureB_short, sep = "\n")),
                                    '')),
                   point.padding = 0.5,
                   label.padding = 0.15,
                   segment.size = 0.25,
                   #box.padding = 1,
                   size = 3) +
  ylab("-log10(adjusted P value)") +
  xlab("Correlation between expression and alteration") +
  theme(legend.title=element_blank()) +
  guides(size = FALSE,
         color = guide_legend(override.aes = list(size=3))) +
  xlim(c(-0.9, 0.9)) +
  scale_fill_manual(values = c("Methylation (METH)" = "#fea719", "Amplification" = "#DC0000FF", "Deletion" = "#3C5488FF", "Nonsynonymous mutation (MUT)" = "black", "Structural variation (SV)" = "#00A087FF"))


p_dlbcl_tcga <- ggplot(dlbcl_tcga_gexp_gnab_cnvr_meth, aes(x = cor, y = log10.FDR)) +
  geom_point(aes(fill = alteration_type, size = log10.FDR), pch = 21, color = "grey20") +#, size=2.5) +
  geom_label_repel(aes(label=ifelse(FDR<=0.005, 
                                    ifelse(grepl("METH", featureB_short), featureB_short, paste(featureA_short, featureB_short, sep = "\n")),
                                    '')),
                   point.padding = 0.5,
                   label.padding = 0.15,
                   segment.size = 0.25,
                   #box.padding = 1,
                   size = 3) +
  ylab("-log10(adjusted P value)") +
  xlab("Correlation between expression and alteration") +
  theme(legend.title=element_blank()) +
  guides(size = FALSE,
         color = guide_legend(override.aes = list(size=3))) +
  xlim(c(-0.9, 0.9)) +
  scale_fill_manual(values = c("Methylation (METH)" = "#fea719", "Amplification" = "#DC0000FF", "Deletion" = "#3C5488FF", "Nonsynonymous mutation (MUT)" = "black", "Structural variation (SV)" = "#00A087FF"))

# ------------------------------------------------------------------

# print plots
pdf("Figure5D_DLBCL_GSE98588_ligand_correlations_volcanoplot.pdf", height = 4, width = 7)
p_dlbcl
dev.off()

pdf("Figure5C_TCGA_AML_ligand_correlations_volcanoplot.pdf", height = 4, width = 8)
p_aml_tcga
dev.off()

pdf("FigureS5F_TCGA_DLBCL_ligand_correlations_volcanoplot.pdf", height = 4, width = 7)
p_dlbcl_tcga
dev.off()

