
# Find DE genes between TGCA AML and DLBCL and plot volcano plot with differential methylation (Figure S5D)

library(circlize)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(circlize)
library(gridExtra)
library(dplyr)
library(limma)
library(readxl)
library(edgeR)

# get gene list
genelist <- as.character(read_excel("ligands.xlsx")$Gene)


# load data

# gene expression
m1_gexp=get(load("/Users/odufva/Documents/Labra/Projects/Hemap_immunology/R analysis/DATA_SHARING/DLBC.rnaseqv2.counts.Rdata"))
m2_gexp=get(load("/Users/odufva/Documents/Labra/Projects/Hemap_immunology/R analysis/DATA_SHARING/LAML.rnaseqv2.counts.Rdata"))

gexp=data.matrix(cbind(m1_gexp, m2_gexp))

# methylation
m1_meth=get(load("/Users/odufva/Documents/Labra/Projects/Hemap_immunology/R analysis/DATA_SHARING/TCGA_AML_meth_probes_genelist.Rdata"))
m2_meth=get(load("/Users/odufva/Documents/Labra/Projects/Hemap_immunology/R analysis/DATA_SHARING/TCGA_DLBCL_meth_probes_genelist.Rdata"))
meth_annot <- meth_annot_extra[meth_annot_extra$nearestTSS %in% genelist & meth_annot_extra$distanceToTSS<1000,]

# merge DLBCL and AML and take mean of probes within 1000f bp of TSS
meth <- merge(m1_meth, m2_meth, by.x = "row.names", by.y = "row.names")
meth <- merge(meth, meth_annot, by.x = "Row.names", by.y = "methProbeIDs")
meth <- aggregate(. ~ meth$nearestTSS, meth[,grepl("TCGA", colnames(meth))], mean)
rownames(meth) <- meth[,1]
meth <- meth[,-1]

# select matching cases
gexp <- gexp[,colnames(gexp) %in% colnames(meth)]
meth <- meth[,colnames(meth) %in% colnames(gexp)]
meth <- meth[,colnames(gexp)]


## DE analysis

# at least 10 samples with cpm>1
keep.exprs <- rowSums(cpm(gexp)>1)>=10
gexp <- gexp[keep.exprs,]

# voom transform
v=voom(gexp+0.01, normalize.method = "quantile")$E

class = c(rep("DLBCL", dim(gexp)[2]))
class[colnames(gexp) %in% colnames(m2_gexp)] = "AML"

design <- model.matrix(~0+class)
colnames(design) <- c("DLBCL","AML")

cont=makeContrasts(DLBCL - AML, levels = design)
fit <- lmFit(v, design)
fit <- contrasts.fit(fit, cont)

fit2 <- eBayes(fit)

res=topTable(fit2, p.value = 0.05, number = 999999)
res_all=topTable(fit2, p.value = 1, number = 999999)

pl=res[res$ID%in%genelist,1]
pl=pl[!pl%in%c("TNFRSF14")]#, "ENTPD1", "ARG1")] # temporary



## Differential methylation analysis
meth_mval <- log2(meth/(1-meth))

class = c(rep("DLBCL", dim(meth)[2]))
class[colnames(meth) %in% colnames(m2_meth)] = "AML"

design <- model.matrix(~0+class)
colnames(design) <- c("DLBCL","AML")

cont=makeContrasts(DLBCL - AML, levels = design)
fit <- lmFit(meth_mval, design)
fit <- contrasts.fit(fit, cont)

fit2 <- eBayes(fit)

res_meth=topTable(fit2, p.value = 1e-12, number = 999999)
res_meth_all=topTable(fit2, p.value = 1, number = 999999)


## prepare data for volcano plot

res_all_df <- res_all
res_all_df$log10.adj.P.Val <- -log10(res_all_df$adj.P.Val)

res_meth_all_df <- res_meth_all
res_meth_all_df$log10.adj.P.Val <- -log10(res_meth_all_df$adj.P.Val)
res_meth_all_df$ID <- rownames(res_meth_all_df)


## combine gexp and meth differential analysis results
res_gexp_meth_df <- merge(res_all_df, res_meth_all_df, by = "ID", suffixes = c("_gexp", "_meth"))

# print supplement table
res_gexp_meth_df_clean <- res_gexp_meth_df %>% arrange(adj.P.Val_gexp) %>% select(-c(log10.adj.P.Val_gexp, log10.adj.P.Val_meth))

colnames(res_gexp_meth_df_clean) <- gsub("AveExpr", "Average", gsub("_gexp", ".expression", gsub("_meth", ".methylation", gsub("ID", "Gene", colnames(res_gexp_meth_df_clean)))))

write.csv(res_gexp_meth_df_clean, "TCGA_AML_DLBCL_ligand_differential_expression_methylation_supplement.csv", quote = F, row.names = F)

# plot volcano plot
pdf("FigureS5D_TCGA_DLBCL_AML_ligand_expression_methylation_volcanoplot.pdf", width = 7, height = 4.5) 
ggplot(res_gexp_meth_df, aes(x = logFC_gexp, y = log10.adj.P.Val_gexp, size = log10.adj.P.Val_meth)) +
  geom_point(color = "#d4d3d1") +
  geom_point(color = ifelse(abs(res_gexp_meth_df$logFC_meth)>1.5 & sign(res_gexp_meth_df$logFC_gexp)==sign(res_gexp_meth_df$logFC_meth), "#fea719", NA)) +#, "#a52a2a")))") +
  geom_label_repel(aes(label=ifelse((abs(logFC_meth)>1.5 & sign(res_gexp_meth_df$logFC_gexp)==sign(res_gexp_meth_df$logFC_meth)) |
                                      (log10.adj.P.Val_gexp>12.5 & logFC_gexp>1.5) |
                                      log10.adj.P.Val_gexp>40, ID, '')),
                   point.padding = 0.5,
                   label.padding = 0.15,
                   segment.size = 0.25,
                   #box.padding = 1,
                   size = 3.5,
                   color = "black") +
  ylab("-log10(adjusted P value)") +
  xlab("Expression fold change (log2)") +
  guides(size = guide_legend(title = "Differential\nmethylation\nadjusted\nP value (-log10)")) +
  scale_color_distiller(palette = "YlOrRd", direction = 2) +
  xlim(c(-10, 10))
dev.off()

## ----------------------------------------------------------

## Differential expression and methylation boxplots (Figure S5E)

# transpose matrices
v_df <- as.data.frame(t(v[rownames(v) %in% genelist,]))
meth_df <- as.data.frame(t(meth))

# merge data frames
df <- merge(v_df, meth_df, by = "row.names", suffixes = c("_gexp", "_meth"))
rownames(df) <- df$Row.names
df$Row.names <- NULL

# add disease annotations
df$cancer <- rep("DLBCL", dim(df)[1])
df$cancer[rownames(df) %in% colnames(m2_gexp)] = "AML"
df$cancer <- factor(df$cancer, levels = c("DLBCL", "AML"))

# boxplot function (continuous methylation color)
boxplot_gexp_meth_continuous <- function(gene){
  
  ggplot(df, aes_string(x = "cancer", y = paste0(gene, "_gexp"), color = paste0(gene, "_meth"), fill = "cancer")) +
    geom_boxplot(outlier.shape = NA, show.legend = FALSE) +
    geom_jitter(width = 0.1) +
    scale_fill_manual(values = c("mediumseagreen", "#c87137")) +
    scale_color_gradientn(colors = c("#a52a2a", "#d4d3d1", "#fea719"), breaks = c(0, 0.5, 1), limits = c(0, 1), guide = "colourbar") +
    ylab("Expression (log2)") +
    xlab("") +
    ggtitle(gene) +
    theme_cowplot() +
    theme(axis.title = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          axis.text.x = element_text(size = 16, angle = 45, hjust = 1),
          plot.title = element_text(size = 18, hjust = 0, face = "italic"),
          plot.margin = unit(c(-0.5,0,0,0.5),"cm")) +
    guides(color = FALSE,
          fill = FALSE)
}

# print Figure S5E plots
p1 <- lapply(as.character(c("PDCD1LG2", "CD80", "CLEC2D", "PVR")), boxplot_gexp_meth_continuous)
m1 <- marrangeGrob(p1, nrow = 2, ncol = 2)
ggsave("FigureS5E_TCGA_ligand_expression_methylation_boxplots.pdf", m1, height = 6, width = 4)

