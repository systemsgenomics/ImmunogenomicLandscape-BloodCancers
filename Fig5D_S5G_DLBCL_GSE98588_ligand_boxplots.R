
# Figure 5D and S5G oxplots of immunomodulatory ligand expression stratified by genetic alterations in DLBCL (GSE98588 Chapuy et al.)

library(ggpubr)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(dplyr)
library(cowplot)
library(readxl)


# load data
load("GSE98588_fm.Rdata")

# create matrix with selected genes
genes = as.character(read_excel("ligands.xlsx")$Gene)

feat <- paste0("N:GEXP:", paste0(genes, ":::::"))

rownames(fm)[grepl("6P", rownames(fm))]

df <- as.data.frame(t(fm)) 


# function for boxplot of altered vs unaltered cases
alt_bplot <- function(GEXP, ALT, TITLE){
  df <- df[,c(GEXP, ALT)] 
  colnames(df) <- c("gexp", "alt") 
  
  ylabel <- gsub("N:GEXP:", "", gsub(":::::", "", GEXP))
  
  alt_short <- gsub("([0-9])P([0-9])", "\\1p\\2", gsub("([0-9])Q([0-9])", "\\1q\\2", gsub("Q$", "q", gsub("P$", "p", gsub("B:.....", "", gsub("_", ".", gsub("SV_", "", gsub("_nonsynonymous:::::", "", gsub(":AMP", " amp ", gsub(":LOSS", " loss", gsub(":DEL", " del ", gsub(":SV", " SV  ", ALT))))))))))))
  alt_short <- ifelse(grepl("GNAB", ALT), paste0(alt_short, " mut "), alt_short)
  
  df <- df %>%
    mutate(gexp = as.numeric(as.character(gexp)),
           alt = gsub("^1", alt_short, 
                      gsub("^0", paste(gsub(".{5}$", "", alt_short), "WT "), alt))) %>%
    mutate(alt = factor(alt, levels = sort(unique(alt))))
  
  ggplot(df, aes(x = alt, y = gexp, fill = alt)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.1, color = "grey20") +
    ylim(c(min(df$gexp)-(max(df$gexp)-min(df$gexp))*0.1, max(df$gexp)+(max(df$gexp)-min(df$gexp))*0.1)) +
    scale_fill_manual(values = brewer.pal(3, "Paired")[c(2,1)]) +
    ylab(substitute(italic(x)~" expression (log2)", list(x = ylabel))) +
    xlab("") +
    ggtitle(TITLE) +
    theme(legend.position="none",
          axis.title = element_text(size = 16),
          axis.text = element_text(size = 16),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(size = 18, hjust = 0.5)) +
    stat_compare_means(aes(label = ..p.signif..),
                       method = "wilcox.test",
                       method.args = list(alternative = "two.sided"),
                       comparisons = list(c(alt_short, paste(gsub(".{5}$", "", alt_short), "WT "))),
                       label.y = max(df$gexp)+(max(df$gexp)-min(df$gexp))*0.05
    ) +
    theme(plot.margin = unit(c(0.5,1,0,0), "cm"))
}

# plot selected alterations
p1 <- alt_bplot("N:GEXP:CD70", "B:GNAB:CD70", expression(paste(italic("CD70"))))
p2 <- alt_bplot("N:GEXP:MICB", "B:CNVR:6P21_33:LOSS", expression(paste(italic("MICB"))))

p6 <- alt_bplot("N:GEXP:TNFRSF14", "B:CNVR:BCL2:SV", expression(paste(italic("TNFRSF14"))))
p7 <- alt_bplot("N:GEXP:CD274", "B:CNVR:9P24_1:AMP", expression(paste(italic("CD274\n(PD-L1)"))))
p8 <- alt_bplot("N:GEXP:PDCD1LG2", "B:CNVR:9P24_1:AMP", expression(paste(italic("PDCD1LG2\n(PD-L2)"))))

pdf("Figure5D_DLBCL_GSE98588_CD70_boxplot.pdf", height = 5, width = 2)
p1 + scale_y_continuous(breaks = c(3,6,9,12), limits = c(2,12.5))
dev.off()

pdf("FigureS5G_DLBCL_GSE98588_MICB_TNFRSF14_CD274_PDCD1LG2_boxplot.pdf", height = 4, width = 10)
grid.arrange(p2, p6, p7, p8, ncol = 4)
dev.off()
