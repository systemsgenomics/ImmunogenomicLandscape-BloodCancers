
# Plots boxplots of cytolytic score startified by genetic alterations in DLBCL (GSE98588 Chapuy et al.) for Figure 3B-C 

library(ggpubr)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(dplyr)
library(cowplot)


# load data
load("GSE98588_fm.Rdata")

# transpose data frame
df <- as.data.frame(t(fm)) 

# function for boxplot of altered vs unaltered cases by subtype
alt_bplot_subtype <- function(GEXP, ALT, ALTNAME, TITLE){
  df <- df[,c(GEXP, ALT, "B:SAMP:COO_byGEP_ABC", "B:SAMP:COO_byGEP_GCB", "B:SAMP:COO_byGEP_Unclassified")] 
  colnames(df) <- c("gexp", "alt", "abc", "gcb", "unclassified") 
  
  ylabel <- gsub("N:GEXP:", "", gsub(":::::", "", GEXP))
  
  df <- df %>%
    filter(abc==1|gcb==1) %>% # remove unclassified
    mutate(gexp = as.numeric(as.character(gexp)),
           alt = gsub("1", ALTNAME, 
                      gsub("0", "WT",alt))) %>%
    mutate(alt = factor(alt, levels = c(ALTNAME, "WT"))) %>%
    mutate(alt_subtype = ifelse(abc==1, paste(alt, "ABC", sep = " "), 
                                ifelse(gcb==1, paste(alt, "GCB", sep = " "), ""))
    ) %>%
    mutate(alt_subtype = factor(alt_subtype, levels = c(paste(ALTNAME, "ABC", sep = " "), "WT ABC", paste(ALTNAME, "GCB", sep = " "), "WT GCB")))
  
  
  ggplot(df, aes(x = alt_subtype, y = gexp, fill = alt_subtype)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.1, color = "grey20") +
    ylim(c(min(df$gexp)-(max(df$gexp)-min(df$gexp))*0.1, max(df$gexp)+(max(df$gexp)-min(df$gexp))*0.1)) +
    scale_fill_manual(values = brewer.pal(6, "Paired")[c(6,5,4,3)]) +
    ylab("Cytolytic score (log2)") +
    xlab("") +
    ggtitle(TITLE) +
    theme_cowplot() +
    theme(legend.position="none",
          axis.title = element_text(size = 16),
          axis.text = element_text(size = 16),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0)) +
    stat_compare_means(aes(label = ..p.signif..),
                       method = "wilcox",
                       method.args = list(alternative = "two.sided"),
                       comparisons = list(c(1, 2), c(3, 4)),
                       label.y = max(df$gexp)+(max(df$gexp)-min(df$gexp))*0.075,
                       label.x.npc = "center"
    ) +
    theme(plot.margin = unit(c(0.5,1,0,1), "cm"))
}

# plot selected alterations
p1 <- alt_bplot_subtype("N:SAMP:CytolyticScore", "B:CNVR:7Q:AMP_GAIN", "7q amp/gain", expression(paste("7q amplification/gain")))
p2 <- alt_bplot_subtype("N:SAMP:CytolyticScore", "B:GNAB:DTX1", "DTX1 mut", expression(paste(italic("DTX1"), " mutation")))

# print pdf
pdf("Figure3BC_DLBCL_cytscore_boxplots.pdf", height = 4.5, width = 6)
grid.arrange(p1, p2, ncol = 2)
dev.off()

