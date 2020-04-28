
# Plot boxplot of numbers of expressed CGAs in CoMMpass myeloma subtypes (Figure 6F)

library(ggpubr)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(dplyr)
library(cowplot)


# load data
load("MM_COMPASS_FM.Rdata")

# only samples with gexp data
fm <- fm[, !is.na(fm["N:GEXP:KRAS",])]

# create df for plotting
df <- as.data.frame(t(fm)) 

# add subtypes
subtypes <- rep("CCND1", length(fm))
subtypes[fm["B:SAMP:cancermap_subtypes_WHSC1_FGFR3_Ig",]==1] <- "FGFR3"
subtypes[fm["B:SAMP:cancermap_subtypes_Hyperdiploid",]==1] <- "Hyperdiploid/gain(11q)"
subtypes[fm["B:SAMP:cancermap_subtypes_Hyperdiploid_amp1q",]==1] <- "Hyperdiploid/gain(1q)"
subtypes[fm["B:SAMP:cancermap_subtypes_MAF_Ig",]==1] <- "MAF"
subtypes[fm["B:SAMP:cancermap_subtypes_TRAF3_Aberrated",]==1] <- "TRAF3"
subtypes[fm["B:SAMP:cancermap_cluster_CGA_Prolif",]==1] <- "CGA/Proliferative"
 
df$subtype <- subtypes

# function for boxplot of different subtypes
subtype_bplot <- function(GEXP, ALT, TITLE, YLAB){
  df <- df[,c(GEXP, ALT)] 
  colnames(df) <- c("gexp", "alt") 
  
  ylabel <- gsub("N:GEXP:", "", gsub(":::::", "", GEXP))
  
  df <- df %>%
    mutate(gexp = as.numeric(as.character(gexp))) %>%
    mutate(alt = factor(alt, levels = sort(unique(alt))))
  
  ggplot(df, aes(x = alt, y = gexp, fill = alt)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, color = "grey20") +
    scale_fill_manual(values = brewer.pal(length(unique(df$alt)), "Set1")) +#,
    ylab(YLAB) +
    xlab("") +
    ggtitle(TITLE) +
    theme(legend.position="none",
          axis.title = element_text(size = 16),
          axis.text = element_text(size = 16),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(size = 18, hjust = 0.5, face = "plain")) +
    theme(plot.margin = unit(c(0.5,1,0,0), "cm"))
}

# plot number of CGAs in subtypes
p <- subtype_bplot("N:SAMP:nCGA", "subtype", "CGAs expressed", "CGAs expressed")

pdf("Figure6F_CoMMpass_CGA_subtype_boxplot.pdf", height = 6, width = 4)
p
dev.off()
