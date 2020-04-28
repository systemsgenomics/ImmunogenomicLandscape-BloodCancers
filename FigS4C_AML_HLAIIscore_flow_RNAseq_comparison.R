
# Plot comparison of AML BM % of HLA-DR+ cells based on flow cytometry vs HLA II score (Figure S4C)

library(dplyr)
library(ggplot2)
library(ggpubr)
library(cowplot)

# read in data
data <- read.table("aml_bm_hla2tscore_comparison.txt", header = TRUE)

# HLA II score vs. flow HLA-DR+ blasts scatter plot

data$HLADRpos_out_of_blasts_per <- data$HLADRpos_out_of_blasts*100

pdf("FigureS4C_AML_HLAIIscore_flow_RNAseq_scatter.pdf", height = 3, width = 3)
ggscatter(data, x = "HLAII_score", y = "HLADRpos_out_of_blasts_per",
          size = 1.5,
          add = "reg.line",  # Add regression line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE # Add confidence interval
) +
  stat_cor(method = "spearman", label.x = 3, label.y = 110, label.sep = "\n") +
  xlab("HLA II score (RNA-seq)") +
  ylab("% HLA-DR+ blasts (flow cytometry)") +
  scale_y_continuous(limits = c(-30,125), breaks = c(0, 25, 50, 75, 100)) +
  xlim(c(3, 10.5))
dev.off()
