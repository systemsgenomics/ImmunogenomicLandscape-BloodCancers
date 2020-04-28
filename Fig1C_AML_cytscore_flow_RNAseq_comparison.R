
# Plot comparison of AML BM T/NK% based on flow cytometry vs cytolytic score (Figure 1C)

library(dplyr)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(gridExtra)

# read in lymphocyte fractions out of total BM cells (gated in FlowJo) 
data <- read.table("aml_bm_tnk_cytscore_comparison.txt", header = TRUE)

# add column with T and NK fraction sum
data$T_NK <- data$T+data$NK

# fractions to percentages
data$T_NK_per <- data$T_NK*100
data$T_per <- data$T*100
data$NK_per <- data$NK*100


# Scatter plot for Figure 1C
pdf("Figure1C_AML_cytolytic_score_flow_scatter.pdf", height = 2.5, width = 2.5)
ggscatter(data, x = "Cytolytic_score", y = "T_NK_per",
            size = 1.5,
            add = "reg.line",  # Add regression line
            add.params = list(color = "blue", fill = "lightgray", size = 1), # Customize reg. line
            conf.int = TRUE # Add confidence interval
  ) +
  stat_cor(method = "spearman", label.x = 2.5, label.y = 25, label.sep = "\n") +
  xlab("Cytolytic score (RNA-seq)") +
  ylab("% T/NK cells (flow cytometry)")
dev.off()


# T and NK cells separately (Figure S1)

t <- ggscatter(data, x = "Cytolytic_score", y = "T_per",
          size = 1.5,
          add = "reg.line",  # Add regression line
          add.params = list(color = "blue", fill = "lightgray", size = 1), # Customize reg. line
          conf.int = TRUE # Add confidence interval
) +
  stat_cor(method = "spearman", label.x = 1.5, label.y = 25, label.sep = "\n") +
  xlab("Cytolytic score\n(RNA-seq)") +
  ylab("% T cells in AML BM\n(flow cytometry)") +
  ggtitle("T cells") +
  theme(plot.title = element_text(hjust = 0.5))

nk <- ggscatter(data, x = "Cytolytic_score", y = "NK_per",
               size = 1.5,
               add = "reg.line",  # Add regression line
               add.params = list(color = "blue", fill = "lightgray", size = 1), # Customize reg. line
               conf.int = TRUE # Add confidence interval
) +
  stat_cor(method = "spearman", label.x = 1.5, label.y = 7, label.sep = "\n") +
  xlab("Cytolytic score\n(RNA-seq)") +
  ylab("% NK cells in AML BM\n(flow cytometry)") +
  ggtitle("NK cells") +
  theme(plot.title = element_text(hjust = 0.5))


# Print pdf for Figure S1E
pdf("FigureS1E_AML_cytolytic_score_flow_TNKseparately_scatter.pdf", height = 2.5, width = 4)
grid.arrange(t, nk, ncol = 2)
dev.off()


