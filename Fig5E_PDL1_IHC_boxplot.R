
# Plot IHC percentages of PD-L1+ cells in different hematological cancers (Figure 5E)

library(data.table)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(dplyr)
library(ggpubr)

# read mIHC data
data <- fread("pdl1_ihc.txt", data.table = F)

# order cancers for plotting

data <- data %>% 
  mutate(Cancer_type = factor(Cancer_type, levels = c("DLBCL (non-GCB)", "DLBCL (GCB)", "pre-B-ALL", "T-ALL", "AML", "CML", "Healthy BM")))
  
# boxplot of all cancers with ABC and GCB separately
p <- ggplot(data, aes(x = Cancer_type, y = PDL1, fill = Cancer_type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, color = "grey20") +
  #geom_jitter(data = df[df$Subtype=="Blast phase CML",], aes(color = Subtype), width = 0.1) +
  scale_size_continuous(range = c(0.1, 2)) +
  scale_y_continuous(limits = c(-0.1,112), breaks = c(0, 25, 50, 75, 100)) +
  ylab("% PD-L1+ cells") +
  xlab("") +
  guides(fill = FALSE) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("DLBCL (GCB)" = "rosybrown", "DLBCL (non-GCB)" = "rosybrown", "AML" = "mediumseagreen", "CML" = "darkgreen", "pre-B-ALL" = "violet", "T-ALL" = "cornflowerblue", "Healthy BM" = "grey50")) +
  labs(color = "") +
  stat_compare_means(ref.group = "Healthy BM",
                     label = "p.signif",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")),
                     method = "wilcox.test",
                     label.y = 112) +
  stat_compare_means(comparisons = list(c("DLBCL (GCB)", "DLBCL (non-GCB)")),
                     label = "p.signif",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")),
                     method = "wilcox.test",
                     label.y = 105)



# print
pdf("Figure5E_PDL1_IHC_boxplot.pdf", width = 3, height = 4.25)
p
dev.off()
