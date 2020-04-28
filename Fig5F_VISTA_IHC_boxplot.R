
# Plot IHC percentages of VISTA+ cells in different hematological cancers (Figure 5F)

library(data.table)
library(ggplot2)
library(cowplot)
library(ComplexHeatmap)
library(dplyr)
library(ggpubr)

# read mIHC data
data <- fread("vista_ihc.txt", data.table = F)

# prepare data frame for plotting
df <- data %>%
  mutate(Disease_2 = gsub("BP CML", "CML", 
                          gsub("PH. B-ALL", "pre-B-ALL", 
                               gsub("HB", "Healthy BM", Disease))),
         VISTA = VISTA.fraction*100) %>%
  mutate(Disease_2 = ifelse(FAB %in% c("4", "5"), "AML M4/M5", 
                            ifelse(Disease_2 == "AML", "AML other", Disease_2))) %>% 
  mutate(Disease_2 = factor(Disease_2, levels = c("AML M4/M5", "AML other", "CML", "pre-B-ALL", "T-ALL", "Healthy BM")))
         
# boxplot
p <- ggplot(df, aes(x=Disease_2, y=VISTA, fill = Disease_2)) +
 geom_boxplot(outlier.shape = NA) +
 geom_jitter(width = 0.1, color = "grey20") +
 scale_size_continuous(range = c(0.1, 2)) +
 ylab("% VISTA+ cells") +
 xlab("") +
 guides(fill = FALSE) +
  theme_cowplot() +
 theme(axis.text.x = element_text(angle=45, hjust=1)) +
 labs(color = "") +
 stat_compare_means(aes(label = ..p.signif..),
   method = "wilcox.test",
   method.args = list(alternative = "two.sided"),
   ref.group = "Healthy BM",
   label.y = 95)

# print
pdf("Figure5F_VISTA_IHC_boxplot.pdf", width = 3, height = 4)
p
dev.off()
