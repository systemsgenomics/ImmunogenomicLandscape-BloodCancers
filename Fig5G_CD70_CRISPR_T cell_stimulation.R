
# Analyze JE6-TPR stimulation experiments with SUDHL5 CD70 knockout co-culture (Figure 5G)

library(data.table)
library(readxl)
library(dplyr)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(RColorBrewer)
library(ggpubr)


# load data
data <- fread("cd70_crispr_flow_cytometry.txt", data.table = F)

# calculate percentages of NFAT/NFkB-activated T cells
# LL (lower left) = double negative
# LR (lower right) = NFAT only
# UL (upper left) = NFkB only
# UR (upper right) = double positive

df <- data %>%
  mutate(condition = paste(TPR, TCS, SUDHL5, sep = "_")) %>%
  mutate(NFAT_NFkB_LL_percentage = 100*Count_of_NFAT_NFkB_LL/Count_of_TPR,
         NFAT_NFkB_LR_percentage = 100*Count_of_NFAT_NFkB_LR/Count_of_TPR,
         NFAT_NFkB_UL_percentage = 100*Count_of_NFAT_NFkB_UL/Count_of_TPR,
         NFAT_NFkB_UR_percentage = 100*Count_of_NFAT_NFkB_UR/Count_of_TPR)


# Barplot of % NFkB only-activated Jurkat-CD27 cells stimulated with TCS-Ctrl (stimulator cells with anti-CD3) + SUDHL5 with sgCD70/Ctrl sgRNAs

PARAMETER = "NFAT_NFkB_UL_percentage"

df_values <- df[df$TPR == "TPR-CD27" & (df$TCS=="TCS-Ctrl"|(df$TCS=="no_TCS"&df$SUDHL5=="no_SUDHL5")), c("TPR", "TCS", "SUDHL5", PARAMETER)]
colnames(df_values) <- c("TPR", "TCS", "SUDHL5", "parameter")
  
  df_values <- df_values %>%
    mutate(condition = ifelse(TCS=="no_TCS"&SUDHL5=="no_SUDHL5", "no stimulation", gsub(" \\+ no SUDHL5", "", gsub("_|\\.", " ", gsub("Ctrl", "Control", paste("anti-CD3 +", SUDHL5))))))
  
  df_stats <- df_values %>%
    group_by(condition) %>% 
    summarize(mean = mean(parameter),
              sd = sd(parameter))
  
  df_plot <- merge(df_values, df_stats) %>%
    mutate(condition = factor(condition, levels = c("no stimulation", "anti-CD3", "anti-CD3 + SUDHL5 sgControl", "anti-CD3 + SUDHL5 sgCD70 1", "anti-CD3 + SUDHL5 sgCD70 2", "anti-CD3 + SUDHL5 sgCD70 3")))
    
  
  p <- ggplot(df_plot, aes(x = condition, y = parameter)) +
    geom_bar(aes(x = condition, y = mean, fill = condition), stat = "identity", position = position_dodge(), color = "grey20") +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),
                  width=0.3,
                  size=0.3,
                  color = "grey20") +
    geom_jitter(aes(x = condition, y = parameter), width = 0.1, color = "grey20", show.legend = F) +
    scale_fill_manual(values = brewer.pal(11, "RdGy")[c(7,8,9,4,3,2)]) +
    ylab("% NFkB-activated\nT cells") +
    xlab("") +
    scale_y_continuous(expand = c(0,0), limits = c(0,31)) +
    guides(fill=FALSE) +
    stat_compare_means(data = df_plot, aes(x = condition, y = PARAMETER, label = ..p.signif..),
                       inherit.aes = F,
                       method = "t.test",
                       comparisons = list(c("anti-CD3", "no stimulation"),
                                          c("anti-CD3 + SUDHL5 sgControl", "no stimulation"),
                                          c("anti-CD3 + SUDHL5 sgCD70 1", "anti-CD3 + SUDHL5 sgControl"),
                                          c("anti-CD3 + SUDHL5 sgCD70 2", "anti-CD3 + SUDHL5 sgControl"),
                                          c("anti-CD3 + SUDHL5 sgCD70 3", "anti-CD3 + SUDHL5 sgControl")),
                       label.y = c(21, 23, 25, 27, 29)) +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle=45, hjust=1),
          legend.title=element_blank())
  
  pdf("Figure5G_CD70_CRISPR_T_cell_stimulation.pdf", width = 2.75, height = 4.5)
  p
  dev.off()
  
  