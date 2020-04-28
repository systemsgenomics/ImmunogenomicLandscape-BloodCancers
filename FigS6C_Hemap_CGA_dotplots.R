GIT_HOME="/research/users/ppolonen/git_home/ImmunogenomicLandscape-BloodCancers/"
source(file.path(GIT_HOME, "common_scripts/visualisation/plotting_functions.R"))

# Plot CGA expression dot plots for (Figure S6C)

library(parallel)
library(gridExtra)
library(ggplot2)
library(ggrepel)
library(cowplot)

# load Hemap gene expression data
data = get(load("data9544_with_gene_symbols.RData"))

# load Hemap annotations
annot = get(load("Hemap_immunology_Annotations.Rdata"))

data=data[annot[,1],]

# plotting functions
FUN_PLOT=function(gene, logicalVectors, namesLV, data=NULL, matrix=NULL, col=NULL, ORDER=F, RANGE=NULL) {
  if(is.null(matrix)&is.null(data))stop("No data to plot, check data/matrix")
  
  GNAME=gsub("N:....:|:::::|DUFVA_", "", gene)
  GNAME=gsub("_", " ", GNAME)
  namesLV=gsub("Cancer_", " ", namesLV)
  
  cols <- read.table("colors_hemap_immunology.tsv", header = TRUE, sep = "\t", comment.char = " ")
  samples <- gsub("LCH", "MDS", gsub("AITL|PTCLNOS|ALCL", "TCL", gsub("FL|MALT", "BCL", gsub("Healthy", "NonCancerHealthy", gsub("_", "", names(logicalVectors))))))

  if(is.null(col)){
    col=as.character(cols[match(samples, cols$sample),2])
  }
  
  
  if(!is.null(matrix)){
    gene2=ifelse(grepl("GEXP", gene), gene, paste0("'N:GEXP:", gene, ":::::'"))
    D=as.numeric(read.delim(pipe(paste0("grep -Fw ", gene2, " ", matrix)), row.names = 1, header=F))
  }
  
  if(!is.null(data)){
    D = as.numeric(data[,colnames(data)%in%gene])
  }
  
  bplot_list=lapply(logicalVectors, function(v){
    D[v]
  })
  names(bplot_list)=gsub("_", " ", namesLV)
  
  if(ORDER){
    ord=sapply(bplot_list, median)
    col=col[order(ord, decreasing = T)]
    bplot_list=bplot_list[order(ord, decreasing = T)]
    
  }
  
  plots=FUNCTION_PLOT_LIST(bplot_list, gene, col, ORDER, RANGE)
  return(plots)
}


FUNCTION_PLOT_LIST=function(bplot_list, GNAME, col, ORDER, RANGE){
  
  df=melt(bplot_list)
  
  df$class <- factor(df[,2], levels = unique(as.character(df[,2])),ordered = TRUE)
  
  df$Expression=as.numeric(as.vector(df[,1]))
  p <- ggplot(data=df, aes(x=class, y=Expression, color=class)) +  
    geom_jitter(width = 0.25, size = 0.1) +
    scale_color_manual(values = col)
  
  p2 <- p +
    
    #theme with white background
    theme_bw() +
    
    # titles
    theme(plot.title = element_text(face="italic", color="black", size=16, hjust=0)) +
    theme(axis.title = element_text(color="black", face=NULL, size=12,angle = 90)) +
    theme(axis.title.y = element_text(size = 14, angle = 90, color="black", face=NULL)) +
    guides(color = FALSE) +
    
    ylab("Expression (log2)") +
    xlab("") +
    labs(title=GNAME) +
    #eliminates background, gridlines, and chart border
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    theme(panel.border= element_blank())+

    #draws x and y axis line
    theme(axis.line = element_line(),
          axis.line.x = element_line(color="black", size = 0.5),
          axis.line.y = element_line(color="black", size = 0.5)) +
    
    # X - axis text
    theme(axis.text.x = element_text(angle=45, hjust=1, color="black", size = 14, face=NULL),
          axis.text.y = element_text(hjust=1, color="black", size = 12, face=NULL))+ 
    
    # if want to limit to range
    if(!is.null(RANGE))scale_y_continuous(breaks=seq(2,14,2), limits = RANGE)
  
  return(p2)
}


boxplots_grid_topcga_subtypes <- function(x){
  p.all=lapply(c("MAGEC1", "MAGEC2", "MORC1", "DSCR8", "MAGEB1", "MAGEB2", "ADAM29", "DMRT1", "SAGE1"), FUN_PLOT, logicalVectors, namesLV=names(logicalVectors), data=data, ORDER=F)
  ggsave(paste0(GENELIST, ".pdf"), do.call(marrangeGrob, append(list(grobs=p.all, nrow=3, ncol=3),list(top=NULL))), width = 350 , height = 250, units = "mm", dpi=250)
}


# make logical vector with cancer subtypes and healthy sample
annot$logicalvector <- annot$Sample.type
annot$logicalvector[annot$Sample.type!="NonCancerHealthy"] <- annot$Category.specifying.subtype[annot$Sample.type!="NonCancerHealthy"]
annot$logicalvector[annot$Sample.type=="NonCancerHealthy"] <- "Healthy"
annot$logicalvector[annot$Category.specifying.lineage.tumor.origin=="AML"] <- "AML"
annot$logicalvector[annot$Category.specifying.lineage.tumor.origin=="MDS"] <- "MDS"
annot$logicalvector[annot$Category.specifying.lineage.tumor.origin=="CLL"] <- "CLL"
annot$logicalvector[annot$Category.specifying.subtype=="AILT"] <- "AITL"
annot$logicalvector[annot$Category.specifying.subtype=="LC"] <- "LCH"

# plot selected genes 
logicalVectors=get.logical(annovector = list(annot$logicalvector), filterv = annot$Sample.type%in%c("Cancer", "Prolif", "NonCancerHealthy"), PREFIX = "")
logicalVectors=logicalVectors[paste0(c("MM", "DLBCL", "MCL", "CHL", "FL", "MALT", "ALCL", "PTCLNOS", "AITL", "T-ALL", "pre-B-ALL", "AML", "MDS", "LCH", "CML", "CLL", "Healthy"), "_")]
GENELIST="FigureS6C_Hemap_CGA_dotplots"
boxplots_grid_topcga_subtypes()
