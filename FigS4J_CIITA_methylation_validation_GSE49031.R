library(limma)
library(edgeR)
library(parallel)
library(gridExtra)

GIT_HOME="/research/users/ppolonen/git_home/"
source(file.path(GIT_HOME, "common_scripts/featurematrix/functions_generate_fm.R"))
source(file.path(GIT_HOME, "common_scripts/visualisation/plotting_functions.R"))

setwd("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/Published_data_figures")

anno = read.delim("cytogenetics_GSE49031.txt", stringsAsFactors=F, header=T)

# extract interesting probes:
meth_probeinfo = read.csv("GPL13534_HumanMethylation450_15017482_v.1.1.csv",skip=7, stringsAsFactors=F, header=T)

# find CIITA locus  chr16:10959400-10965899
find=meth_probeinfo$CHR=="16"&meth_probeinfo$MAPINFO>10959400&meth_probeinfo$MAPINFO<11020097
locus_interest=meth_probeinfo[find,]

write.table(locus_interest$IlmnID, "probenames.txt", sep="\t", quote=F, row.names=F, col.names=F, append=F)
system("grep -f probenames.txt GSE49031_processed.txt > GSE49031_ciita.txt")

# can start here
header=unlist(strsplit(readLines("GSE49031_processed.txt", n = 1),split = "\t"))
meth = read.delim("GSE49031_ciita.txt", stringsAsFactors=F, header=F, row.names=1)
colnames(meth)=header[2:length(header)]

meth=meth[,grep("Average Beta", colnames(meth))]
colnames(meth)=gsub(" Average Beta", "", colnames(meth))

meth_subset=meth[,na.omit(match(anno[,1], colnames(meth)))]
anno=anno[match(colnames(meth_subset), anno[,1]),]


A=meth_subset[,anno$X11q23.MLL==1&!anno$T.ALL==1]
B=meth_subset[,anno$t.12.21.==1&!anno$T.ALL==1]
C=meth_subset[,anno$t.1.19.==1&!anno$T.ALL==1]
D=meth_subset[,anno$dic.9.20.==1&!anno$T.ALL==1]
E=meth_subset[,anno$t.9.22.==1&!anno$T.ALL==1]
G=meth_subset[,anno$T.ALL==1]

result=cbind(rownames(A), locus_interest$UCSC_RefGene_Name, 
             locus_interest$UCSC_RefGene_Accession, 
             locus_interest$UCSC_RefGene_Group,
             locus_interest$UCSC_CpG_Islands_Name,
             paste0(locus_interest$CHR, ":", locus_interest$MAPINFO,"-", locus_interest$MAPINFO+nchar(locus_interest$AlleleA_ProbeSeq))
)

d=data.frame(t(meth_subset), anno$Subtype)
library(reshape)

plotd=melt(d)

plotd=plotd[plotd[,2]=="cg04945379",]

logicalVectors=lapply(unique(anno$Subtype), function(a)anno$Subtype%in%a)
names(logicalVectors)=unique(anno$Subtype)

genelist=rownames(A)

plot.boxplot(gene = "cg04945379", logicalVectors = logicalVectors,data = meth_subset,order = T,spread = T)

# Check gene expression to confirm its higher:
data_gexp=get(load("reRMA210915.RData"))

library(org.Hs.eg.db)
## Bimap interface:
x <- org.Hs.egSYMBOL
# Get the gene symbol that are mapped to an entrez gene identifiers
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])

rownames(data_gexp)=xx[match(rownames(data_gexp), names(xx))]

# colnames(data)=gsub("_ALL*.*","", colnames(data))
colnames(data_gexp)=paste0("ALL_", gsub("GSM*.*_|.CEL.gz","", colnames(data_gexp)))

data2=data_gexp[,na.omit(match(anno[,1], colnames(data_gexp)))]

anno2=anno[match(colnames(data2), anno[,1]),]

logicalVectors=lapply(unique(anno2$Subtype), function(a)anno2$Subtype%in%a)
names(logicalVectors)=unique(anno2$Subtype)

genelist=c("CIITA", grep("HLA-D", rownames(data2), value=T))

annotv=colnames(data2)[colnames(data2)%in%colnames(meth_subset)[meth_subset["cg04945379",]>0.6]]
methv=meth_subset["cg04945379",match(colnames(data2), colnames(meth_subset))]

# also HLAII
dat_a3=data2[rownames(data2)%in%c("HLA-DMA","HLA-DMB","HLA-DPA1","HLA-DPB1","HLA-DRA","HLA-DRB1"),]

dat3=2^dat_a3+0.01
gm3=log2(t(apply(dat3, 2, gm_mean)))
rownames(gm3)="HLAIIScore"

data2=rbind(data2, gm3)

genelist=c("HLAIIScore", "CIITA")
p.all=lapply(genelist, plot.boxplot, logicalVectors = logicalVectors, data = data2, order.bl = T,spread = T, sample.color.continuous = as.numeric(methv), outlier.size = 2)
p.all=append(p.all, lapply(genelist, plot.boxplot, logicalVectors = logicalVectors, data = data2, order.bl = T,spread = T,  outlier.size = 2))
ggsave("FigS4J_GSE49031_CIITA_methylation.pdf", do.call(marrangeGrob, append(list(grobs=p.all, nrow=4, ncol=2),list(top=NULL))), width = 260 , height = 297, units = "mm", dpi=150, limitsize = FALSE)
