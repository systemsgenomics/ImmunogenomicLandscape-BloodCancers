source('/home/users/ppolonen/scripts/Projects/Dufva_helsinki_2016/CIBERSORT.R')

# Options:
# i)  perm = No. permutations; set to >=100 to calculate p-values (default = 0)
# ii) QN = Quantile normalization of input mixture (default = TRUE)
# Input: signature matrix and mixture file, formatted as specified at http://cibersort.stanford.edu/tutorial.php
# Output: matrix object containing all results and tabular data written to disk 'CIBERSORT-Results.txt'
# License: http://cibersort.stanford.edu/CIBERSORT_License.txt

setwd("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/")

# read annot for filtering
annot = get(load("Hemap_immunology_Annotations.Rdata"))

#read in data, filter
X <- read.table('LM22.txt',header=T,sep="\t",row.names=1,check.names=F)
Y=get(load("/research/groups/sysgen/PROJECTS/HEMAP/dat2figs/data/data9544/data9544_FINAL.RData"))
Y=Y[,colnames(Y)%in%annot$GSM.identifier..sample.]

# results will go to WD
results <- CIBERSORT(X,Y, perm=100, QN=F)

# Also MCP counter
library(MCPcounter)
ExampleEstimates=MCPcounter.estimate(Y,featuresType="HUGO_symbols")
rownames(ExampleEstimates)=paste0("MCP_", gsub(" ", "_", rownames(ExampleEstimates)))

save(ExampleEstimates, file="MCP_counter_data.Rdata")


#****************************** TCGA AML data *******************************************
setwd("/home/work/public/biowhat/hemap_submission/dat2figs/feat_matrices/TCGA_AML")

#read in data
X <- read.table('/home/users/ppolonen/scripts/Projects/Dufva_helsinki_2016/LM22.txt',header=T,sep="\t",row.names=1,check.names=F)
Y <- read.table('/home/work/public/biowhat/petri/ISB_data/TCGA/laml/laml.unc.edu__illuminahiseq_rnaseqv2__rnaseqv2.test.txt',header=T,sep="\t",row.names=1,check.names=F)

# count is high enough in 25% of the samples
filt=rowSums(Y>1) >= dim(Y)[2]/4

Y=Y[filt,]

# variance filter, low variance
vars <- apply(Y, 1, var, na.rm = T)

Y=Y[vars > quantile(vars, 0.15, na.rm = T), ]


# make data continuous
library(limma)
v<-voom(Y,plot=TRUE)

Y=v$E

results <- CIBERSORT(X,Y, perm=100, QN=F)

colnames(results)=paste0("N:SAMP:CIBERSORT_", gsub(" |-", "_", colnames(results)), ":::::")
results_fm=t(results)
rownames(results_fm)

library(MCPcounter)
ExampleEstimates=MCPcounter.estimate(Y,featuresType="HUGO_symbols")
rownames(ExampleEstimates)=paste0("MCP_", gsub(" ", "_", rownames(ExampleEstimates)))
save(ExampleEstimates, file="/research/groups/sysgen/PROJECTS/HEMAP/dat2figs/feat_matrices/TCGA_AML/MCP_counter_results.Rdata")



load("/home/work/public/biowhat/hemap_submission/dat2figs/feat_matrices/TCGA_AML/TCGA_AML_fm_cancermap_clusters_070915.Rdata")
matrix=rbind(matrix, results_fm)
save(matrix, file="/home/work/public/biowhat/hemap_submission/dat2figs/feat_matrices/")


#****************************** TCGA DLBCL data *******************************************
setwd("/research/groups/sysgen/PROJECTS/HEMAP/dat2figs/feat_matrices/TCGA_DLBCL/")

#read in data
X <- read.table('/home/users/ppolonen/scripts/Projects/Dufva_helsinki_2016/LM22.txt',header=T,sep="\t",row.names=1,check.names=F)
Y_org <- read.table('/home/groups/biowhat/public_data/TCGA/firehose/RNAseq_counts/gdac.broadinstitute.org_DLBC.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2016012800.0.0/DLBC.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt',header=T,sep="\t",row.names=1,check.names=F, stringsAsFactors = F)

# process firehose
Y=data.matrix(Y_org[-1,Y_org[1,]%in%"raw_count"])
colnames(Y)=substr(colnames(Y), 1, 15)
rownames(Y)=gsub("\\|.*", "", rownames(Y))
Y=Y[!rownames(Y)%in%"?",]

# filter data
library(edgeR)
keep <- rowSums(cpm(Y) > 1) >= round(dim(Y)[2]*0.15)
Y=Y[keep,]

# make data continuous
library(limma)
Y<-voom(Y,plot=TRUE)$E

# variance filter, low variance
# vars <- apply(Y, 1, var, na.rm = T)
# Y=Y[vars > quantile(vars, 0.15, na.rm = T), ]
# 
results <- CIBERSORT(X,Y, perm=100, QN=F)
colnames(results)=paste0("N:SAMP:DUFVA_CIBERSORT_", gsub(" |-", "_", colnames(results)), ":::::")
results_fm=t(results)
rownames(results_fm)

save(results_fm, file="/research/groups/sysgen/PROJECTS/HEMAP/dat2figs/feat_matrices/TCGA_DLBCL/CIBERSORT_results.Rdata")


# MCP counter
library(MCPcounter)
ExampleEstimates=MCPcounter.estimate(Y,featuresType="HUGO_symbols")
rownames(ExampleEstimates)=paste0("N:SAMP:DUFVA_MCP_", gsub(" ", "_", rownames(ExampleEstimates)), ":::::")

save(ExampleEstimates, file="/research/groups/sysgen/PROJECTS/HEMAP/dat2figs/feat_matrices/TCGA_DLBCL/MCP_counter_results.Rdata")