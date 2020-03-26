GIT_HOME="/research/users/ppolonen/git_home/ImmunogenomicLandscape-BloodCancers/"
source(file.path(GIT_HOME, "common_scripts/visualisation/plotting_functions.R"))

library(parallel)
library(org.Hs.eg.db)
library(seqinr)
library(reshape2)

setwd("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/Published_data_figures")

# load annotations
annot = get(load("Hemap_immunology_Annotations.Rdata"))

# Genes that are expressed in normal cells:
normal_res=get(load("dufva_mw_hgt_pval_normal_expressed.Rdata"))
OE_cancer=get(load("dufva_mw_hgt_pval_cancer_overexpressed.Rdata"))

# T/NK failed genes
tnk_hgt=get(load("dufva_mw_hgt_pval_all_TNKspecific_genes.Rdata"))
tnk_genes=tnk_hgt[tnk_hgt$name%in%"TNK"&tnk_hgt$adj.pvalue>3&tnk_hgt$FoldChange>1.5,1]

normal_res=normal_res[!normal_res[,1]%in%tnk_genes,]
OE_cancer=OE_cancer[!OE_cancer[,1]%in%tnk_genes,]

# profiles
profile=get(load("mixtureM_profile.Rdata"))
profile[profile==-1] = 0
profile=profile[,colnames(profile)%in%annot$GSM.identifier..sample.]

# filter using hemap genes:
# annotations
load("data9544_with_gene_symbols.RData")
data = t(data[rownames(data)%in%annot$GSM.identifier..sample.,])

# take only high expressed into account
profile[data<5]=0

#*******************************************************************************
#*******************************************************************************
#************************** genes expressed/overexpressed in normal cells**************************
genes_normals=unique(normal_res[normal_res$adj.pvalue>2&normal_res$FoldChange>0&!grepl("_CD34", normal_res$name),1])

#************************** Genes overexpressed in CD34+ cancer **********************
# Genes expressed in CD34 cells:
# filter out StemCell expressed genes:
genes_CD34_2=unique(OE_cancer[OE_cancer$adj.pvalue>2&OE_cancer$FoldChange>1.5&grepl("_CD34", OE_cancer$name),])
genes_CD34=genes_CD34_2[!grepl("StemCell|CD8+Tcell|NKCell", genes_CD34_2$fails),]

#************************** Genes overexpressed in cancer **********************
overexpressed_cancer=OE_cancer[OE_cancer$adj.pvalue>2&OE_cancer$fails=="",]

# 1st condition, high in cancer (at least 5% patients), but low in normals
rm=unique(normal_res[normal_res$pvalue>=2,1]) # low expression in normals
overexpressed_cancer=overexpressed_cancer[!overexpressed_cancer[,1]%in%rm,]

#  2nd condition, expressed in normals, but highly expressed in cancer
overexpressed_cancer2=OE_cancer[OE_cancer$adj.pvalue>2&OE_cancer$FoldChange_max>=1.25&OE_cancer$fails=="",]
overexpressed_cancer2=overexpressed_cancer2[!overexpressed_cancer2[,1]%in%overexpressed_cancer[,1],]

# 3rd condition, low in both, but expressed in portion of cancer samples, has to be over 5% of the patients
overexpressed_cancer3=OE_cancer[OE_cancer$adj.pvalue>50&OE_cancer$FoldChange_max<1.25&!OE_cancer$fails==""&!OE_cancer[,1]%in%normal_res[normal_res$adj.pvalue>2,1],]
overexpressed_cancer3=overexpressed_cancer3[!(overexpressed_cancer3[,1]%in%overexpressed_cancer[,1]|overexpressed_cancer3[,1]%in%overexpressed_cancer2[,1]),]

#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
rm=unique(normal_res[normal_res$pvalue>=3,1]) # not expressed in normal sample category
allnormals=rownames(profile)[rowSums(profile[,annot$MAINCLASS%in%"NonCancerHealthy"]==1)>200] # might be missed from hgt test if expressed in lot of the samples

expressed_cancer=unique(OE_cancer[OE_cancer$adj.pvalue>2,1])
expressed_cancer=expressed_cancer[!(expressed_cancer%in%rm|expressed_cancer%in%allnormals)]

#****************** Genes that are expressed in over 5% of the pure cancer samples ******************
# CD34 positive, also add for the test
Sys.setlocale('LC_ALL','C')
cd34=get.logical(annovector = list(annot$colorClass), filterv = grepl("CD34", annot$Sample.isolation, ignore.case = T), PREFIX = "CD34")
cd34=cd34[!names(cd34)%in%c("Erythroid_CD34", "StemCell_CD34")]

# can use here all samples
subclass=get.logical(annovector = list(annot$subclasses), filterv = annot$Sample.type%in%c("Prolif", "Cancer"))
tbly=get.logical(annovector = list(annot$tbLY), filterv = annot$Sample.type%in%c("Prolif", "Cancer"))
tbly=tbly[names(tbly)%in%c("Lymphoma_BCL_CHL ", "Lymphoma_BCL_MCL", "Lymphoma_TCL_PTCLNOS", "Lymphoma_BCL_FL", "Lymphoma_BCL_MALT", "Lymphoma_TCL_AILT", "Lymphoma_TCL_ALCL", "Lymphoma_BCL_DLBCL_PMBL", "Lymphoma_TCL_CTCL_MF")]

list_cancers=lapply(c(subclass, tbly, cd34), list)
list_cancers=unlist(list_cancers, recursive=F)
list_cancers=list_cancers[sapply(list_cancers, sum)>20]

# also must be expressed in at least 5% of the patients, high enough
FUN_TEST <- function(gene, logicalVectors) {

  D = as.numeric(profile[rownames(profile)%in%gene,])

  test=sapply(logicalVectors, function(v){
    dd=D[v]
    res=(sum(dd[dd==1])/length(dd))
  })
  answ=data.frame(gene, "percentage"=as.numeric(test), "name"=names(test), stringsAsFactors = F)

  return(answ)
}
# filter out genes, expressed in any cancer over 5%
number_expressed=do.call(rbind, mclapply(rownames(profile), FUN_TEST, list_cancers, mc.cores=4))

expressed_cancer2=unique(number_expressed[as.numeric(number_expressed$percentage)>=0.05,1])
expressed_cancer3=unique(number_expressed[as.numeric(number_expressed$percentage)>=0.10,1])


#************************ read GTEX database data ************************
gtex=read.delim("GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct", header=T, stringsAsFactors = F, skip=2)

# process matrix
n=gtex[!duplicated(gtex[,2]),2]
gtex=gtex[!duplicated(gtex[,2]),3:dim(gtex)[2]]
rownames(gtex)=n

gtex=gtex[rownames(gtex)%in%rownames(profile),]

# test where the genes in hemap are expressed in GTEX
d=gtex[,colnames(gtex)%in%"Testis", drop=F]
d2=gtex[,!colnames(gtex)%in%c("Testis", "Ovary"), drop=F]

d=d[as.numeric(d[,1])>2,,drop=F]
filt=rownames(d2)[rowSums(d2>0.25)==0]
expressed_testis=rownames(d)[rownames(d)%in%filt]

# test where the genes in hemap are expressed in GTEX
d3=gtex
filt=rownames(d3)[rowSums(d3>0.25)==0]
gtex_notexpressed=rownames(d3)[rownames(d3)%in%filt]


#****************************************************************

# genelists to compare
genes_normals               # genes expressed in normal cells, filter out, unless overexpressed in cancer
expressed_testis            # germline genes
genes_CD34                  # these are good candidates
expressed_cancer            # these have to be expressed in cancer, not expressed in normals
expressed_cancer2           # expressed quite highly in at least 5%
expressed_cancer3           # expressed quite highly in at least 10%
overexpressed_cancer        # low normal expressed genes, must be over 10%
overexpressed_cancer2       # overexpressed genes
overexpressed_cancer3       # expressed in a portion of cancer cells
gtex_notexpressed           # not expressed in GTEX

# remove genes that are overexpressed from normals
genes_normals=genes_normals[!genes_normals%in%overexpressed_cancer|genes_normals%in%overexpressed_cancer2]

# remove normals from testis list
expressed_testis=expressed_testis[!expressed_testis%in%genes_normals]

# expressed_cancer is only interesting, if they are testis antigens
expressed_cancer_testis=expressed_cancer[expressed_cancer%in%expressed_testis&expressed_cancer%in%expressed_cancer2]

testis_df=number_expressed[number_expressed[,1]%in%expressed_cancer_testis&number_expressed[,2]>0.05,]

# overexpressed highly, or not expressed in normals and in 10% patients
overexpressed_cancer=overexpressed_cancer[overexpressed_cancer%in%expressed_cancer3,]

# These genes are expresed anywhere else
cancer_overexpressed_notGTEX=rbind(overexpressed_cancer[overexpressed_cancer[,1]%in%gtex_notexpressed,], overexpressed_cancer2[overexpressed_cancer2[,1]%in%gtex_notexpressed,])

geneLS=list(expressed_testis, genes_CD34[,1], expressed_cancer_testis, overexpressed_cancer[,1], overexpressed_cancer2[,1], cancer_overexpressed_notGTEX[,1])
geneLS2=list(expressed_testis, genes_CD34, expressed_cancer_testis, overexpressed_cancer)

# We can rename our list vectors
names(geneLS) <- c("expressed_testis", "genes_CD34", "expressed_cancer_testis", "overexpressed_cancer","overexpressed_cancer2", "cancer_overexpressed_notGTEX")
names(geneLS2) <- c("expressed_testis", "genes_CD34", "expressed_cancer_testis", "overexpressed_cancer")

a.df=rbind(overexpressed_cancer, overexpressed_cancer2, cancer_overexpressed_notGTEX, genes_CD34)

# remove linc and AS RNA
a.df=a.df[!grepl("LINC|-AS", a.df[,1]),]
testis_df=testis_df[!grepl("LINC|-AS", testis_df[,1]),]

# CCLE filter if  gene is not expressed:
ccle=get(load("CCLE_RNASEQ_SYMBOL_RPKM.Rdata"))

# CCLE annot:
annot=read.delim("CCLE_sample_info_file_2012-10-18.txt", header=T, stringsAsFactors = F)
find=intersect(annot$CCLE.name, colnames(ccle))
ccle=ccle[,match(find, colnames(ccle))]
annot=annot[match(find, annot$CCLE.name),]
ccle=ccle[,grepl("HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", colnames(ccle))]

ccle=ccle[rownames(ccle)%in%testis_df[,1],]

# at least 6 cell lines express over 1 rpkm
keep <- rowSums(ccle>0.5)>=5
unique(testis_df[,1])[!unique(testis_df[,1])%in%unique(testis_df[testis_df[,1]%in%rownames(ccle)[keep],1])]
testis_df=testis_df[testis_df[,1]%in%rownames(ccle)[keep],]

write.table(a.df, "antigen_df.txt", quote = F, row.names = F, col.names = T, sep="\t")
write.table(testis_df, "t.antigen_df.txt", quote = F, row.names = F, col.names = T, sep="\t")