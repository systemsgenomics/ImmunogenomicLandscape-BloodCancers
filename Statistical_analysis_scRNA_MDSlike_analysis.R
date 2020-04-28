GIT_HOME="/research/users/ppolonen/git_home/ImmunogenomicLandscape-BloodCancers/"
source(file.path(GIT_HOME, "common_scripts/scRNA/functions.scRNA.analysis.R"))
source(file.path(GIT_HOME, "common_scripts/visualisation/plotting_functions.R"))

library(Matrix)
library(Seurat)
library(data.table)
library(ComplexHeatmap)
library(circlize)
library(parallel)
library(ggplot2)

setwd("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/Published_data_figures")

# analyze FIMM and Galen AML and find costim associations to certain clusters.
co.stim=data.table::fread("costim_ligands_final.txt", data.table = F)[,c(1,3,5)]
co.stimR=data.table::fread("costim_ligands_final.txt", data.table = F)
co.stimR.f=c(unlist(strsplit(co.stimR$`Receptor gene`, ", ")), "LAG3")
co.stim=rbind(co.stim, c("FGL1", "LAG3", "Inhibitory"))

load("AML_Galen_scRNA.Rdata")
galen=scmat

load("FIMM_AML_scRNA.Rdata")

# get more detailed annotations for T and NK cells:
NK=get(load("Yang_HCA_AML_cluster_cell.Rdata"))
Tcell=get(load("Szabo_HCA_AML_cluster_cell.Rdata"))

anno=as.character(Idents(scmat))

# add the new annotation:
labels=lapply(unique(NK), function(s)names(NK)[NK%in%s])
names(labels)=unique(NK)

anno[names(Idents(scmat))%in%labels[[1]]]=paste("Mature NK cells")

# add the new annotation:
labels=lapply(unique(Tcell), function(s)names(Tcell)[Tcell%in%s])
names(labels)=unique(Tcell)

anno[names(Idents(scmat))%in%labels[[4]]]=paste("Cytotoxic CD8+ Tem")
anno[names(Idents(scmat))%in%labels[[1]]]=paste("Cytokine CD8+ Tem")

Idents(scmat)=factor(anno, levels=c(levels(Idents(scmat)), "Mature NK cells", "Cytotoxic CD8+ Tem", "Cytokine CD8+ Tem"))

# only take cells that need to be computed:
scmat2=scmat[,grepl("HSC|MPP|MEP|GMP|CMP|Monocyte|Mature NK cells|Cytotoxic CD8|Cytokine CD8", Idents(scmat))&scmat[["batch"]][,1]%in%c("5897", "3667", "5249")]

#************************************************** run cellphoneDB using these categories, focusing on interactions between blasts and NK/CD8 cells **************************************************
samples=unique(scmat2[["batch"]][,1])

# cellPhoneDB analysis (TableS3 tab for the results):
run=lapply(samples, function(s){
  results=run.cellphonedb(scmat2[,scmat2[["batch"]][,1]%in%s], celltype = NULL, out=file.path(getwd(),"cellphonedb", paste0(s, "_FIMM_AML_0.1")), threshold = 0.1, cores = 6)
})

#******************************************************************************************************************************************************************************************************


#************************************************** run DE analysis for each cluster and above defined NK-CD8 cells **************************************************

library(future)
plan("multiprocess", workers = 8)

# only do the analysis on most abundant clusters
scmat[["MDSlike"]]=ifelse(scmat[["batch"]][,1]%in%c("5897", "3667", "5249"), "MDS-like", "other")

addname=c(paste("MDS-like", levels(Idents(scmat))), paste("other", levels(Idents(scmat))))

Idents(scmat)=factor(paste(scmat[["MDSlike"]][,1], Idents(scmat)), levels=addname[addname%in%paste(scmat[["MDSlike"]][,1], Idents(scmat))])

markers.all.up=FindAllMarkers(scmat,  features = rownames(scmat),only.pos = T, logfc.threshold = 0.1, min.cells.group = 300)

markers.all.filt=markers.all.up[markers.all.up$p_val_adj<1e-10,]
markers.all.filt=markers.all.filt[markers.all.filt$avg_logFC>0.25,]
write.table(markers.all.filt[,c(7,6,2,3,4,1,5)], "TableS3_Significant_genes_MDSlike_scRNA_all.txt", quote = F, row.names = F, col.names = T, sep="\t")
save(markers.all.filt, file="markers.all.filt.allgenes.Rdata")

#******************************************************************************************************************************************************************************************************

# pathway analysis
load("Combined_pathway_signatures_2017_filtered_robust.Rdata")
genesets=listA[grepl("HALLMARKS|REACTOME|BIOCARTA|PID|KEGG|WIKIPW", names(listA))]

group=unique(markers.all.filt$cluster)

res.pw.up=mclapply(group, function(g){
  up=markers.all.filt$gene[markers.all.filt$cluster%in%g]
  
  if(!length(up))return(NULL)
  
  pw.up=do.call(rbind, mclapply(seq(genesets), function(i){
    if(sum(genesets[[i]]%in%rownames(scmat))==0)return(NULL)
    r=fisher.2x2(lv1 = rownames(scmat)%in%genesets[[i]], lv2 = rownames(scmat)%in%up, alternative = "greater", usefast = F)
    data.frame("Name"=names(genesets)[i], r, "genes.signif"=paste(up[up%in%genesets[[i]]], collapse=","), stringsAsFactors = F)
  }, mc.cores=1))
  
  pw.up$cluster=g
  pw.up$FDR=p.adjust(pw.up$p, "BH")
  pw.up=pw.up[order(pw.up$FDR, decreasing = F),]
  
  pw.up=pw.up[pw.up$FDR<0.05,]
  
  if(dim(pw.up)[1]==0)return(NULL)
  
  return(pw.up)
}, mc.cores=8)

res.pw.all=do.call(rbind, res.pw.up)

write.table(res.pw.all[,c(1,5,3,2,6,4)], "TableS3_MDSlike_DE_pathways_up.txt", quote = F, row.names = F, sep="\t")
