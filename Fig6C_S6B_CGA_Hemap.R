GIT_HOME="/research/users/ppolonen/git_home/common_scripts/"
source(file.path(GIT_HOME, "visualisation/plotting_functions.R"))
source(file.path(GIT_HOME, "statistics/functions_statistics.R"))
library(RColorBrewer)
library(survival)
library(ComplexHeatmap)
library(circlize)

setwd("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/Published_data_figures")

t.df = read.delim("t.antigen_df.txt", stringsAsFactors=F, header=T)
t.df=t.df[order(t.df[,3]),]

annot = get(load("Hemap_immunology_Annotations.Rdata"))

# profiles
profile=get(load("mixtureM_profile.Rdata"))
data=t(get(load("data9544_with_gene_symbols.RData")))
data=data[,colnames(data)%in%annot$GSM.identifier..sample.]
profile2=profile[,colnames(profile)%in%annot$GSM.identifier..sample.]

# same filtering as when generating gene lists
profile2[profile2==-1] = 0
profile2[data<5]=0

# data:
genelist=unique(t.df[,1])

# rank patients by number of testis antigens expressed
expressed_testis_num=colSums(profile2[rownames(profile2)%in%genelist,])

# filtering, certain phenotypes, go through one by one.
logicalVectors=get.logical(annovector = list(annot$subclasses), filterv = !(annot$disease=="NonCancer"|annot$colorClass=="CellLine"))
logicalVectors=logicalVectors[c(6,7,3,5,11,9,8,10)]

annot$tbLY[annot2$tbLY%in%c("Lymphoma_BCL_DLBCL_ABC", "Lymphoma_BCL_DLBCL_GCB")]="Lymphoma_BCL_DLBCL"
logicalVectors2=get.logical(annovector = list(annot$tbLY), filterv = !(annot$disease=="NonCancer"|annot$colorClass=="CellLine")&annot$colorClass%in%c("BCL", "TCL"))
logicalVectors2=logicalVectors2[sapply(logicalVectors2, sum)>25]
logicalVectors=append(logicalVectors, logicalVectors2[c(1,7,9,5,8, 2:4)])
names(logicalVectors)=gsub("Prolif_Myeloproliferative_|_LC$|Lymphoma_BCL_|Lymphoma_TCL_|Cancer_", "", names(logicalVectors))
names(logicalVectors)=gsub("Myeloma", "MM", names(logicalVectors))

logicalVectors=logicalVectors[c(8:16,1:7)]

perc=sapply(seq(logicalVectors), function(i){
  lv=logicalVectors[[i]]
  n=names(logicalVectors)[i]
  
  signif(rowSums(profile2[rownames(profile2)%in%genelist, lv])/sum(lv),2)*100
})

colnames(perc)=names(logicalVectors)

# summarize each disease
data_anno=do.call(rbind, lapply(seq(logicalVectors), function(i){
  lv=logicalVectors[[i]]
  n=names(logicalVectors)[i]
  
  df=data.frame("GSMid"=annot$GSM.identifier..sample.[lv], "Disease"=n, "Number_Antigens_expressed"=expressed_testis_num[lv], "Cytolytic"=as.numeric(annot$CytScore)[lv], "HLAI"=as.numeric(annot$HLAIScore)[lv], "HLAII"=as.numeric(annot$HLAIIScore)[lv], "Age"=annot$AGE[lv], "OS"=annot$OS_Time[lv], "PFS"=annot$PFS_Time[lv], "Gender"=annot$GENDER[lv], stringsAsFactors = F)
  df=df[order(df$Number_Antigens_expressed),]
  
}))

# sort based on antigens expressed and disease
plot_data=data[match(genelist, rownames(data)), match(data_anno[,1], colnames(data))]
plot_data2=profile2[match(genelist, rownames(profile2)), match(data_anno[,1], colnames(profile2))]

ha = HeatmapAnnotation(df=data_anno[,2,drop=F], Number.Antigens = anno_barplot(data_anno$Number_Antigens_expressed, axis = T, gp = gpar(fill = "indianred")), height = unit(3, "inch"))
ha2 = HeatmapAnnotation(df = data_anno[,c(3:10)], height = unit(4, "inch"))

plot_data=t(scale(t(plot_data)))
plot_data[plot_data<(-1)]=-1
col_fun = colorRamp2(c(-4, 1, 2, 4), c("#579ec9", "white", "#ca4c41", "#67001f"))

ord=rowSums(plot_data2)

plot_data=plot_data[order(ord, decreasing = T),]
plot_data2=plot_data2[order(ord, decreasing = T),]
perc=perc[match(rownames(plot_data), rownames(perc)),]

colv=list(c("MM"="#924064", "DLBCL"="#d46c24","FL"="#c08e94","MCL"="#aa561e","MALT"="#c08e94","CHL"="#c08e94","ALCL"="#0a0091", "PTCLNOS"="#0a0091","AILT"="#0a0091","T-ALL"="#708dd0", "pre-B-ALL"="#d680c8", "AML"="#4ec173", "MDS"="#ff9d53", "LCH"="#ffea90", "CML"="#266c34", "CLL"="#8648ad"))

names(colv)="Disease"

ha = HeatmapAnnotation(df=data_anno[,2,drop=F], "n.CGA" = anno_barplot(data_anno$Number_Antigens_expressed, axis = T, gp = gpar(fill = "indianred")), height = unit(20, "mm"), col=colv)

data_anno=do.call(rbind, lapply(seq(logicalVectors), function(i){
  lv=logicalVectors[[i]]
  n=names(logicalVectors)[i]
  
  df=data.frame("GSMid"=annot$GSM.identifier..sample.[lv], "Disease"=n, "Number_Antigens_expressed"=expressed_testis_num[lv], "Cytolytic"=as.numeric(annot$CytScore)[lv], "HLAI"=as.numeric(annot$HLAIScore[lv])[lv], "HLAII"=as.numeric(annot$HLAIIScore)[lv], "Age"=annot$AGE[lv], "OS"=annot$OS_Time[lv], "PFS"=annot$PFS_Time[lv], "Gender"=annot$GENDER[lv], stringsAsFactors = F)
  df=df[order(df$Number_Antigens_expressed),]
  
}))

dp=t(table(data_anno$Disease, data_anno$Number_Antigens_expressed>=4))[2,]/table(data_anno$Disease)*100
dp=dp[match(names(logicalVectors), names(dp))]

topgenes=unlist(apply(perc, 2, function(v){
  v=v[order(v, decreasing = T)]
  v=names(v[v>0])
  if(length(v)>3)v=v[1:3]
  return(v)
}))


ha2 = HeatmapAnnotation(df=unique(data_anno[,2,drop=F]),">=4 CGA" = anno_barplot(data.frame(as.numeric(dp)), axis = T, gp = gpar(fill = colv$Disease), height = unit(10, "mm")), height = unit(20, "mm"), col=colv)

pdf("Fig6C_FigS6B_HEMAP_CGA_heatmap.pdf", height = 8, width = 10)
Heatmap(plot_data, top_annotation = ha, name = "# CGA", col = col_fun, cluster_columns = F, cluster_rows = F, show_column_names = F, row_names_gp = gpar(fontsize = 8), width = unit(50, "mm"), height = unit(2*dim(plot_data)[1], "mm"))
Heatmap(perc[unique(topgenes),], top_annotation = ha2, name = "% CGA", cluster_columns = F, cluster_rows = F, show_column_names = T, rect_gp = gpar(col = "white", lty = 1, lwd = 0.5), row_names_gp = gpar(fontsize = 8),column_names_gp = gpar(fontsize = 8), width = unit(50, "mm"), height = unit(1.75*dim(plot_data)[1], "mm"), col=colorRamp2(c(0,5,10,15,20,50,90), c("white", "#FCBBA1", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D", "#99000D")))
dev.off()


