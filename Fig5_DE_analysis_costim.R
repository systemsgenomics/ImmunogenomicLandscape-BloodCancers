GIT_HOME="/research/users/ppolonen/git_home/ImmunogenomicLandscape-BloodCancers/"
source(file.path(GIT_HOME, "common_scripts/visualisation/plotting_functions.R"))
source(file.path(GIT_HOME, "common_scripts/statistics/functions_statistics.R"))
source(file.path(GIT_HOME, "common_scripts/statistics/statistics_wrappers.R"))
source(file.path(GIT_HOME, "common_scripts/pathway_analysis/functions.GSEA.R"))
source(file.path(GIT_HOME, "common_scripts/scRNA/functions.scRNA.analysis.R"))
source(file.path(GIT_HOME, "common_scripts/statistics/useful_functions.R"))

setwd("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/Published_data_figures")

# analyze FIMM and Galen AML and find costim associations to certain clusters.
co.stim=data.table::fread("costim_ligands_final.txt", data.table = F)[,c(1,3,5)]
co.stim=rbind(co.stim, c("FGL1", "LAG3", "Inhibitory"))
co.stimR.f=unlist(strsplit(co.stim$`Receptor gene`, ", "))

co.stim.all=data.table::fread("costim_ligands_final.txt", data.table = F)
ligand.name=gsub(" \\(\\)", "", paste0(co.stim.all$Gene, " (", co.stim.all$`Common name`, ")"))
Receptor.name=gsub(" \\(\\)", "", paste0(co.stim.all$`Receptor gene`, " (", co.stim.all$`Receptor common name`, ")"))
Receptor.name[co.stim.all$`Receptor gene`==co.stim.all$`Receptor common name`]=co.stim.all$`Receptor gene`[co.stim.all$`Receptor gene`==co.stim.all$`Receptor common name`]

rec=strsplit(co.stim$`Receptor gene`, ", ")
co.stimR.func=unique(do.call(rbind, lapply(seq(rec), function(i)cbind(rec[[i]], co.stim$`Immune checkpoint function`[i]))))

files=list.files(path = ".", "subtypes.Rdata")
names(files)=gsub("_subtypes.Rdata", "", files)

# plot subtypes
wrapper.de.analysis=function(i, files){
  load(files[i])
  library(ggplot2)
  
  if(!is.null(dim(coordinates.subtype))){
    subtype=factor(coordinates.subtype$subtype)
    
    coordinates.subtype=coordinates.subtype[order(subtype),]
    subtype=subtype[order(subtype)]
    
    SIZE=ifelse(dim(gexp)[2]>500, 0.5, 2)
    
    # # plot colorv for all:
    # colorv = data.frame("subtype"=c("PML-RARA", "CBFB-MYH11", "Progenitor-like", "MDS-like", "CEBPA", "Monocyte-like","Monocyte-like-MLL", "RUNX1-RUNX1T1"), "color"=c("#aff558", "#dfe4c3", "#d3c684", "#fac75d", "#8a9979", "#d4bcc5", "brown", "#edb2c6"))
    # 
    # colorv=colorv[match(unique(subtype), colorv[,1]),2]
    # 
    # plot.scatter(x=coordinates.subtype$x, y = coordinates.subtype$y, group = subtype, namev = subtype, main = names(files)[i],colorv =as.character(colorv), rasterize = F, width = 70*2, height = 74*2, SIZE = SIZE)
    
    # plot colorv for all:
    # filt=as.character(subtype)%in%c("Ph", "KMT2A", "TCF3-PBX1","Hyperdiploid","ETV6-RUNX1", "PAX5 P80R", "PAX5alt","DUX4", "ZNF384", "MEF2D", "Other")
    # 
    # coordinates.subtype=coordinates.subtype[filt,]
    # subtype=subtype[filt]
      
    plot.scatter(x=coordinates.subtype$x, y = coordinates.subtype$y, group = subtype, namev = subtype, main = names(files)[i], rasterize = F, width = 70*2, height = 74*2, SIZE = SIZE)
    
  }else{
    return(NULL)
  }
}
plot.list=lapply(seq(files), wrapper.de.analysis, files)
plot.list=plot.list[!sapply(plot.list, is.null)]
# 74mm * 99mm per panelwidth = 77, height = 74
figure <-multipanelfigure:: multi_panel_figure(width = 170*2, height = ceiling(length(plot.list)/2)*75*2+175, rows = ceiling(length(plot.list)/2)+1, columns = 2, panel_label_type = "none", unit = "mm", row_spacing = unit(1, "mm"), column_spacing = unit(4, "mm"))

for(i in seq(plot.list)){
  figure <- multipanelfigure::fill_panel(figure,  plot.list[[i]])
}

name="FigS5A_FigS3D_Fig6E_FigS6F_subtypes"
multipanelfigure::save_multi_panel_figure(figure, filename=paste0(name, "_scatterplot.pdf"), limitsize = FALSE)

files=files[!files%in%"PecanALL_subtypes.Rdata"]

#************************************** Differential expression ****************************************

wrapper.de.analysis=function(i, files, genelist){
  load(files[i])
  
  if(!is.null(dim(coordinates.subtype))){
    subtype=coordinates.subtype$subtype
  }else{
    subtype=coordinates.subtype
  }
  
  # make lv of the subtype:
  lv=get.logical(list(subtype))
  
  res=wrapper.wilcoxtest(genelist[genelist%in%rownames(gexp)], data = gexp, logicalVectors = lv, ALTERNATIVE = "two.sided", adj.method = "BH", CORES = 8, prettynum = F)
  
  res$Name=names(files)[i]
  # res=res[res$FDR<0.05,]
  return(res)
}

genelist=co.stim[,1]

res=lapply(seq(files), wrapper.de.analysis, files, genelist)

genelist=co.stimR.f

resR=lapply(seq(files), wrapper.de.analysis, files, genelist)

wrapper.immunoscore.analysis=function(i, files, genesets){
  load(files[i])
  
  if(!is.null(dim(coordinates.subtype))){
    subtype=coordinates.subtype$subtype
  }else{
    subtype=coordinates.subtype
  }
  
  # make lv of the subtype:
  lv=get.logical(list(subtype))
  
  # scores
  gm.objects=do.call(rbind, lapply(seq(genesets), function(i){
    dat3=2^gexp[rownames(gexp)%in%genesets[[i]],,drop=F]+0.01
    gm=log2(t(apply(dat3, 2, gm_mean))) # done to normalized values
    rownames(gm)=names(genesets)[i]
    return(gm)
  }))
  
  
  res=wrapper.wilcoxtest(rownames(gm.objects), data = gm.objects, logicalVectors = lv, ALTERNATIVE = "two.sided", adj.method = "BH", CORES = 8, prettynum = F)
  
  res$Name=names(files)[i]
  # res=res[res$FDR<0.05,]
  return(res)
}
# compute cytolytic and HLA scores
genesets=list(HLAIScore=c("B2M", "HLA-A", "HLA-B","HLA-C"), HLAIIScore=c("HLA-DMA","HLA-DMB","HLA-DPA1","HLA-DPB1","HLA-DRA","HLA-DRB1"), CytolyticScore=c("GZMA", "PRF1", "GNLY", "GZMH", "GZMM"))

res.immunoscores=lapply(seq(files), wrapper.immunoscore.analysis, files, genesets)

# diseases to test:
AML=grep("AML", files)
ALL=grep("ALL", files)
MM=grep("MM", files)
DLBCL=grep("DLBCL", files)

# all results
AML.res=do.call(rbind, res[AML])
ALL.res=do.call(rbind, res[ALL])
ALL.res$Gene[ALL.res$Gene=="VSIR"]="C10orf54"
MM.res=do.call(rbind, res[MM])
DLBCL.res=do.call(rbind, res[DLBCL])

# all results
AML.resR=do.call(rbind, resR[AML])
ALL.resR=do.call(rbind, resR[ALL])
MM.resR=do.call(rbind, resR[MM])
DLBCL.resR=do.call(rbind, resR[DLBCL])

# all results
AML.res.immunoscores=do.call(rbind, res.immunoscores[AML])
ALL.res.immunoscores=do.call(rbind, res.immunoscores[ALL])
MM.res.immunoscores=do.call(rbind, res.immunoscores[MM])
DLBCL.res.immunoscores=do.call(rbind, res.immunoscores[DLBCL])

get.summary.df=function(gene, res.df, stars=T){
  a=sapply(unique(res.df$Group1), function(subtype){
    mean(res.df$FC[res.df$Gene%in%gene&res.df$Group1%in%subtype])
  })
  b=sapply(unique(res.df$Group1), function(subtype){
    survcomp::combine.test(res.df$FDR[res.df$Gene%in%gene&res.df$Group1%in%subtype], method = "z.transform")
  })
  
  pval=-log10(b)
  if(stars){
    pval[b>0.05]=0
    pval[b<0.05]=1
    pval[b<0.01]=2
    pval[b<0.001]=3
    pval[b<1e-5]=4
    pval[b<1e-16]=6
  }
  
  data.frame("variable.1"=log2(a), "variable.2"=pval, "features"=gene, "id"=unique(res.df$Group1))
}

genes.costim=c('CD80', 'CD86', 'CD274', 'PDCD1LG2', 'CD276','C10orf54', 'ICOSLG', 'TNFRSF14', 'PVR','PVRL2', 'PVRL3', 'LGALS9', 'CD48', 'CD58', 'SLAMF6', 'SLAMF7', 'CD84', 'LY9', 'TNFSF4', 'TNFSF8', 'TNFSF9', 'TNFSF18', 'TNFSF15', 'CD70', 'BTN1A1', 'BTN2A2', 'BTN3A1', 'BTNL2', 'BTNL8', 'CD72', 'CD200',  'MICA', 'MICB', 'ULBP1', 'ULBP2', 'ULBP3','RAET1E', 'CLEC2B', 'CLEC2D', 'VTCN1', 'HHLA2', 'IDO1', 'IDO2', 'TDO2', 'NT5E', 'ENTPD1', 'ARG1')

cat(ligand.name[match(genes.costim, co.stim$Gene)], sep="\n")
cat(Receptor.name[match(genes.costim, co.stim$Gene)], sep="\n")
cat(co.stim$`Immune checkpoint function`[match(genes.costim, co.stim$Gene)], sep="\n")

# summarize all the values costim:
df.AML=do.call(rbind, lapply(genes.costim, get.summary.df, AML.res))
df.ALL=do.call(rbind, lapply(genes.costim, get.summary.df, ALL.res))
df.MM=do.call(rbind, lapply(genes.costim, get.summary.df, MM.res))
df.DLBCL=do.call(rbind, lapply(genes.costim, get.summary.df, DLBCL.res))

# all together:
all.subtypes=c(c("MDS-like", "Progenitor-like", "Monocyte-like", "Monocyte-like-MLL", "CEBPA", "RUNX1-RUNX1T1", "CBFB-MYH11", "PML-RARA"),c("Ph", "KMT2A", "TCF3-PBX1","Hyperdiploid","ETV6-RUNX1", "PAX5 P80R", "PAX5alt","DUX4", "ZNF384", "MEF2D", "Other"),  c("ABC", "GCB"), c("Hyperdiploid_gain1q", "Hyperdiploid_gain11q", "CCND1_Ig", "MAF_Ig", "WHSC1_FGFR3_Ig", "TRAF3_Aberrated"))
all.df=rbind(df.AML, df.ALL, df.MM, df.DLBCL)
all.df=all.df[all.df$id%in%all.subtypes,]
all.df=all.df[order(match(as.character(all.df$id), all.subtypes)),]

# summarize all the values, costim R:
df.AML.R=do.call(rbind, lapply(unique(co.stimR.f), get.summary.df, AML.resR))
df.ALL.R=do.call(rbind, lapply(unique(co.stimR.f), get.summary.df, ALL.resR))
df.MM.R=do.call(rbind, lapply(unique(co.stimR.f), get.summary.df, MM.resR))
df.DLBCL.R=do.call(rbind, lapply(unique(co.stimR.f), get.summary.df, DLBCL.resR))

# summarize all the values, immunoscores:
df.AML.immunoscores=do.call(rbind, lapply(rev(unique(AML.res.immunoscores[,1])), get.summary.df, AML.res.immunoscores))
df.ALL.immunoscores=do.call(rbind, lapply(rev(unique(AML.res.immunoscores[,1])), get.summary.df, ALL.res.immunoscores))
df.MM.immunoscores=do.call(rbind, lapply(rev(unique(AML.res.immunoscores[,1])), get.summary.df, MM.res.immunoscores))
df.DLBCL.immunoscores=do.call(rbind, lapply(rev(unique(AML.res.immunoscores[,1])), get.summary.df, DLBCL.res.immunoscores))

all.df.immunoscores=rbind(df.AML.immunoscores, df.ALL.immunoscores, df.MM.immunoscores, df.DLBCL.immunoscores)
all.df.immunoscores=all.df.immunoscores[all.df.immunoscores$id%in%all.subtypes,]
all.df.immunoscores=all.df.immunoscores[order(match(as.character(all.df.immunoscores$id), all.subtypes)),]

all.df$features=factor(all.df$features, levels=unique(all.df$features))
all.df$id=factor(all.df$id, levels=unique(all.df$id))
all.df.immunoscores$features=factor(all.df.immunoscores$features, levels=unique(all.df.immunoscores$features))
all.df.immunoscores$id=factor(all.df.immunoscores$id, levels=unique(all.df.immunoscores$id))

plot.dat=rbind(all.df.immunoscores, all.df)

# Final:
pdf("FigS5B_dotplot.pdf", width = 6, height = 6)
plot.DotPlot.df(data.plot = plot.dat, name.variable.1 = "Fold-Change (log2)", name.variable.2 = "FDR (-log10)", cols = c("blue", "white","red"), col.min = -2, col.max = 2, scale.min = 1, scale.max = 6, dot.scale = 2.5, number.legend.points = 6, fontsize = 7)
dev.off()


df.AML=df.AML[order(-df.AML$variable.2, df.AML$id),]
df.AML$features=factor(df.AML$features, levels = unique(as.character(df.AML$features))[order(df.AML$variable.1)])
df.AML$id=factor(df.AML$id, levels = unique(as.character(df.AML$id)))

# find for each subtype genes that are upregulated/downregulated:
find=ALL.res$FDR<0.05&abs(log2(ALL.res$FC))>0.15
ALL.res=ALL.res[find,]

# find for each subtype genes that are upregulated/downregulated:
find=MM.res$FDR<0.05&abs(log2(MM.res$FC))>0.15
MM.res=MM.res[find,]

# find for each subtype genes that are upregulated/downregulated:
find=DLBCL.res$FDR<0.05&abs(log2(DLBCL.res$FC))>0.15
DLBCL.res=DLBCL.res[find,]

# write tables for supplement:
write.table(AML.res, file = "TableS5_AML_costim_DEgenes.txt", col.names = T, row.names = F, quote = F, sep="\t")
write.table(ALL.res, file = "TableS5_ALL_costim_DEgenes.txt", col.names = T, row.names = F, quote = F, sep="\t")
write.table(MM.res, file = "TableS5_MM_costim_DEgenes.txt", col.names = T, row.names = F, quote = F, sep="\t")
write.table(DLBCL.res, file = "TableS5_DLBCL_costim_DEgenes.txt", col.names = T, row.names = F, quote = F, sep="\t")


# diseases to test:
AML=grep("AML", files)
ALL=grep("ALL", files)
MM=grep("MM", files)
DLBCL=grep("DLBCL", files)

# find for each subtype genes that are upregulated/downregulated:
find=AML.resR$FDR<0.05&abs(log2(AML.resR$FC))>0.15
AML.resR=AML.resR[find,]

# find for each subtype genes that are upregulated/downregulated:
find=ALL.resR$FDR<0.05&abs(log2(ALL.resR$FC))>0.15
ALL.resR=ALL.resR[find,]

# find for each subtype genes that are upregulated/downregulated:
find=MM.resR$FDR<0.05&abs(log2(MM.resR$FC))>0.15
MM.resR=MM.resR[find,]

# find for each subtype genes that are upregulated/downregulated:
find=DLBCL.resR$FDR<0.05&abs(log2(DLBCL.resR$FC))>0.15
DLBCL.resR=DLBCL.resR[find,]

# write tables for supplement:
# write.table(AML.resR, file = "AML_costimR_DEgenes.txt", col.names = T, row.names = F, quote = F, sep="\t")
# write.table(ALL.resR, file = "ALL_costimR_DEgenes.txt", col.names = T, row.names = F, quote = F, sep="\t")
# write.table(MM.resR, file = "MM_costimR_DEgenes.txt", col.names = T, row.names = F, quote = F, sep="\t")
# write.table(DLBCL.resR, file = "DLBCL_costimR_DEgenes.txt", col.names = T, row.names = F, quote = F, sep="\t")



# similar analysis to immunoscores:
# find for each subtype genes that are upregulated/downregulated:
find=AML.res.immunoscores$FDR<0.05&abs(log2(AML.res.immunoscores$FC))>0.15
AML.res.immunoscores=AML.res.immunoscores[find,]

# find for each subtype genes that are upregulated/downregulated:
find=ALL.res.immunoscores$FDR<0.05&abs(log2(ALL.res.immunoscores$FC))>0.15
ALL.res.immunoscores=ALL.res.immunoscores[find,]

# find for each subtype genes that are upregulated/downregulated:
find=MM.res.immunoscores$FDR<0.05&abs(log2(MM.res.immunoscores$FC))>0.15
MM.res.immunoscores=MM.res.immunoscores[find,]

# find for each subtype genes that are upregulated/downregulated:
find=DLBCL.res.immunoscores$FDR<0.05&abs(log2(DLBCL.res.immunoscores$FC))>0.15
DLBCL.res.immunoscores=DLBCL.res.immunoscores[find,]

# write tables for supplement:
# write.table(AML.res.immunoscores, file = "AML_immunoscores_DEgenes.txt", col.names = T, row.names = F, quote = F, sep="\t")
# write.table(ALL.res.immunoscores, file = "ALL_immunoscores_DEgenes.txt", col.names = T, row.names = F, quote = F, sep="\t")
# write.table(MM.res.immunoscores, file = "MM_immunoscores_DEgenes.txt", col.names = T, row.names = F, quote = F, sep="\t")
# write.table(DLBCL.res.immunoscores, file = "DLBCL_immunoscores_DEgenes.txt", col.names = T, row.names = F, quote = F, sep="\t")

extract.res=function(de.genes, geneannot){
  subtypes=de.genes$Group1
  
  de.genes.up=de.genes$FC>1
  de.genes.down=de.genes$FC<1
  
  res=lapply(unique(subtypes), function(type){
    
    counts.up=table(de.genes$Gene[subtypes%in%type&de.genes.up])
    counts.down=table(de.genes$Gene[subtypes%in%type&de.genes.down])
    
    res.up=list()
    res.down=list()
    
    # up
    if(sum(counts.up>1)>0){
      significant=names(counts.up)[counts.up>1]
      funct=geneannot[match(significant, geneannot[,1]),2]
      
      res.up=data.frame(significant, funct, direction="upregulated", stringsAsFactors = F)
      res.up=res.up[order(res.up$funct),]
    }
    
    # down
    if(sum(counts.down>1)>0){
      significant=names(counts.down)[counts.down>1]
      funct=geneannot[match(significant, geneannot[,1]),2]
      
      res.down=data.frame(significant, funct, direction="downregulated", stringsAsFactors = F)
      res.down=res.down[order(res.down$funct),]
    }
    
    res=rbind(res.up, res.down)   
  })
  
  names(res)=unique(subtypes)
  
  return(res)
}

# all results:
data.AML=extract.res(AML.res[AML.res$FDR<0.05&abs(log2(AML.res$FC))>0.15,], co.stim[,c(1,3)])
# sapply(names(data.AML),function (x) write.table(data[[x]], file=paste(x, "_AML_costim.txt", sep=""), quote = F, sep="\t", col.names = T, row.names = F))

data.ALL=extract.res(ALL.res[ALL.res$FDR<0.05&abs(log2(ALL.res$FC))>0.15,], co.stim[,c(1,3)])
# sapply(names(data.ALL)[!sapply(data.ALL, length)==0],function (x) write.table(data[[x]], file=paste(x, "_ALL_costim.txt", sep=""), quote = F, sep="\t", col.names = T, row.names = F))

data.MM=extract.res(MM.res, co.stim[,c(1,3)])
# sapply(names(data.MM),function (x) write.table(data[[x]], file=paste(x, "_MM_costim.txt", sep=""), quote = F, sep="\t", col.names = T, row.names = F))

data.DLBCL=extract.res(DLBCL.res[DLBCL.res$FDR<0.05&abs(log2(DLBCL.res$FC))>0.15,], co.stim[,c(1,3)])
# sapply(names(data.DLBCL),function (x) write.table(data[[x]], file=paste(x, "_DLBCL_costim.txt", sep=""), quote = F, sep="\t", col.names = T, row.names = F))

data.AML.R=extract.res(AML.resR, co.stimR.func)
# sapply(names(data.AML.R)[names(data.AML.R)%in%"MDS-like"],function (x) write.table(data[[x]], file=paste(x, "_AML_costimR.txt", sep=""), quote = F, sep="\t", col.names = T, row.names = F))

data.ALL.R=extract.res(ALL.resR[ALL.resR$FDR<0.05&abs(log2(ALL.resR$FC))>0.15,], co.stimR.func)
# sapply(names(data)[-15],function (x) write.table(data[[x]], file=paste(x, "_ALL_costimR.txt", sep=""), quote = F, sep="\t", col.names = T, row.names = F))

data=extract.res(MM.resR, co.stimR.func)
# sapply(names(data),function (x) write.table(data[[x]], file=paste(x, "_MM_costimR.txt", sep=""), quote = F, sep="\t", col.names = T, row.names = F))

data.DLBCL.R=extract.res(DLBCL.resR, co.stimR.func)
# sapply(names(data.DLBCL.R),function (x) write.table(data[[x]], file=paste(x, "_DLBCL_costimR.txt", sep=""), quote = F, sep="\t", col.names = T, row.names = F))

data=extract.res(AML.res.immunoscores, co.stim)
# sapply(names(data),function (x) write.table(data[[x]], file=paste(x, "_AML_immunoscores.txt", sep=""), quote = F, sep="\t", col.names = T, row.names = F))

data=extract.res(ALL.res.immunoscores, co.stim)
# sapply(names(data),function (x) write.table(data[[x]], file=paste(x, "_ALL_immunoscores.txt", sep=""), quote = F, sep="\t", col.names = T, row.names = F))

data=extract.res(MM.res.immunoscores, co.stim)
# sapply(names(data),function (x) write.table(data[[x]], file=paste(x, "_MM_immunoscores.txt", sep=""), quote = F, sep="\t", col.names = T, row.names = F))

data=extract.res(DLBCL.res.immunoscores, co.stim)
# sapply(names(data),function (x) write.table(data[[x]], file=paste(x, "_DLBCL_immunoscores.txt", sep=""), quote = F, sep="\t", col.names = T, row.names = F))


inhibitory=c("CytolyticScore", "HLAIIScore","HLAIScore", co.stim[grepl("Inhibitory|unknown", co.stim$`Immune checkpoint function`),1], co.stimR.func[grepl("Inhibitory|unknown", co.stimR.func[,2]),1])
stimulatory=c("CytolyticScore", "HLAIIScore","HLAIScore", co.stim[grepl("Stimulatory", co.stim$`Immune checkpoint function`),1], co.stimR.func[grepl("Stimulatory", co.stimR.func[,2]),1])
stimulatory=stimulatory[!stimulatory%in%"CTLA4"]

# AML.interesting:
subtype.order=c("MDS-like", "Progenitor-like", "Monocyte-like", "Monocyte-like-MLL", "CEBPA", "RUNX1-RUNX1T1", "CBFB-MYH11", "PML-RARA")
df.AML.pick=rbind(df.AML,df.AML.R, df.AML.immunoscores)
df.AML.pick$features=as.character(df.AML.pick$features)
df.AML.pick$id=as.character(df.AML.pick$id)
interesting.feats=unique(c("CytolyticScore", "HLAIIScore","PVRL3", "SLAMF7","HHLA2","CD274", "TNFRSF14", "ARG1","PVR","CD86", "C10orf54", "CD276","ENTPD1","NT5E", "CD48", "ICOSLG","TNFSF8", "LY9",  "BTN3A1","CD200",  "ULBP1", "LGALS9", "CD70", "MICB", "TNFSF9", "ULBP1", "ULBP3","CLEC2D","TNFSF9","CD84","LY9","CD58","CD72","LGALS9", "CD86", "CD48", "CLEC2B", "LY9", "TNFSF15", "TNFSF8","CD276", "PVR", "BTNL8", "CD58",
                    "LAG3", "SLAMF7","TIGIT","TMIGD2", "PDCD1", "CD160"))

# take all
interesting.feats.L=unlist(sapply(data.AML[subtype.order], function(d)d$significant[d$direction=="upregulated"]))
interesting.feats.R=unlist(sapply(data.AML.R[names(data.AML.R)%in%"MDS-like"], function(d)d$significant[d$direction=="upregulated"]))

interesting.feats=c("CytolyticScore", "HLAIIScore", interesting.feats.L, interesting.feats.R)

df.AML.pick=df.AML.pick[df.AML.pick$features%in%interesting.feats,]
df.AML.pick=df.AML.pick[order(match(df.AML.pick$id, subtype.order)),]
df.AML.pick=df.AML.pick[order(match(df.AML.pick$features, interesting.feats)),]
df.AML.pick$features=factor(df.AML.pick$features, levels=unique(df.AML.pick$features))
df.AML.pick$id=factor(df.AML.pick$id, levels=unique(df.AML.pick$id))

# ALL.interesting:
subtype.order=c("Ph", "KMT2A", "TCF3-PBX1","Hyperdiploid","ETV6-RUNX1", "PAX5 P80R", "PAX5alt","DUX4", "ZNF384", "MEF2D", "Other")
df.ALL.pick=rbind(df.ALL,df.ALL.R, df.ALL.immunoscores)
df.ALL.pick=df.ALL.pick[df.ALL.pick$id%in%subtype.order,]
df.ALL.pick$features=as.character(df.ALL.pick$features)
df.ALL.pick$id=as.character(df.ALL.pick$id)
interesting.feats=unique(c("CytolyticScore", "HLAIIScore","HLAIScore", 
                           "CD200", "NT5E", "TNFRSF14", "BTN3A1", "CLEC2B", "TNFSF4", "TNFSF8",
                           "CLEC2D", "CD58", "CD72","TNFSF8", "TNFSF9",
                           "CLEC2D","CD48", "CD58", "CD72","CD84",
                           "CD200", "CD86", "ICOSLG", "MICA",
                           "BTN2A2", "CD200", "NT5E", "BTN3A1", "TNFSF4",
                           "TNFRSF14", "CD48", "LY9", "TNFSF4",
                           "NT5E", "TNFRSF14", "CLEC2D", "CD48", "CD70", "CLEC2B", "LY9",
                           "BTN2A2", "LGALS9", "TNFRSF14", "CD84", "LY9", "TNFSF4", "TNFSF9",
                           "ENTPD1", "LGALS9", "TNFRSF14", "CD84", "CLEC2B", "TNFSF9",
                           "NT5E", "CLEC2D", "CD48", "CD72", "CD84",
                           "CD274", "CD80"
                           ))

interesting.feats.L=unlist(sapply(data.ALL[subtype.order], function(d)d$significant[d$direction=="upregulated"]))
interesting.feats.R=unlist(sapply(data.ALL.R[subtype.order], function(d)d$significant[d$direction=="upregulated"]))
interesting.feats=c("CytolyticScore", "HLAIIScore","HLAIScore", interesting.feats.L, interesting.feats.R)

df.ALL.pick=df.ALL.pick[df.ALL.pick$features%in%interesting.feats,]
df.ALL.pick=df.ALL.pick[order(match(df.ALL.pick$id, subtype.order)),]
df.ALL.pick=df.ALL.pick[order(match(df.ALL.pick$features, interesting.feats)),]
df.ALL.pick$features=factor(df.ALL.pick$features, levels=unique(df.ALL.pick$features))
df.ALL.pick$id=factor(df.ALL.pick$id, levels=unique(df.ALL.pick$id))

# MM.interesting:
subtype.order=c("Hyperdiploid_gain1q", "Hyperdiploid_gain11q", "MAF_Ig", "WHSC1_FGFR3_Ig", "TRAF3_Aberrated", "CCND1_Ig")
df.MM.pick=rbind(df.MM,df.MM.R, df.MM.immunoscores)
df.MM.pick$features=as.character(df.MM.pick$features)
df.MM.pick$id=as.character(df.MM.pick$id)
interesting.feats=unique(c("HLAIIScore","HLAIScore", "MICA", "MICB", 
                           "CD274","LGALS9", "CD200",  "CD48", "CD86", "ICOSLG", "PVRL3"))
interesting.feats=unlist(sapply(data.MM[subtype.order], function(d)d$significant[d$direction=="upregulated"]))
interesting.feats=c("HLAIIScore","HLAIScore", interesting.feats)

df.MM.pick=df.MM.pick[df.MM.pick$features%in%interesting.feats,]
df.MM.pick=df.MM.pick[order(match(df.MM.pick$id, subtype.order)),]
df.MM.pick=df.MM.pick[order(match(df.MM.pick$features, interesting.feats)),]
df.MM.pick$features=factor(df.MM.pick$features, levels=unique(df.MM.pick$features))
df.MM.pick$id=factor(df.MM.pick$id, levels=unique(df.MM.pick$id))

# DLBCL.interesting:
df.DLBCL.pick=rbind(df.DLBCL,df.DLBCL.R, df.DLBCL.immunoscores)
df.DLBCL.pick$features=as.character(df.DLBCL.pick$features)
df.DLBCL.pick$id=as.character(df.DLBCL.pick$id)
interesting.feats=unique(c("CytolyticScore", "HLAIIScore","HLAIScore","CD274","TNFRSF14", "ENTPD1","SLAMF7",  "ICOSLG", "LY9","CD86",
                           "PDCD1", "BTLA","CD160",  "SLAMF7", "LAG3"))
# take all
interesting.feats.L=unlist(sapply(data.DLBCL[c("ABC", "GCB")], function(d)d$significant[d$direction=="upregulated"]))
interesting.feats.R=unlist(sapply(data.DLBCL.R[c("ABC", "GCB")], function(d)d$significant[d$direction=="upregulated"]))

interesting.feats=c("CytolyticScore", "HLAIIScore", interesting.feats.L, interesting.feats.R)

df.DLBCL.pick=df.DLBCL.pick[df.DLBCL.pick$features%in%interesting.feats&df.DLBCL.pick$id%in%c("ABC", "GCB"),]
df.DLBCL.pick=df.DLBCL.pick[order(match(df.DLBCL.pick$id, c("ABC", "GCB"))),]
df.DLBCL.pick=df.DLBCL.pick[order(match(df.DLBCL.pick$features, interesting.feats)),]
df.DLBCL.pick$features=factor(df.DLBCL.pick$features, levels=unique(df.DLBCL.pick$features))
df.DLBCL.pick$id=factor(df.DLBCL.pick$id, levels=unique(df.DLBCL.pick$id))

# Most interesting shown in main figure, rest in FigS5B. These were significant in at least 2 data sets. Also checked expression in scRNA for various targets.
pdf("Fig5B_AML.pdf", height = 3.5, width = 4.25)
plot.DotPlot.df(data.plot = df.AML.pick[df.AML.pick$features%in%c("CytolyticScore", "HLAIIScore", "PVRL3", "SLAMF7", "HHLA2", "CD274", "TNFRSF14", "ARG1", "PVR", "CD86", "C10orf54", "CD276", "ENTPD1", "NT5E", "CLEC2B", "CD84", "CD48", "LAG3", "TIGIT", "SLAMF7", "TMIGD2", "PDCD1", "CD160", "CD2", "KLRF1"),], name.variable.1 = "Fold-Change (log2)", name.variable.2 = "FDR (-log10)", cols = c("blue", "white","red"), col.min = -2, col.max = 2, scale.min = 1, scale.max = 6, dot.scale = 2.5, number.legend.points = 6, fontsize = 8)
dev.off()

pdf("Fig5B_ALL.pdf", height = 1.5, width = 4.25)
plot.DotPlot.df(data.plot = df.ALL.pick[df.ALL.pick$features%in%c("CD274", "C10orf54", "NT5E"),], name.variable.1 = "Fold-Change (log2)", name.variable.2 = "FDR (-log10)", cols = c("blue", "white","red"), col.min = -2, col.max = 2, scale.min = 1, scale.max = 6, dot.scale = 2.5, number.legend.points = 6, fontsize = 8)
dev.off()

pdf("Fig5B_MM.pdf", height = 2.25, width = 3.8)
plot.DotPlot.df(data.plot = df.MM.pick[df.MM.pick$features%in%c("CytolyticScore", "HLAIIScore","MICA", "MICB","CD274","LGALS9", "CD48", "CD86", "ICOSLG", "PVRL3"),], name.variable.1 = "Fold-Change (log2)", name.variable.2 = "FDR (-log10)", cols = c("blue", "white","red"), col.min = -2, col.max = 2, scale.min = 1, scale.max = 6, dot.scale = 2.5, number.legend.points = 6, fontsize = 8)
dev.off()

pdf("Fig5B_DLBCL.pdf", height = 3.5, width = 3.5)
# all interesting:
plot.DotPlot.df(data.plot = df.DLBCL.pick, name.variable.1 = "Fold-Change (log2)", name.variable.2 = "FDR (-log10)", cols = c("blue", "white","red"), col.min = -2, col.max = 2, scale.min = 1, scale.max = 6, dot.scale = 3, number.legend.points = 6)
plot.DotPlot.df(data.plot = df.DLBCL.pick, name.variable.1 = "Fold-Change (log2)", name.variable.2 = "FDR (-log10)", cols = c("blue", "white","red"), dot.scale = 3)
dev.off()