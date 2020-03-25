#'
#'
plot.scRNA.features.complexHeatmap=function(s.anno, features, name,width=300, cores){
  log=mclapply(seq(dim(s.anno)[1]), plot.scRNA.features, s.anno,features, name, mc.cores = cores, width=width)
}

plot.scRNA.features=function(i, s.anno, features, name, add.annotation=NULL, width=300){
  # eccite PBMC dataset:
  load(s.anno[i,2])
  
  #HSC, myeloid, B-solu, T-solu, erythroi
  # STEP1 find the cluster order, and set them based on lineage:
  clusters=Idents(scmat)
  if(s.anno[i,3]!=""){
    order.clu=as.numeric(unlist(strsplit(s.anno[i,3], ",|, | ,"))) # MODIFY
    
    ord.names=names(clusters[order(match(clusters,order.clu))])
    
  }else{
    
    order.clu=levels(clusters)
    
    # order samples within cluster by umap component with more spread
    # this way there is still information within cluster, for heatmaps
    coords=Embeddings(scmat[["umap"]])
    fit.cluster=levels(Idents(scmat))
    
    ord.names=unlist(lapply(fit.cluster, function(clu){
      d=coords[Idents(scmat)%in%clu,]
      ind=which.max(c(max(d[,1])-min(d[,1]), max(d[,2])-min(d[,2]))) # take the axis with most spread
      names(sort(d[,ind]))
    }))
    
  }
  
  
  
  features=data.frame(features, stringsAsFactors = F)
  
  if(!grepl("^B:|^N:", features[1,1])){
    features[,1]=paste0("N:GEXP:", features[,1])
    allfeats=rownames(fm)[!rowSums(fm>0)==0]
    
  }else{
    allfeats=gsub("N:GEXP:|N:RPPA:", "", rownames(fm)[!rowSums(fm>0)==0])
  }
  
  # text annot created if more columns than 1.
  if(dim(features)[2]>1){
    features=features[features[,1]%in%allfeats,]
    feats=features[,1]
    t.annot=apply(features[,2:dim(features)[2], drop=F], 1, paste, collapse=" ")
    names(t.annot)=feats
  }else{
    t.annot=NULL
    feats=features[features[,1]%in%allfeats,1]
  }
  
  # make a data frame for annotations on top of the heatmap
  annotdf=data.frame("cluster"=clusters,"HLAIScore"=as.numeric(fm["N:SAMP:HLAIScore",]), "HLAIIScore"=as.numeric(fm["N:SAMP:HLAIIScore",]), "CytolyticScore"=as.numeric(fm["N:SAMP:CytolyticScore",]))
  
  if(!is.null(add.annotation))annotdf$annotation=add.annotation
  
  name=gsub("__", "_", gsub("[[:punct:]]", "_", paste(name, s.anno[i,1], sep="_")))
  
  # plot using a complexheatmap package using a featurematrix format and custom plotting function
  r1=plot.complexHM.fm(feats = feats, text_annot = t.annot, fm.f = fm, annotdf = annotdf, feats.barplot = c("HLAIScore", "HLAIIScore", "CytolyticScore"), order_columns=ord.names, order_rows=F, split.columns = T, plotting.param="/research/groups/sysgen/PROJECTS/students/HEMAP_IMMUNOLOGY/minna_work/scRNA_genelist_plotting/plotting_param_fm.txt", NAME=name, WIDTH = width)
}

plot.gene.seurat=function(s.anno, feature, name,cols=c("grey75", "red"), max=1, min=0, cores=1){
  plot.feature=function(i, feature){
    # plot proportions as a barplot next to 2D points
    load(s.anno[i,2])
    
    if(!feature%in%rownames(scmat))return(NULL)
    
    p <- FeaturePlot(scmat, features = feature, cols = cols, max.cutoff = max, min.cutoff = min, order = T, combine = FALSE)
    p=lapply(p, `+`, labs(title = s.anno[i,1]))
    
    if(i>1){
      for(i in 1:length(p)) {
        p[[i]] <- p[[i]] + NoLegend() + NoAxes()
      }
    }else{
      p[[1]] <- p[[1]] + NoAxes()
    }
    
    cowplot::plot_grid(plotlist = p)
  }
  
  # plot here lineage proportion per patient
  plot.list=mclapply(seq(dim(s.anno)[1]), plot.feature, feature, mc.cores=cores)
  plot.list=plot.list[!sapply(plot.list, is.null)]
  
  if(!length(plot.list))return(NULL)
  
  # 74mm * 99mm per panelwidth = 77, height = 74
  figure <-multipanelfigure:: multi_panel_figure(width = 170, height = ceiling(length(plot.list)/2)*75+75, rows = ceiling(length(plot.list)/2)+1, columns = 2, panel_label_type = "none", unit = "mm", row_spacing = unit(1, "mm"), column_spacing = unit(4, "mm"))
  
  for(i in seq(plot.list)){
    figure <- multipanelfigure::fill_panel(figure,  plot.list[[i]])
  }
  
  multipanelfigure::save_multi_panel_figure(figure, filename=paste0(name, "_fm_scatterplot.pdf"), limitsize = FALSE)
}

plot.fm.feature=function(s.anno, feature, name, max=3, min=0, cores=1){
  plot.feature=function(i, feature){
    # plot proportions as a barplot next to 2D points
    load(s.anno[i,2])
    
    scmat[[feature]]=fm[rownames(fm)%in%feature,]
    
    p <- FeaturePlot(scmat, features = feature, cols = c("grey75", "red"), max.cutoff = max, min.cutoff = min, order = T, combine = FALSE)
    p=lapply(p, `+`, labs(title = s.anno[i,1]))
    
    if(i>1){
      for(i in 1:length(p)) {
        p[[i]] <- p[[i]] + NoLegend() + NoAxes()
      }
    }else{
      p[[1]] <- p[[1]] + NoAxes()
    }
    
    return(p)
  }
  
  # plot here lineage proportion per patient
  plot.list=unlist(mclapply(seq(dim(s.anno)[1]), plot.feature, feature, mc.cores=cores), recursive = F)
  
  # 74mm * 99mm per panelwidth = 77, height = 74
  figure <-multipanelfigure:: multi_panel_figure(width = 170, height = ceiling(length(plot.list)/2)*75+75, rows = ceiling(length(plot.list)/2)+1, columns = 2, panel_label_type = "none", unit = "mm", row_spacing = unit(1, "mm"), column_spacing = unit(4, "mm"))
  
  for(i in seq(plot.list)){
    figure <- multipanelfigure::fill_panel(figure,  plot.list[[i]])
  }
  
  # 74mm * 99mm per panelwidth = 77, height = 74
  figure <-multipanelfigure:: multi_panel_figure(width = 170, height = ceiling(length(plot.list)/2)*75+75, rows = ceiling(length(plot.list)/2)+1, columns = 2, panel_label_type = "none", unit = "mm", row_spacing = unit(1, "mm"), column_spacing = unit(4, "mm"))
  
  for(i in seq(plot.list)){
    figure <- multipanelfigure::fill_panel(figure,  plot.list[[i]])
  }
  
  multipanelfigure::save_multi_panel_figure(figure, filename=paste0(name, "_fm_scatterplot.pdf"), limitsize = FALSE)
}

plot.several.fm.feature=function(s.anno.object, features, name, min=0, max=5, cores=5){
  load(s.anno.object)
  
  features=features[features%in%rownames(fm)]
  
  plot.feature=function(i, features){
    # plot proportions as a barplot next to 2D points
    
    scmat[[features[i]]]=fm[rownames(fm)%in%features[i],]
    
    p <- FeaturePlot(scmat, features = features[i], cols = c("grey75", "red"), max.cutoff = max, min.cutoff = min, order = T, combine = FALSE)
    p=lapply(p, `+`, labs(title = features[i]))
    
    if(i>1){
      for(i in 1:length(p)) {
        p[[i]] <- p[[i]] + NoLegend() + NoAxes()
      }
    }else{
      p[[1]] <- p[[1]] + NoAxes()
    }
    return(p)
  }
  
  # plot here lineage proportion per patient
  plot.list=unlist(mclapply(seq(features), plot.feature, features, mc.cores=cores), recursive = F)
  
  # 74mm * 99mm per panelwidth = 77, height = 74
  figure <-multipanelfigure:: multi_panel_figure(width = 170, height = ceiling(length(plot.list)/2)*75+75, rows = ceiling(length(plot.list)/2)+1, columns = 2, panel_label_type = "none", unit = "mm", row_spacing = unit(1, "mm"), column_spacing = unit(4, "mm"))
  
  for(i in seq(plot.list)){
    figure <- multipanelfigure::fill_panel(figure,  plot.list[[i]])
  }
  
  multipanelfigure::save_multi_panel_figure(figure, filename=paste0(name, "_fm_scatterplot.pdf"), limitsize = FALSE)
}

plot.seurat.ident=function(s.anno, colorv, name, cores=1){
  plot.feature=function(i){
    # plot proportions as a barplot next to 2D points
    load(s.anno[i,2])
    
    cells=gsub(" \\d+$", "", levels(Idents(scmat)))   
    
    p=DimPlot(scmat, group.by = c("ident"), cols = colorv[match(cells, colorv[,1]),2], ncol = 1, label = T, repel = T) + NoLegend() + NoAxes()
  }
  
  # plot here lineage proportion per patient
  plot.list=mclapply(seq(dim(s.anno)[1]), plot.feature, mc.cores=cores)
  
  # 74mm * 99mm per panelwidth = 77, height = 74
  figure <-multipanelfigure:: multi_panel_figure(width = 170, height = ceiling(length(plot.list)/2)*75+75, rows = ceiling(length(plot.list)/2)+1, columns = 2, panel_label_type = "none", unit = "mm", row_spacing = unit(1, "mm"), column_spacing = unit(4, "mm"))
  
  for(i in seq(plot.list)){
    figure <- multipanelfigure::fill_panel(figure,  plot.list[[i]])
  }
  
  multipanelfigure::save_multi_panel_figure(figure, filename=paste0(name, "_fm_scatterplot.pdf"), limitsize = FALSE)
}


#'
plot.multi.scatter.scRNA=function(s.anno, scmat=NULL, colors.group, identityVector.samples.list=NULL, seurat.feature=NULL, name="scRNA", cores=7, SIZE=0.1, rasterize=T, width=77, height = 74, text.size=10){
  # tools:
  library(Matrix)
  library(Seurat)
  source("/research/users/ppolonen/git_home/common_scripts/visualisation/plotting_functions.R")
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  library(parallel)
  library(ggplot2)
  
  
  # plot here lineage proportion per patient
  plot.list=mclapply(seq(dim(s.anno)[1]), function(i){
    # plot proportions as a barplot next to 2D points
    load(s.anno[i,2])
    coords=Embeddings(scmat[["umap"]])
    identityVector=unlist(strsplit(s.anno$Type[i], ", |,"))
    clusters=unlist(strsplit(s.anno$Clusters[i], ", |,"))
    clusters.samples=as.character(Idents(scmat))
    
    if(is.null(identityVector.samples.list)){
      identityVector.samples=clusters.samples
      
      for(j in seq(clusters)){
        identityVector.samples[clusters.samples%in%clusters[j]]=identityVector[j]
      }
    }else{
      identityVector.samples=as.character(identityVector.samples.list[[i]]) # must be in a list
    }
    
    if(!is.null(seurat.feature)){
      identityVector.samples=as.character(scmat[[seurat.feature]][,1]) # must be in a list
    }
    
    name.sample=s.anno$Sample[i]
    
    # take the colors:
    colorv=colors.group[match(unique(identityVector.samples), colors.group[,1]),2]
    
    a=plot.scatter(x = coords[,1], y = coords[,2], group = identityVector.samples, colorv = colorv, namev=name, main = name, add.proportions = T, SIZE = SIZE, add.legend=F, rasterize=rasterize, width=width, height = height, text.size=text.size)
    return(a)
  }, mc.cores=cores)
  
  # add legend last
  add.legend=get.only.legend(group = colors.group[,1], colorv = colors.group[,2])
  
  nrows=ceiling(length(plot.list)/2)+1
  # 74mm * 99mm per panelwidth = 77, height = 74
  figure <-multipanelfigure:: multi_panel_figure(width = 2*width+width*0.2, height = ceiling(length(plot.list)/2)*height+120, rows = nrows, columns = 2, panel_label_type = "none", unit = "mm", row_spacing = unit(1, "mm"), column_spacing = unit(4, "mm"))
  
  rowi=1
  for(i in seq(plot.list)){
    figure <- multipanelfigure::fill_panel(figure,  plot.list[[i]])
  }
  
  figure <- multipanelfigure::fill_panel(figure,row = nrows,column = 1:2,  add.legend)
  multipanelfigure::save_multi_panel_figure(figure, filename=paste0(name, "_fm_scatterplot.pdf"), limitsize = FALSE)
}

plot.scatter.seurat=function(scmat, colors.group,clusters="", identityVector="",  identityVector.samples.list=NULL, seurat.feature=NULL, lv=NULL, name="scRNA", cores=7, SIZE=0.1, rasterize=T, width=77, height=74, text.size=10, add.proportions=T, add.density=F, add.labels=F){
  # tools:
  library(Matrix)
  library(Seurat)
  source("/research/users/ppolonen/git_home/common_scripts/visualisation/plotting_functions.R")
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  library(parallel)
  library(ggplot2)
  
  coords=Embeddings(scmat[["umap"]])
  clusters.samples=as.character(Idents(scmat))
  
  if(is.null(identityVector.samples.list)){
    identityVector.samples=clusters.samples
    
    for(j in seq(clusters)){
      identityVector.samples[clusters.samples%in%clusters[j]]=identityVector[j]
    }
  }else{
    identityVector.samples=as.character(identityVector.samples.list[[i]]) # must be in a list
  }
  
  if(!is.null(seurat.feature)){
    identityVector.samples=as.character(scmat[[seurat.feature]][,1]) # must be in a list
    
    ind=order(match(identityVector.samples, colors.group[,1]))
    
    # order everything:
    identityVector.samples=identityVector.samples[ind]
    coords=coords[ind,]
    lv=lv[ind]
  }
  
  # take the colors:
  colorv=colors.group[match(unique(identityVector.samples), colors.group[,1]),2]
  
  plot.list=list(plot.scatter(x = coords[,1], y = coords[,2], group = identityVector.samples, colorv = colorv, lv = lv, namev=name, main = name, add.proportions = add.proportions, add.density = add.density, add.labels = add.labels, SIZE = SIZE, add.legend=T, rasterize=rasterize, width=width, height = height, text.size=text.size))

  figure <-multipanelfigure:: multi_panel_figure(width = width+50, height = height, rows = 1, columns = 1, panel_label_type = "none", unit = "mm", row_spacing = unit(1, "mm"), column_spacing = unit(4, "mm"))
  
  for(i in seq(plot.list)){
    figure <- multipanelfigure::fill_panel(figure,  plot.list[[i]])
  }
  
  multipanelfigure::save_multi_panel_figure(figure, filename=paste0(name, "_fm_scatterplot.pdf"), limitsize = FALSE)
}


#'
#'
compute.FDR.scRNA=function(i, s.anno, nperm=1000, geneset, nCores=6){
  library(AUCell)
  library(Matrix)
  library(Seurat)
  
  load(s.anno[i,2])
  
  # take only genes that are in the matrix
  geneset=geneset[geneset%in%rownames(scmat)]
  
  # Random
  geneset.random=lapply(seq(nperm), function(i, genes){
    sample(rownames(scmat), length(geneset))
  })
  
  names(geneset.random)=paste0("Random",seq(nperm))
  
  geneset.observed=list("MDS_signature"=geneset)
  
  data_n <- data.matrix(GetAssayData(scmat, assay = "RNA"))
  
  cells_rankings <- AUCell_buildRankings(data_n, nCores=nCores, plotStats=TRUE)
  AUC.random <- AUCell_calcAUC(geneset.random, cells_rankings, nCores = nCores, nrow(cells_rankings)*0.05)
  AUC.observed <- AUCell_calcAUC(geneset.observed, cells_rankings, nCores = nCores, nrow(cells_rankings)*0.05)
  
  # Alternative: assign cells according to the 'automatic' threshold
  cells_assignment <- AUCell_exploreThresholds(AUC.observed, plotHist=T, assignCells=TRUE)
  
  # how many times random higher than observed, divide by the number of permutations == eFDR
  y=as.numeric(getAUC(AUC.observed))
  x=getAUC(AUC.random)
  
  FDR=sapply(seq(y), function(i){
    sum(x[,i]>y[i])/nperm
  })
  
  empirical.FDR=data.frame(t(FDR), check.names = F)
  rownames(empirical.FDR)=paste0(names(geneset.observed), "_eFDR_", nperm, "_", length(geneset))
  colnames(empirical.FDR)=colnames(scmat)
  
  a=scmat[["SingleR.label"]]
  table(a[rownames(a)%in%cells_assignment$MDS_signature$assignment,])
  table(a[empirical.FDR<0.05,])
  
  scmat[["empirical.FDR"]]=as.numeric(empirical.FDR)
  DimPlot(scmat, group.by = "empirical.FDR")
  
  scmat[["y"]]=as.numeric(y)>0.04
  DimPlot(scmat, group.by = "y")
  
  return(list(empirical.FDR, y))
}

get.only.legend=function(colorv, group, position="right", direction="vertical"){
  library(ggplot2)
  
  dat=data.frame(group, colorv)
  levels(dat$group)=group
  myColors <- colorv
  names(myColors) <- levels(group)
  colScale <- scale_colour_manual(name = "group",values = myColors)
  dat$x=dim(dat)[1]
  dat$y=dim(dat)[1]
  
  legend <- cowplot::get_legend(ggplot(dat, aes(x, y, colour = group)) + theme_bw() +
                                  geom_point(size=0.6) + colScale +theme(legend.position = position, legend.direction = direction, legend.title = element_blank(), 
                                                                         legend.key=element_blank(), legend.box.background=element_blank(),
                                                                         legend.text=element_text(size = 10, face = "bold", family="Helvetica"))+
                                  theme(plot.margin=unit(c(1,1,1,1), "mm"))+
                                  guides(colour = guide_legend(override.aes = list(size=3))))
  return(legend)
}


statistical.comparisons.scRNA=function(i, s.anno, genelist.statistical.analysis=NULL, test.use="wilcox"){
  # plot proportions as a barplot next to 2D points
  load(s.anno[i,2])
  coords=Embeddings(scmat[["umap"]])
  identityVector=gsub(" ", "", unlist(strsplit(s.anno$Type[i], ", |,")))
  clusters=gsub(" ", "", unlist(strsplit(s.anno$Clusters[i], ", |,")))
  clusters.samples=as.character(Idents(scmat))
  identityVector.samples=clusters.samples
  
  for(j in seq(clusters)){
    identityVector.samples[clusters.samples%in%clusters[j]]=identityVector[j]
  }
  
  name.sample=s.anno$Sample[i]
  
  if(is.null(genelist.statistical.analysis))genelist.statistical.analysis=rownames(scmat)
  
  # Find the DE genes per factor and within factor:
  statistical.analysis=do.call(rbind, mclapply(unique(identityVector.samples), function(id){
    clusters2test=unique(clusters.samples[identityVector.samples%in%id])
    
    # clusters2test vs Rest:
    markers=FindMarkers(scmat, ident.1 = clusters2test, group.by = "seurat_clusters", test.use = test.use,  logfc.threshold =0.01, only.pos = T, features = genelist.statistical.analysis)
    
    #DE in clusters within category
    cluster.markers=FindAllMarkers(scmat[,clusters.samples%in%clusters2test], test.use = test.use, logfc.threshold = 0.01,  only.pos = T, features = genelist.statistical.analysis)
    
    # make genelists with annotation:
    markers$test=paste0(id, " vs. Rest")
    markers$gene=rownames(markers)
    cluster.markers$test=paste0(id, " ", cluster.markers$cluster, " vs. Rest")
    cluster.markers=cluster.markers[,!colnames(cluster.markers)%in%"cluster"]
    df.statistics=plyr::rbind.fill(list(markers, cluster.markers))
    df.statistics$test.used=test.use
    return(df.statistics)
  }))
  
  return(statistical.analysis)
}


#' @param seurat.object Seurat class object or data matrix that can be converted to Seurat object. Can also be a list of data matrices (RNA+citeseq etc.), Must be named by assay type with four uppercase letters (PROT, etc) (also found from seurat object with this name). 
#' @param regress.cell.label Character vector with as many columns as seurat.object Contains batch information for each cell. Uses fastMNN/SCTtransform to batch correct the data (usually different samples or experiments).
#' @param check.pcs Plot jackstraw and elbowplot to optimize the number of PCA componens
#' @param plot.umap Umap plot for clusters and for sample ID, if batch corrected also batch clustering
#' @param resolution Resolution parameter for clustering
#' @param nFeature.min Filter cells with < nFeature.min
#' @param nFeature.max Filter cells with > nFeature.max
#' @param percent.mitoDNA Filter cells with mitochondrial genes > percent.mitoDNA
#' @param singleR Annotate each cell using singleR built in annotations. Mainly uses Encode+blueprint annotations.
#' @param cores Number of cores to use for Seurat/singleR/MNNcorrect
#' @param add.gm.score Adds geneset score (using geometric mean) to fm, must be a named list of genes. Geneset name is added to fm and seurat object.
#' @param ... Pass parameters to Seurat tools
sc.data.analysis=function(scmat, name="scRNA_object", nr.pcs=30, regress.cell.label=NULL, batch.correction.method=c("SCTtransform","MNNcorrect"), auto.order=F, normalize.input=T, check.pcs=T, plot.umap=T, resolution=0.5, nFeature.min=0, nFeature.max=10000, percent.mitoDNA=25, singleR=T, singleR.reference=c("blueprint_encode", "HPCA"), cores=6, add.gm.score=NULL, ...){
  
  # determine computational resources:
  future::plan("multiprocess", workers = cores)
  options(future.globals.maxSize= 5e+09)
  
  # set method for batch correction
  batch.correction.method=match.arg(batch.correction.method)
  
  # list of data matrix, can be citeseq data:
  if(class(scmat)=="list"){
    cat("Input is list, splitting to data matrix by type and adding to Seurat", sep="\n\n")
    
    list.additional=lapply(2:length(scmat), function(i)CreateAssayObject(scmat[[i]]))
    names(list.additional)=names(scmat)[2:length(scmat)]
    scmat <- CreateSeuratObject(scmat[[1]])
  }else{
    list.additional=NULL
    # if class is not seurat object
    if(!class(scmat)=="Seurat"){
      cat("Input is data matrix, converting to Seurat object", sep="\n\n")
      scmat <- CreateSeuratObject(scmat)
    }else{
      cat("Input is Seurat object", sep="\n\n")
    }
  }
  DefaultAssay(object = scmat) <- "RNA"
  
  # is it a vector of pcs or single value?
  if(length(nr.pcs)==1)nr.pcs=seq(nr.pcs)
  
  # QC
  scmat=PercentageFeatureSet(scmat, pattern = "^MT-", col.name = "percent.mt")
  r1=VlnPlot(scmat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  ggplot2::ggsave(r1, filename = paste0(name, "_QC.pdf"), width = 12)
  
  # Filtering
  nFeature_RNA <- FetchData(object = scmat, vars = "nFeature_RNA")
  pct.mt <- FetchData(object = scmat, vars = "percent.mt")
  filt=which(x = nFeature_RNA > nFeature.min & nFeature_RNA < nFeature.max & pct.mt < percent.mitoDNA)
  cat(paste("Dimension before/after filtering:", paste(paste(dim(scmat), collapse="x"), paste(dim(scmat[, filt]), collapse="x"), collapse = "/")), sep="\n\n")
  scmat=scmat[, filt]
  
  if(normalize.input){
    assay="SCT"
    scmat=SCTransform(scmat, variable.features.n = 3000) #,...)
    
    # seurat standard processing
    scmat=RunPCA(scmat, npcs = 100)
    reduction="pca"
  }else{
    assay="RNA"
    
    # assume that the data is normalized
    scmat <- FindVariableFeatures(object = scmat)
    scmat <- ScaleData(object = scmat)
    scmat=RunPCA(scmat, npcs = 100)
    reduction="pca"
  }
  
  # set default assay to use:
  DefaultAssay(object = scmat) <- assay
  
  # use SCTransform to remove batch effect, used in umap and clustering
  # else use SCT to normalize data
  if(!is.null(regress.cell.label)&batch.correction.method=="SCTtransform"){
    batch.df=data.frame("batch"=regress.cell.label[filt])
    rownames(batch.df)=colnames(scmat)
    scmat[["batch"]] = batch.df
    scmat=SCTransform(scmat,vars.to.regress="batch", variable.features.n = 3000)
    
    # seurat standard processing
    scmat=RunPCA(scmat, npcs = 100)
    
    name=paste0(name, "_", batch.correction.method)
    reduction="pca"
  }
  
  # optimize pcs number:
  if(check.pcs){
    
    # setting here original RNA to define PCAs, SCT not working yet, wait for a fix
    # also the issue here how to deal with the mnnCorrect?
    # DefaultAssay(object = scmat) <- "RNA"
    # scmat <- NormalizeData(scmat,display.progress = FALSE)
    # scmat <- FindVariableFeatures(scmat,do.plot = F,display.progress = FALSE)
    # scmat <- ScaleData(scmat)
    # scmat=RunPCA(scmat, npcs = 100)
    
    # determine how many pcs should be used:
    r3=ElbowPlot(scmat, ndims = 100, reduction = "pca")
    ggplot2::ggsave(r3, filename = paste0(name, "_elbow.pdf"))
    
    # find optimal number of PCs to use:
    scmat <- JackStraw(scmat, num.replicate = 100, dims = 100, reduction = "pca")
    scmat <- ScoreJackStraw(scmat, dims = 1:100, reduction = "pca")
    r4=JackStrawPlot(scmat, dims = 1:100)
    ggplot2::ggsave(r4, filename = paste0(name, "_Jackstraw.pdf"), width = 12)
    
    cat("PCs checking done", sep="\n\n")
  }
  
  # use MNN to remove batch effect, used in umap and clustering:
  if(!is.null(regress.cell.label)&batch.correction.method=="MNNcorrect"){
    library(BiocParallel)
    library(batchelor)
    
    cat("FastMNN batch correction...", sep="\n\n")
    
    if(is.factor(regress.cell.label)&!auto.order){
      order=levels(regress.cell.label)
      cat("FastMNN batch correction order:", sep="\n\n")
      cat(order, sep="-->")
      cat("", sep="\n\n")
    }
    
    if(!is.factor(regress.cell.label)&!auto.order){
      order=unique(regress.cell.label)
      cat("FastMNN batch correction order:", sep="\n\n")
      cat(order, sep="-->")
      cat("", sep="\n\n")
    }
    
    batch.df=data.frame("batch"=regress.cell.label[filt])
    rownames(batch.df)=colnames(scmat)
    scmat[["batch"]] = batch.df
    
    # http://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/fast_mnn.html
    scmat <- NormalizeData(scmat)
    scmat <- FindVariableFeatures(scmat)
    
    object.list = SplitObject(scmat, split.by = "batch")[order(order)]
    
    scmat <- SeuratWrappers::RunFastMNN(object.list, reduction.name = "mnnCorrect", auto.order=auto.order, BPPARAM=MulticoreParam(cores))
    
    # set defaults 
    reduction="mnnCorrect"
    cat("FastMNN batch correction... done", sep="\n\n")
    
  }
  
  # Basic seurat processing with umap:
  scmat=RunUMAP(scmat, dims = nr.pcs, reduction = reduction)
  scmat=FindNeighbors(scmat, dims = nr.pcs, reduction = reduction)
  scmat=FindClusters(scmat, resolution = resolution)
  
  # add celltype automatically
  if(singleR){
    cat("Annotating with singleR", sep="\n\n")
    
    if(!file.exists(paste0(name, "_singler_object.Rdata"))){
      library("SingleR")
      
      counts <- GetAssayData(scmat, assay = assay, slot = "counts")
      
      if(singleR.reference=="blueprint_encode"){
        if(dim(counts)[2]<100000){
          singler = CreateSinglerObject(counts, annot = NULL, clusters=Idents(scmat), ref.list = list(blueprint_encode), project.name=name, fine.tune = T,do.signatures = T, numCores=cores)
        }else{
          cat("SingleR in batches", sep="\n\n")
          singler = CreateBigSingleRObject.custom(counts, annot = NULL, clusters=Idents(scmat),  N=30000, ref.list = list(blueprint_encode), xy = Embeddings(scmat[["umap"]]), project.name=name, fine.tune = T, do.signatures = T, numCores=cores)
        }
      }
      
      if(singleR.reference=="HPCA"){
        if(dim(counts)[2]<100000){
          singler = CreateSinglerObject(counts, annot = NULL, clusters=Idents(scmat), ref.list = list(hpca), project.name=name, fine.tune = T,do.signatures = T, numCores=cores)
        }else{
          cat("SingleR in batches", sep="\n\n")
          singler = CreateBigSingleRObject.custom(counts, annot = NULL, clusters=Idents(scmat), N=30000, ref.list = list(hpca), project.name=name, xy = Embeddings(scmat[["umap"]]), fine.tune = T,do.signatures = T, numCores=cores)
        }
      }
      
      # add original identifiers
      singler$meta.data$orig.ident = scmat@meta.data$orig.ident # the original identities, if not supplied in 'annot'
      
      ## if using Seurat v3.0 and over use:
      singler$meta.data$xy = scmat@reductions$umap@cell.embeddings # the tSNE coordinates
      singler$meta.data$clusters = scmat@active.ident # the Seurat clusters (if 'clusters' not provided)
      
      # add annot from hpca
      # singler2 = CreateSinglerObject(counts, annot = NULL, clusters=Idents(scmat), ref.list = list(hpca), project.name=name, fine.tune = T,do.signatures = T, numCores=cores)
      # singler$singler[[1]]$SingleR.single$labels[singler2$singler[[1]]$SingleR.single$labels%in%"Pre-B_cell_CD34-"]="pre-B-cells"
      # singler$singler[[1]]$SingleR.single$labels[singler2$singler[[1]]$SingleR.single$labels%in%"Pro-B_cell_CD34+"]="pro-B-cells"
      # singler$singler[[1]]$SingleR.clusters$labels[singler2$singler[[1]]$SingleR.clusters$labels%in%"Pre-B_cell_CD34-"]="pre-B-cells"
      # singler$singler[[1]]$SingleR.clusters$labels[singler2$singler[[1]]$SingleR.clusters$labels%in%"Pro-B_cell_CD34+"]="pro-B-cells"
      
      # these names are not intuitive, replace them, mentioned that these are actually naive http://comphealth.ucsf.edu/SingleR/SupplementaryInformation2.html:
      singler$singler[[1]]$SingleR.single$labels[singler$singler[[1]]$SingleR.single$labels%in%"CD4+ T-cells"]="CD4+ Tn"
      singler$singler[[1]]$SingleR.single$labels[singler$singler[[1]]$SingleR.single$labels%in%"CD8+ T-cells"]="CD8+ Tn"
      singler$singler[[1]]$SingleR.clusters$labels[singler$singler[[1]]$SingleR.clusters$labels%in%"CD4+ T-cells"]="CD4+ Tn"
      singler$singler[[1]]$SingleR.clusters$labels[singler$singler[[1]]$SingleR.clusters$labels%in%"CD8+ T-cells"]="CD8+ Tn"
      
      save(singler, file=paste0(name, "_singler_object.Rdata"))
      
    }else{
      library("SingleR")
      load(paste0(name, "_singler_object.Rdata"))
      cat(paste0("Loaded existing SingleR object (", paste0(name, "_singler_object.Rdata"), "), remove/rename it if you want to re-compute."), sep="\n\n")
      
    }
    # can be added as cluster id
    cluster=singler$singler[[1]]$SingleR.clusters$labels
    cluster.main=singler$singler[[1]]$SingleR.clusters.main$labels
    
    # order based on cell type and add celltype to cluster:
    lineage=c("HSC","MPP","CLP","CMP","GMP","MEP","Monocytes","DC","Macrophages","Macrophages M1","Macrophages M2","CD4+ Tn","CD4+ T-cells", "CD4+ Tcm", "CD4+ Tem","CD8+ Tn","CD8+ T-cells","CD8+ Tcm","CD8+ Tem","NK cells","Tregs","naive B-cells","Memory B-cells","Class-switched memory B-cells","Plasma cells","Endothelial cells","Neutrophils","Eosinophils","Fibroblasts","Smooth muscle","Erythrocytes","Megakaryocytes")
    lineage=unique(c(lineage,blueprint_encode$types, hpca$types))
    lineage=lineage[lineage%in%singler$singler[[1]]$SingleR.single$labels]
    
    scmat[["SingleR.label.main"]]=singler$singler[[1]]$SingleR.single.main$labels
    scmat[["SingleR.label.cluster"]]=paste(singler$singler[[1]]$SingleR.single$labels, Idents(scmat))
    scmat[["SingleR.label"]]=factor(singler$singler[[1]]$SingleR.single$labels, levels=lineage[lineage%in%unique(singler$singler[[1]]$SingleR.single$labels)])
    
    identityVector.samples=as.character(Idents(scmat))
    clusters.samples=Idents(scmat)
    
    for(j in seq(cluster)){
      identityVector.samples[clusters.samples%in%rownames(cluster)[j]]=paste0(cluster[j], ":", rownames(cluster)[j])
    }
    
    # order and cluster identity;
    cluster=cluster[order(match(cluster[,1],lineage)),,drop=F]
    cat(paste(levels(Idents(scmat)), collapse=","), paste(cluster, collapse=","), sep="\t\t")   
    
    scmat[["SingleR.cluster"]]=gsub(":.*.", "", Idents(scmat)) # just the "cluster label" per cell
    
    # order seurat clusters too:
    Idents(scmat)=factor(identityVector.samples, levels=paste0(cluster, ":", rownames(cluster)))
    
    r4=DimPlot(object = scmat, group.by = "SingleR.label", reduction = 'umap' , label = TRUE,  pt.size = 0.5, repel = T)  + NoLegend() + NoAxes()
    ggplot2::ggsave(r4, filename = paste0(name, "_UMAP_singleR.pdf"), width = 6, height = 6)
    
    cat("\nSingleR done", sep="\n\n")
    
    
    # new singleR:
    # counts <- GetAssayData(scmat, assay = "RNA", slot = "counts")
    # 
    # library(SingleR)
    # library(scater)
    # 
    # ref=NovershternHematopoieticData()	
    # 
    # common <- intersect(rownames(counts), rownames(ref))
    # ref <- ref[common,]
    # counts <- counts[common,]
    # counts.log <- logNormCounts(counts)
    # 
    # singler <- SingleR(test = counts.log, ref = ref, labels = ref$label.main, numCores=cores)
    # 
    # # order based on cell type and add celltype to cluster:
    # lineage=c("HSC","MPP","CLP","CMP","GMP","MEP","Monocytes","DC","Macrophages","Macrophages M1","Macrophages M2","CD4+ Tn","CD4+ T-cells", "CD4+ Tcm", "CD4+ Tem","CD8+ Tn","CD8+ T-cells","CD8+ Tcm","CD8+ Tem","NK cells","Tregs","naive B-cells","Memory B-cells","Class-switched memory B-cells","Plasma cells","Endothelial cells","Neutrophils","Eosinophils","Fibroblasts","Smooth muscle","Erythrocytes","Megakaryocytes")
    # lineage=unique(c(lineage,ref$label.main))
    # 
    # scmat[["SingleR.label.main"]]=singler$first.labels
    # scmat[["SingleR.label"]]=factor(singler$labels, levels=lineage[lineage%in%unique(singler$labels)])
    # 
    # r4=DimPlot(object = scmat, group.by = "SingleR.label", reduction = 'umap' , label = TRUE,  pt.size = 0.5, repel = T)  + NoLegend() + NoAxes()
    # ggplot2::ggsave(r4, filename = paste0(name, "_UMAP_singleR.pdf"), width = 6, height = 6)
    # 
    # cat("\nSingleR done", sep="\n\n")
    
  }
  
  # plot clusters and samples
  if(plot.umap){
    if(!is.null(regress.cell.label)){
      r6=DimPlot(scmat, group.by = c("batch"), ncol = 2, label = T, repel = T)  + NoLegend() + NoAxes()
      ggplot2::ggsave(r6, filename = paste0(name, "_UMAP_batch.pdf"), width = 6, height = 6)
      r6=DimPlot(scmat, group.by = c("ident"), ncol = 2, label = T, repel = T)  + NoLegend() + NoAxes()
      ggplot2::ggsave(r6, filename = paste0(name, "_UMAP_clusters.pdf"), width = 6, height = 6)
    }else{
      r6=DimPlot(scmat, group.by = c("ident"), ncol = 1, label = T, repel = T)  + NoLegend() + NoAxes()
      ggplot2::ggsave(r6, filename = paste0(name, "_UMAP_clusters.pdf"), width = 6, height = 6)
    }
  }
  
  # add immunoscores to FM format data matrix
  fm.f <- GetAssayData(scmat)
  
  # add gm based score:
  add.scores=list(HLAIScore=c("B2M", "HLA-A", "HLA-B","HLA-C"), HLAIIScore=c("HLA-DMA","HLA-DMB","HLA-DPA1","HLA-DPB1","HLA-DRA","HLA-DRB1"), CytolyticScore=c("GZMA", "PRF1", "GNLY", "GZMH", "GZMM"))
  if(!is.null(add.gm.score))add.scores=append(add.scores, add.gm.score)
  
  gm.objects=do.call(rbind, lapply(seq(add.scores), function(i){
    dat3=fm.f[rownames(fm.f)%in%add.scores[[i]],]
    gm=t(apply(dat3, 2, gm_mean)) # done to normalized values
    rownames(gm)=names(add.scores)[i]
    return(gm)
  }))
  
  # also add to seurat object:
  for(i in seq(add.scores)){
    scmat[[names(add.scores)[i]]] <- gm.objects[i,]
  }
  
  # gene expression and scores to fm:
  rownames(fm.f)=paste0("N:GEXP:", rownames(fm.f))
  rownames(gm.objects)=paste0("N:SAMP:", rownames(gm.objects))
  fm=rbind(gm.objects, fm.f)
  
  # DE analysis:
  markers.all=FindAllMarkers(scmat, only.pos = T)#, ...)
  
  # go through each data matrix, add to seurat and fm if more matrices
  if(!is.null(list.additional)){
    
    for(i in seq(list.additional)){
      mat.add=list.additional[[i]]
      mat.add=subset(mat.add,cells = colnames(scmat))
      scmat[[names(list.additional)[i]]] <- mat.add
      scmat <- NormalizeData(scmat, assay = names(list.additional)[i], normalization.method = "CLR") # test scttransform too
      scmat <- ScaleData(scmat, assay = names(list.additional)[i])
      
      add.data.normalized <- data.matrix(GetAssayData(scmat, assay = names(list.additional)[i], slot = "data"))
      rownames(add.data.normalized)=paste0("N:",names(list.additional)[i], ":", rownames(add.data.normalized))
      fm=rbind(fm, add.data.normalized)
    }
  }
  
  save(list=c("fm", "scmat", "markers.all"), file=paste0(name, "_scRNA.Rdata"))
  
  return(paste0(name, "_scRNA.Rdata"))
}

# make sure  data is on a linear non-log scale
# check 0 values, if a lot between 0-1 better to use add.one. If counts, remove works too quite well
# zero fix methods: https://arxiv.org/pdf/1806.06403.pdf
gm_mean=function(x, na.rm=TRUE, zero.handle=c("remove", "add.one", "zero.propagate")){
  zero.handle=match.arg(zero.handle)
  if(any(x < 0, na.rm = TRUE)){
    warning("Negative values produced NaN - is the data on linear - non-log scale?")
    return(NaN)
  }
  
  if(zero.handle=="remove"){
    return(exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x)))
  }
  
  if(zero.handle=="add.one"){
    return(exp(mean(log(x+1), na.rm = na.rm))-1)
  }
  
  if(zero.handle=="zero.propagate"){
    if(any(x == 0, na.rm = TRUE)){
      return(0)
    }
    return(exp(mean(log(x), na.rm = na.rm)))
  }
  
}

get.cluster.label.singleR=function(singler){
  
  # how to add this?
  cluster=singler$singler[[1]]$SingleR.clusters$labels
  cluster.main=singler$singler[[1]]$SingleR.clusters.main$labels
  
  # order based on cell type:
  # cat(paste0("'", colors.group[,1], "'"), sep=",")
  lineage=c('HSC','MPP','CLP','CMP','GMP','MEP','Monocytes','DC','Macrophages','Macrophages M1','Macrophages M2','CD4+ T-cells','CD8+ T-cells','Tregs','T_cell:gamma-delta','T_cell:CD4+_Naive','T_cell:CD8+_naive','T_cell:Treg:Naive','CD4+ Tcm','CD4+ Tem','CD8+ Tcm','CD8+ Tem','NK cells','naive B-cells','Memory B-cells','B-cells','Class-switched memory B-cells','Plasma cells','Endothelial cells','Neutrophils','Eosinophils','Fibroblasts','Smooth muscle','Erythroblast','Erythrocytes','Megakaryocytes')
  
  cluster=cluster[order(match(cluster[,1],lineage)),,drop=F]
  
  # order and cluster identity;
  cat(paste(rownames(cluster), collapse=","), paste(cluster, collapse=","), sep="\t")
  cat("", sep="\n")
  return(cluster)
}

# run cellphoneDB
run.cellphonedb=function(scmat, celltype, out=getwd(), cores=10, threshold=0.05){
  
  if(!dir.exists(out))dir.create(out, recursive = T)
  
  if(is.null(celltype)){
    celltype=gsub(" ", "_", as.character(Idents(scmat)))
  }
  
  setwd(out)
  # input1:
  df=data.frame("Cell"=gsub("\\.", "_", colnames(scmat)), "cell_type"=celltype, stringsAsFactors = F)
  df$cell_type[scmat[["SingleR.label"]]=="Tregs"]="Tregs"
  
  meta=file.path(out, "meta.tsv")
  write.table(df, meta, row.names = F, col.names = T, quote = F, sep="\t")
  
  # input2:
  counts <- data.matrix(GetAssayData(scmat, assay = "SCT", slot = "data"))
  countsfile=tempfile()
  
  # gene symbol was not working...
  library("org.Hs.eg.db") # remember to install it if you don't have it already
  ensambleid <- mapIds(org.Hs.eg.db, keys = rownames(counts), keytype = "SYMBOL", column="ENSEMBL")
  rownames(counts)=ensambleid
  
  counts=counts[!is.na(rownames(counts)),]
  
  countsfile=file.path(out, "counts.tsv")
  
  write.table(t(c("Gene", gsub("\\.", "_", colnames(counts)))), countsfile, row.names = F, col.names = F, quote = F, sep="\t")
  write.table(counts, countsfile, row.names = T, col.names = F, quote = F, sep="\t", append=T)
  
  system(paste0("cellphonedb method statistical_analysis ", meta, " ", countsfile, " --iterations=1000 --threads=", cores, " --threshold=", threshold))
  
  # make heatmap of interactions:
  system(paste0("cellphonedb plot heatmap_plot ", meta))
  
  deconvoluted=data.table::fread(file.path(out, "out/deconvoluted.txt"), data.table = F)
  pvalues=data.table::fread(file.path(out, "out/pvalues.txt"), data.table = F)
  significant_means.txt=data.table::fread(file.path(out, "out/significant_means.txt"), data.table = F)
  means=data.table::fread(file.path(out, "out/means.txt"), data.table = F)
  
  unlink(countsfile)
  unlink(meta)
  
  return(list("deconvoluted"=deconvoluted, "pvalues"=pvalues, "significant_means"=significant_means.txt, "means"=means))
}

# found bug in reference list, was missing, also min genes was 200 not 0 causing errors...
CreateBigSingleRObject.custom=function (counts, annot = NULL, project.name, xy, clusters, N = 10000,
                                        min.genes = 0, technology = "10X", species = "Human", citation = "",
                                        ref.list = list(), normalize.gene.length = F, variable.genes = "de",
                                        fine.tune = T, reduce.file.size = T, do.signatures = F, do.main.types = T,
                                        temp.dir = getwd(), numCores = SingleR.numCores){
  n = ncol(counts)
  s = seq(1, n, by = N)
  dir.create(paste0(temp.dir, "/singler.temp/"), showWarnings = FALSE)
  for (i in s) {
    print(i)
    A = seq(i, min(i + N - 1, n))
    
    singler = CreateSinglerObject(counts[, A], annot = annot[A], ref.list = ref.list,
                                  project.name = project.name, min.genes = min.genes,
                                  technology = technology, species = species, citation = citation,
                                  do.signatures = do.signatures, clusters = NULL, numCores = numCores)
    
    
    save(singler, file = paste0(temp.dir, "/singler.temp/",
                                project.name, ".", i, ".RData"))
  }
  singler.objects.file <- list.files(paste0(temp.dir, "/singler.temp/"),
                                     pattern = "RData", full.names = T)
  singler.objects = list()
  for (i in 1:length(singler.objects.file)) {
    load(singler.objects.file[[i]])
    singler.objects[[i]] = singler
  }
  singler = SingleR.Combine.custom(singler.objects, order = colnames(counts),expr = counts, clusters = clusters, xy = xy)
  singler
}

# bug also in this...
SingleR.Combine.custom=function (singler.list, order = NULL, clusters = NULL, expr = NULL, xy = NULL) 
{
  singler = c()
  singler$singler = singler.list[[1]]$singler
  for (j in 1:length(singler.list[[1]]$singler)) {
    singler$singler[[j]]$SingleR.cluster = c()
    singler$singler[[j]]$SingleR.cluster.main = c()
    singler$singler[[j]]$SingleR.single$clusters = c()
  }
  singler$meta.data = singler.list[[1]]$meta.data
  singler$meta.data$clusters = c()
  singler$meta.data$xy = c()
  singler$meta.data$data.sets = rep(singler$meta.data$project.name, 
                                    length(singler$meta.data$orig.ident))
  for (i in 2:length(singler.list)) {
    for (j in 1:length(singler$singler)) {
      if (singler.list[[i]]$singler[[j]]$about$RefData != 
          singler.list[[1]]$singler[[j]]$about$RefData) {
        stop("The objects are not ordered by the same reference data.")
      }
      singler$singler[[j]]$about$Organism = c(singler$singler[[j]]$about$Organism, 
                                              singler.list[[i]]$singler[[j]]$about$Organism)
      singler$singler[[j]]$about$Citation = c(singler$singler[[j]]$about$Citation, 
                                              singler.list[[i]]$singler[[j]]$about$Citation)
      singler$singler[[j]]$about$Technology = c(singler$singler[[j]]$about$Technology, 
                                                singler.list[[i]]$singler[[j]]$about$Technology)
      singler$singler[[j]]$SingleR.single$labels = rbind(singler$singler[[j]]$SingleR.single$labels, 
                                                         singler.list[[i]]$singler[[j]]$SingleR.single$labels)
      if (!is.null(singler$singler[[j]]$SingleR.single$labels1)) {
        singler$singler[[j]]$SingleR.single$labels1 = rbind(singler$singler[[j]]$SingleR.single$labels1, 
                                                            singler.list[[i]]$singler[[j]]$SingleR.single$labels1)
      }
      singler$singler[[j]]$SingleR.single$scores = rbind(singler$singler[[j]]$SingleR.single$scores, 
                                                         singler.list[[i]]$singler[[j]]$SingleR.single$scores)
      singler$singler[[j]]$SingleR.single.main$labels = rbind(singler$singler[[j]]$SingleR.single.main$labels, 
                                                              singler.list[[i]]$singler[[j]]$SingleR.single.main$labels)
      if (!is.null(singler$singler[[j]]$SingleR.single.main$labels1)) {
        singler$singler[[j]]$SingleR.single.main$labels1 = rbind(singler$singler[[j]]$SingleR.single.main$labels1, 
                                                                 singler.list[[i]]$singler[[j]]$SingleR.single.main$labels1)
      }
      singler$singler[[j]]$SingleR.single.main$scores = rbind(singler$singler[[j]]$SingleR.single.main$scores, 
                                                              singler.list[[i]]$singler[[j]]$SingleR.single.main$scores)
      singler$singler[[j]]$SingleR.single$cell.names = c(singler$singler[[j]]$SingleR.single$cell.names, 
                                                         singler.list[[i]]$singler[[j]]$SingleR.single$cell.names)
      singler$singler[[j]]$SingleR.single.main$cell.names = c(singler$singler[[j]]$SingleR.single.main$cell.names, 
                                                              singler.list[[i]]$singler[[j]]$SingleR.single.main$cell.names)
      if (!is.null(singler$singler[[j]]$SingleR.single.main$pval)) {
        singler$singler[[j]]$SingleR.single.main$pval = c(singler$singler[[j]]$SingleR.single.main$pval, 
                                                          singler.list[[i]]$singler[[j]]$SingleR.single.main$pval)
      }
      if (!is.null(singler$singler[[j]]$SingleR.single$pval)) {
        singler$singler[[j]]$SingleR.single$pval = c(singler$singler[[j]]$SingleR.single$pval, 
                                                     singler.list[[i]]$singler[[j]]$SingleR.single$pval)
      }
    }
    singler$meta.data$project.name = paste(singler$meta.data$project.name, 
                                           singler.list[[i]]$meta.data$project.name, sep = "+")
    singler$meta.data$orig.ident = c(singler$meta.data$orig.ident, 
                                     singler.list[[i]]$meta.data$orig.ident)
    singler$meta.data$data.sets = c(singler$meta.data$data.sets, 
                                    rep(singler.list[[i]]$meta.data$project.name, length(singler.list[[i]]$meta.data$orig.ident)))
  }
  for (j in 1:length(singler$singler)) {
    if (!is.null(order)) {
      singler$singler[[j]]$SingleR.single$labels = singler$singler[[j]]$SingleR.single$labels[order, 
                                                                                              ]
      if (!is.null(singler$singler[[j]]$SingleR.single$labels1)) {
        singler$singler[[j]]$SingleR.single$labels1 = singler$singler[[j]]$SingleR.single$labels1[order, 
                                                                                                  ]
      }
      singler$singler[[j]]$SingleR.single$scores = singler$singler[[j]]$SingleR.single$scores[order, 
                                                                                              ]
      singler$singler[[j]]$SingleR.single$cell.names = singler$singler[[j]]$SingleR.single$cell.names[order]
      singler$singler[[j]]$SingleR.single.main$labels = singler$singler[[j]]$SingleR.single.main$labels[order, 
                                                                                                        ]
      if (!is.null(singler$singler[[j]]$SingleR.single.main$labels1)) {
        singler$singler[[j]]$SingleR.single.main$labels1 = singler$singler[[j]]$SingleR.single.main$labels1[order, 
                                                                                                            ]
      }
      singler$singler[[j]]$SingleR.single.main$scores = singler$singler[[j]]$SingleR.single.main$scores[order, 
                                                                                                        ]
      singler$singler[[j]]$SingleR.single.main$cell.names = singler$singler[[j]]$SingleR.single.main$cell.names[order]
      if (!is.null(singler$singler[[j]]$SingleR.single$pval)) {
        singler$singler[[j]]$SingleR.single$pval = singler$singler[[j]]$SingleR.single$pval[order]
        singler$singler[[j]]$SingleR.single.main$pval = singler$singler[[j]]$SingleR.single.main$pval[order]
      }
    }
  }
  if (!is.null(clusters) && !is.null(expr)) {
    for (j in 1:length(singler$singler)) {
      if (is.character(singler$singler[[j]]$about$RefData)) {
        ref = get(tolower(singler$singler[[j]]$about$RefData))
      }
      else {
        ref = singler$singler[[j]]$about$RefData
      }
      singler$singler[[j]]$SingleR.clusters = SingleR("cluster", 
                                                      expr, ref$data, types = ref$types, clusters = factor(clusters), 
                                                      sd.thres = ref$sd.thres, genes = "de", fine.tune = T)
      singler$singler[[j]]$SingleR.clusters.main = SingleR("cluster", 
                                                           expr, ref$data, types = ref$main_types, clusters = factor(clusters), 
                                                           sd.thres = ref$sd.thres, genes = "de", fine.tune = T)
    }
    singler$meta.data$clusters = clusters
    if (!is.null(xy)) {
      singler$meta.data$xy = xy
    }
  }
  singler
}