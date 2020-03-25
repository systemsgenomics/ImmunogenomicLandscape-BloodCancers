
wrapper.GSEA=function(i, logicalVectors, data, cls.vector, GENESETS, datatype=NULL, clsname="", WD=getwd(), OUTDIR=getwd()){
  if(is.null(datatype)){stop("specify data type (B, N)")}
  
  # filter matrix and cls.vector
  data=data[,logicalVectors[[i]]&!is.na(logicalVectors)]
  cls.vector=cls.vector[logicalVectors[[i]]&!is.na(logicalVectors)]
  dataname=names(logicalVectors)[i]

  command=run.GSEA(data, cls.vector, GENESETS, datatype, dataname, clsname, WD, OUTDIR)
}

# function to go through
run.GSEA=function(data, cls.vector, GENESETS, datatype=NULL, dataname="", clsname="", WD=getwd(), OUTDIR=getwd(), GSEA_home="~/gsea-3.0.jar"){
  
  clsname=gsub("^_|_$", "", gsub("[[:punct:]]", "_", clsname)) # no special allowed here
  
  if(is.null(datatype)){stop("specify data type (B, N)")}
  
  # PATHS, create if needed
  if(!dir.exists(WD)){dir.create(WD, recursive = T); cat("Working directory made", sep="\n\n")}
  if(!dir.exists(OUTDIR)){dir.create(OUTDIR, recursive = T); cat("Output directory made", sep="\n\n")}
  setwd(WD)
  
  if(datatype=="N"){
    data=data[,order(as.numeric(cls.vector), decreasing=T)]
    cls.vector=cls.vector[order(as.numeric(cls.vector), decreasing=T)]
    CLUSTER_COMPARISON="continuous_phenotype"
  }
  if(datatype=="B"){
    data=data[,order(as.numeric(cls.vector), decreasing=T)]
    cls.vector=cls.vector[order(as.numeric(cls.vector), decreasing=T)]
    CLUSTER_COMPARISON=paste(unique(cls.vector), collapse="_versus_")
  }
  
  print("data sorted based on cls vector")
  #*********************** Need gct? ***************************
  
  NAME_OUT2=paste(dataname, gsub("-", "_", clsname), sep="_")
  GEXP=paste(WD, "/", NAME_OUT2, ".gct", sep="")
  
  print("making gct file")
  
  gsea.write.gct(data, GEXP)
  
  print("gct ready")
  
  #****************************** genesets ***********************************
  if(is.list(GENESETS)){
    writeGMT(GENESETS, "geneset_temp.gmt")
    print("gmt made from named list")
    GENESETS="geneset_temp.gmt"
  }
  # no need to do anything
  
  #*************************************************************************************
  
  if(datatype=="N"){
    print("this is numeric feature class, making class file")
    
    A1=paste("#numeric", sep = "\t")
    A2=paste("#feat", sep = "\t")
    A3=t(data.frame(as.character(cls.vector)))
    CLS=paste0(OUTDIR,"/",dataname, "_", CLUSTER_COMPARISON, "_", gsub("-", "_", clsname), ".cls")
    write(A1, file = CLS, append = F) # All .cls files require this dummy header line.
    write(A2, file = CLS, append = T)
    write.table(A3, file = CLS, append = T, quote = F, sep = "\t",
                na = "", row.names = F, col.names = F)
    PARAMS=" -collapse false -mode Max_probe -norm meandiv -nperm 1000 -permute sample -rnd_type no_balance -scoring_scheme weighted -metric Pearson -sort real -order descending -include_only_symbols true -make_sets true -median false -num 100 -plot_top_x 25 -rnd_seed timestamp -save_rnd_lists false -set_max 500 -set_min 5 -zip_report false -gui false"
  }
  
  if(datatype=="B"){
    print("this is binary/two group feature class")
    
    # make vector with labels
    cls.vector=cls.vector[!is.na(cls.vector), drop=F]
    
    A1=paste(length(cls.vector), 2, 1, sep = "\t")
    A2=paste("#", "1", "0", sep = "\t")
    A3=t(data.frame(as.character(cls.vector)))
    CLS=paste0(OUTDIR,"/", dataname, "_", CLUSTER_COMPARISON, "_", gsub("-", "_", clsname), "_temp.cls")
    write(A1, file = CLS, append = F) # All .cls files require this dummy header line.
    write(A2, file = CLS, append = T)
    write.table(A3, file = CLS, append = T, quote = F, sep = "\t",
                na = "", row.names = F, col.names = F)
    
    PARAMS=" -collapse false -mode Max_probe -norm meandiv -nperm 1000 -permute sample -rnd_type no_balance -scoring_scheme weighted -metric Signal2Noise -sort real -order descending -include_only_symbols true -make_sets true -median false -num 100 -plot_top_x 25 -rnd_seed timestamp -save_rnd_lists false -zip_report false -gui false"
    
  }
  
  RUN=paste0("java -Xmx5000m -cp ", GSEA_home, " xtools.gsea.Gsea -res ", GEXP,
             " -cls ", CLS, 
             " -gmx ", GENESETS,
             " -rpt_label ", NAME_OUT2, "_", CLUSTER_COMPARISON,
             " -out ", OUTDIR, " ", PARAMS)
}

gsea.write.gct <- function(exprMat, gctFn){
  nGenes = nrow(exprMat)
  nConds = ncol(exprMat)
  write("#1.2", file = gctFn, append = F) # All .gct files require this dummy header line.
  write(paste(nGenes, nConds, sep = "\t"), file = gctFn, append = T)
  write(paste("Name", "Description", paste(colnames(exprMat), collapse = "\t"), sep = "\t"), file = gctFn, append = T)
  # The second column of the .gct file, "Description", is filled out with "na"'s. 
  rownames(exprMat) = paste(rownames(exprMat), "na", sep = "\t") # Append "\tna" to every gene clsname. 
  write.table(exprMat, file = gctFn, append = T, quote = F, sep = "\t",
              na = "", row.names = T, col.names = F)
}

writeGMT <- function #Create a gmt (gene matrix transposed) file
### Createss a gmt (gene matrix transposed) file such as those
### provided by mSigDB or geneSigDB, from an R list object.
### Function by Levi Waldron.
(object,
 ### R list object that will be converted to GMT file.  Each element
 ### should contain a vector of gene names, and the names of the
 ### elements will used for the gene set names
 fname
 ### Output file name for .gmt file
){
  if (class(object) != "list") stop("object should be of class 'list'")
  if(file.exists(fname)) unlink(fname)
  for (iElement in 1:length(object)){
    write.table(t(c(make.names(rep(names(object)[iElement],2)),object[[iElement]])),
                sep="\t",quote=FALSE,
                file=fname,append=TRUE,col.names=FALSE,row.names=FALSE)
  }
  ### Called for the effect of writing a .gmt file
}


GSEA = function(gene_list, genesets, pval) {
  set.seed(54321)
  library(dplyr)
  library(gage)
  library(fgsea)
  
  if ( any( duplicated(names(gene_list)) )  ) {
    warning("Duplicates in gene names")
    gene_list = gene_list[!duplicated(names(gene_list))]
  }
  if  ( !all( order(gene_list, decreasing = TRUE) == 1:length(gene_list)) ){
    warning("Gene list not sorted")
    gene_list = sort(gene_list, decreasing = TRUE)
  }
  
  # if not list:
  if(is.character(genesets)){
    genesets = fgsea::gmtPathways(genesets)
  }
  
  fgRes <- fgsea::fgsea(pathways = genesets, 
                        stats = gene_list,
                        minSize=15,
                        maxSize=600,
                        nperm=10000) %>% 
    as.data.frame() %>% 
    dplyr::filter(padj < !!pval)
  #print(dim(fgRes))
  
  ## Filter FGSEA by using gage results. Must be significant and in same direction to keep 
  gaRes = gage::gage(gene_list, gsets=genesets, same.dir=TRUE, set.size =c(15,600))
  
  ups = as.data.frame(gaRes$greater) %>% 
    tibble::rownames_to_column("Pathway") %>% 
    dplyr::filter(!is.na(p.geomean) & q.val < pval ) %>%
    dplyr::select("Pathway")
  
  downs = as.data.frame(gaRes$less) %>% 
    tibble::rownames_to_column("Pathway") %>% 
    dplyr::filter(!is.na(p.geomean) & q.val < pval ) %>%
    dplyr::select("Pathway")
  
  #print(dim(rbind(ups,downs)))
  keepups = fgRes[fgRes$NES > 0 & !is.na(match(fgRes$pathway, ups$Pathway)), ]
  keepdowns = fgRes[fgRes$NES < 0 & !is.na(match(fgRes$pathway, downs$Pathway)), ]
  
  ### Collapse redundant pathways
  Up = fgsea::collapsePathways(keepups, pathways = genesets, stats = gene_list,  nperm = 500, pval.threshold = 0.05)
  Down = fgsea::collapsePathways(keepdowns, genesets, gene_list,  nperm = 500, pval.threshold = 0.05) 
  
  fgRes = fgRes[ !is.na(match(fgRes$pathway, 
                              c( Up$mainPathways, Down$mainPathways))), ] %>% 
    arrange(desc(NES))
  fgRes$pathway = stringr::str_replace(fgRes$pathway, "GO_" , "")
  
  fgRes$Enrichment = ifelse(fgRes$NES > 0, "Up-regulated", "Down-regulated")
  filtRes = rbind(head(fgRes, n = 10),
                  tail(fgRes, n = 10 ))
  g = ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
    geom_segment( aes(reorder(pathway, NES), xend=pathway, y=0, yend=NES)) +
    geom_point( size=5, aes( fill = Enrichment),
                shape=21, stroke=2) +
    scale_fill_manual(values = c("Down-regulated" = "dodgerblue",
                                 "Up-regulated" = "firebrick") ) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title="GSEA - Biological Processes") + 
    theme_minimal()
  
  output = list("Results" = fgRes, "Plot" = g)
  return(output)
}