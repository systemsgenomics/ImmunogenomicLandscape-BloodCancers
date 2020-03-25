# ******************************************************************************
#       FUNCTION: GSVA run, obtain eFDR | 12/02/2013 | Petri PÃ¶lÃ¶nen
# ******************************************************************************

library(GSVA)
library(parallel)
library(snow)

# ******************************************************
# Running notes: create file, named anything*EXPNAME*anything*ENTITYDF.txt
#                This file is needed to run in multiple folders
#                contains at least NAME, ENTNAME columns
#*******************************************************


#******************** Experiment, output processing, modify ***********************************

# This should mach your entity data frame partially
# EXPNAME="AgilentG4502A_07"
# EXPNAME="IlluminaHiSeq_RNASeqV2_geneExp"
# EXPNAME="illuminahiseq_rnaseqv2"
EXPNAME="data9544"

# ENT2 contains names for entities, create this file, contains NAME, ENTNAME columns, which are needed
ENT2 <- read.delim(list.files(,paste(EXPNAME,"*.*ENTITYDF.txt", sep=""), ignore.case=T)[1], header=T, stringsAsFactors=F)

# limit to certain samples
LIMITTO=F

# make fm suitable features
FM_PROCESSING=F

PANCAN=F # all in ENTITYDF are run
CORES=3 # parallel for PANCAN
TYPE="data9544" # if PANCAN is F, specify type, this can be any NAME in ENTITYDF

#******************* GSVA parameters *******************************
METHOD="gsva" # ssgsea, plage, zscore
RNASEQ=F
PARALLEL_SZ=10 # parallel run, important for several geneset testing

FDR=T # estimate significance, get table of p-values, this takes long time

# if FDR TRUE
SIMULATION="_GSVA_permutation"
NUMSIM=100
PERMUTATE_SUBJECT=F
PERMUTATE_GENE=T
PERMUTATE_GENE_LEAN=F # compute for different number of genes
BOOTSTRAP=F
PVAL=T

#******************* Pathways to run ******************************

# which pathway to run
SINGLEPATHW="Combined_pathway_drug_signatures_2017.gmt"
SINGLEPATHW="CLL_genesets"

GENESETS="~/genesets/"
GENESETS="/research/groups/sysgen/PROJECTS/HEMAP/HEMAP/dat2figs/revision_2018/"

#******************************* After this, you should not need to change anything **************
mainDir=getwd()

FUNPATHW=function(PATHW, matrix, NAME, entname){

  #**************************** load genesets *******************************************
  PATHW=gsub(".gmt", "", SINGLEPATHW)

  if(LIMITTO){
    FILE="AML_GSE13159,GSE13204_BHSNE_mean-shift.txt"
    A=read.delim(FILE, sep="\t", header=T, stringsAsFactors=F)
    FILE=gsub(".txt", "", FILE)
    LIMIT=A$ID
    matrix=matrix[,colnames(matrix)%in%LIMIT]
  }else{
    FILE="all_samples"
  }

  file=paste(GENESETS, PATHW, ".gmt", sep="")

  if(file.exists(paste0(GENESETS, PATHW, "_listA_tempfile.Rdata"))){
    load(paste0(GENESETS, PATHW, "_listA_tempfile.Rdata"))
    print("previously computed geneset list loaded")
    print("if you want to re-calculate, rm file:")
    print(paste0(GENESETS, PATHW, "_listA_tempfile.Rdata"))
  }else{

    # Geneset list
    Onc.pathways=read.delim(file, stringsAsFactors = FALSE, header=F, col.names = paste("V",1:max(count.fields(file, sep = '\t'), na.rm = T)), fill = TRUE)

    # Make list
    listA=mclapply(1:length(Onc.pathways[,1]), function(i){A=as.character(Onc.pathways[i,3:length(Onc.pathways),])
    B=A[!A==""&!A=="NA"]}, mc.cores=PARALLEL_SZ)

    names(listA) <- Onc.pathways[,1]

    save(listA, file=paste0(GENESETS, PATHW, "_listA_tempfile.Rdata"))
    names(listA) <- Onc.pathways[,1]

  }

  if(FM_PROCESSING){
    colnames(matrix)=gsub("\\.", "-", substr(colnames(matrix), 1, 15))
    names(listA) <- paste("N:SAMP:", gsub(":", "", names(listA)), "_", "GSVA", ":::::", sep="")

  }

  if(length(listA)==1){
    listA=append(listA, listA)
  }

  # remove duplicates just in case....
  listA=listA[!duplicated(names(listA))]

  A=unlist(mclapply(listA, function(pw)sum(rownames(matrix)%in%pw), mc.cores=3))

  listA=listA[!A<5&!A>500]
  print(paste0("filtering out ", sum(A<5), " genesets, size too small"))
  print(paste0("filtering out ", sum(A>500), " genesets, size too large"))

  print(paste(length(listA), "genesets after filtering"))

  #************************************** GSVA run ******************************************************

  if(file.exists(paste(NAME,  "_", FILE, "_", PATHW, "_GSVA.Rdata", sep=""))){

    # load matrix
    load(paste(NAME, "_", FILE, "_", PATHW, "_GSVA.Rdata", sep=""))
    listA=listA[names(listA)%in%rownames(gsva_es)]

    print("previously calculated GSVA matrix loaded")
    print("if you want to re-calculate, rm file:")
    print(paste(NAME, "_", FILE, "_", PATHW, "_GSVA.Rdata", sep=""))

  }else{
    #****************** Gene Set Variation Analysis function --> matrix with enrichment scores for genesets in each sample *****************

    if(BOOTSTRAP){
      if(RNASEQ){
        gsva_es <- gsva(as.matrix(matrix), method=METHOD, listA, mx.diff=F, tau=0.25, verbose=T, rnaseq=TRUE, min.sz=5, max.sz=500, no.bootstraps=NUMSIM, parallel.sz=PARALLEL_SZ)
      }else{
        gsva_es <- gsva(as.matrix(matrix), method=METHOD, listA, mx.diff=F, tau=0.25, verbose=T, rnaseq=F, min.sz=5, max.sz=500, no.bootstraps=NUMSIM, parallel.sz=PARALLEL_SZ)
      }

      # write permutations
      write.table(gsva_es$bootstrap$p.vals.sign, file=paste(NAME,"_", PATHW, "_bootstrap_pvals_GSVA.txt", sep=""), sep="\t", col.names=T, row.names=T, quote=FALSE)
      save(gsva_es, file=paste(NAME, "_", PATHW, "_bootstrap_object_GSVA.Rdata", sep=""))

      gsva_es=gsva_es$es.obs

    }else{

      if(RNASEQ){
        gsva_es <- gsva(as.matrix(matrix), listA, mx.diff=F, method=METHOD, tau=0.25, verbose=T, rnaseq=TRUE, parallel.sz=PARALLEL_SZ)$es.obs
      }else{
        gsva_es <- gsva(as.matrix(matrix), listA, mx.diff=F, method=METHOD, tau=0.25, verbose=T, rnaseq=F, parallel.sz=PARALLEL_SZ)$es.obs

      }
      listA=listA[names(listA)%in%rownames(gsva_es)]
    }


    # make fm ready feature table
    if(FM_PROCESSING){
      write.table(t(c("N:SAMP", as.character(colnames(gsva_es)))), file=paste(NAME, "_", FILE, "_", PATHW, "_GSVA.tsv", sep=""), sep="\t", col.names=F, row.names=F, quote=FALSE, append=F)
      write.table(gsva_es, file=paste(NAME, "_", FILE, "_", PATHW, "_GSVA.tsv", sep=""), sep="\t", col.names=F, row.names=T, quote=FALSE, append=T)
      write.table(gsva_es, file=paste(NAME, "_", FILE, "_", PATHW, "_GSVA.txt", sep=""), sep="\t", col.names=T, row.names=T, quote=FALSE)

    }else{
      write.table(gsva_es, file=paste(NAME, "_", FILE, "_", PATHW, "_GSVA.txt", sep=""), sep="\t", col.names=T, row.names=T, quote=FALSE)
    }

    save(gsva_es, file=paste(NAME, "_", FILE, "_", PATHW, "_GSVA.Rdata", sep=""))

    print("GSVA done: ")
    print(NAME)

  }

  #**************************************************** eFDR ********************************************************

  if(FDR){
    # eFDR my way
    if(PERMUTATE_SUBJECT){
      PERMTYPE="sampleperm"

      # Sample permutation
      for(i in 1:length(listA)){

        FUN_SIM=function(id, i){
          bootstrap.percent=0.632
          n.samples=dim(matrix)[2]
          bootstrap.nsamples <- floor(bootstrap.percent * n.samples)
          ind=sample(n.samples, bootstrap.nsamples, replace=T)
          sim_es <- gsva(matrix[,ind], listA[i], mx.diff=T, verbose=F, parallel.sz=PARALLEL_SZ)$es.obs
          rownames(sim_es)=paste("per", id, sep="")
          if(any((i/(NUMSIM/10))==(1:10))){
            print(paste("permutation ", i, "/", NUMSIM," done", sep=""))
          }
          flush.console()
          return(sim_es)
        }

        simulations=mclapply(1:NUMSIM, FUN_SIM, mc.cores=CORES, i)

        simulations_comb=do.call("rbind", simulations)

        write.table(simulations_comb, file=paste(NAME, "_", PATHW, "_", names(listA)[i], "_", SIMULATION, "_sampleperm.txt", sep=""), sep="\t", col.names=T, row.names=T, quote=FALSE)

        print(paste("permutation done for ", names(listA)[i], sep=""))
      }
    }


    if(PERMUTATE_GENE){
      PERMTYPE="geneperm"

      print("permuting gene labels:")
      print(PERMTYPE)

      # calculate for this sequence permutations, SIM_SEQUENCE2 is automaticly finding necessary permutation steps
      SIM_SEQUENCE=c(5:500)

      # Make list
      list_genesets=mclapply(SIM_SEQUENCE, function(geneset_l){perm_genesets=lapply(1:NUMSIM, function(i, ...){random.list=rownames(matrix[sample(length(matrix[,1]), geneset_l), ])})
      names(perm_genesets)=paste("per", 1:NUMSIM, "_geneset_length_", geneset_l, sep="")
      return(perm_genesets)}, mc.cores=PARALLEL_SZ)

      list_genesets=unlist(list_genesets, recursive = F)


      if(RNASEQ){
        sim_es <- gsva(as.matrix(matrix), method=METHOD, list_genesets, mx.diff=F, tau=0.25, verbose=T, rnaseq=T, parallel.sz=PARALLEL_SZ)$es.obs
      }else{
        sim_es <- gsva(as.matrix(matrix), method=METHOD, list_genesets, mx.diff=F, tau=0.25, verbose=T, parallel.sz=PARALLEL_SZ)$es.obs

      }

      save(sim_es, file=paste(NAME, "_", FILE, "_", PATHW, "_GSVA_genepermutations_lean_all.Rdata", sep=""))

      sim_df=lapply(SIM_SEQUENCE, function(geneset_l){find=paste("_geneset_length_", geneset_l, "$", sep="")
      sim_es[grep(find, rownames(sim_es)),]})

      save(sim_df, file=paste(NAME, "_", FILE, "_", PATHW, "_GSVA_genepermutations_lean.Rdata", sep=""))

    }

    if(PERMUTATE_GENE_LEAN){
      PERMTYPE="geneperm_lean"

      print("permuting gene labels:")
      print(PERMTYPE)

      # calculate for this sequence permutations, SIM_SEQUENCE2 is automaticly finding necessary permutation steps
      SIM_SEQUENCE1=c(5:20, 25, 30, 40, 50, 75, 100, 200, 300, 400, 500)

      geneset_length=unlist(lapply(listA, length))

      # round to interval
      iround <- function(x, interval){
        interval[ifelse(x < min(interval), 1, findInterval(x, interval))]
      }

      SIM_SEQUENCE2=unique(iround(geneset_length, sort(SIM_SEQUENCE1)))

      if(length(SIM_SEQUENCE2)<length(SIM_SEQUENCE1)){
        SIM_SEQUENCE=SIM_SEQUENCE2
      }else{
        SIM_SEQUENCE=SIM_SEQUENCE1
      }

      # Make list
      list_genesets=mclapply(SIM_SEQUENCE, function(geneset_l){perm_genesets=lapply(1:NUMSIM, function(i, ...){random.list=rownames(matrix[sample(length(matrix[,1]), geneset_l), ])})
      names(perm_genesets)=paste("per", 1:NUMSIM, "_geneset_length_", geneset_l, sep="")
      return(perm_genesets)}, mc.cores=PARALLEL_SZ)

      list_genesets=unlist(list_genesets, recursive = F)

      if(RNASEQ){
        sim_es <- gsva(as.matrix(matrix), method=METHOD, list_genesets, mx.diff=F, tau=0.25, verbose=T, rnaseq=T, parallel.sz=PARALLEL_SZ)$es.obs
      }else{
        sim_es <- gsva(as.matrix(matrix), method=METHOD, list_genesets, mx.diff=F, tau=0.25, verbose=T, parallel.sz=PARALLEL_SZ)$es.obs
      }


      save(sim_es, file=paste(NAME, "_", FILE, "_", PATHW, "_GSVA_genepermutations_lean.Rdata", sep=""))

      sim_df=lapply(SIM_SEQUENCE, function(geneset_l){find=paste("_geneset_length_", geneset_l, "$", sep="")
      sim_es[grep(find, rownames(sim_es)),]})

      save(sim_df, file=paste(NAME, "_", FILE, "_", PATHW, "_GSVA_genepermutations_lean.Rdata", sep=""))

    }
    if(PVAL){

      print("computing eFDR for observed vs permuted score")

      # find intersect
      find=intersect(names(listA), rownames(gsva_es))
      gsva_es=gsva_es[rownames(gsva_es)%in%find,]
      listA=listA[names(listA)%in%find]

      geneset_length=unlist(lapply(listA, length))

      # round to interval
      iround <- function(x, interval){
        interval[ifelse(x < min(interval), 1, findInterval(x, interval))]
      }

      rounded_geneset=iround(geneset_length, sort(SIM_SEQUENCE))

      # go through each geneset scores, find corresponding geneset length and compute eFDR
      FUN_SIM=function(i){
        set_length=rounded_geneset[i] # geneset length
        permutation_df=sim_df[[which(SIM_SEQUENCE%in%set_length)]] # which permuation is accessed
        observed_scores=as.numeric(gsva_es[i,]) # observed scores

        #*********** go through observed scores and compute FDR *****************
        m=unlist(lapply(1:length(observed_scores), function(j){
          observed_score=observed_scores[j]
          test=observed_score>0
          if(test){
            length(which(permutation_df[,j]>observed_score))/length(permutation_df[,1])
          }else{
            length(which(permutation_df[,j]<observed_score))/length(permutation_df[,1])
          }
        }))
        #************************************************************************

        return(signif(m,2))
      }

      m=mclapply(1:length(rounded_geneset), FUN_SIM, mc.cores=PARALLEL_SZ)

      m_comb=as.matrix(do.call("rbind", m))

      colnames(m_comb)=colnames(gsva_es)

      rownames(m_comb)=rownames(gsva_es)

      write.table(m_comb, file=paste(NAME, "_", PATHW, "_GSVA_", PERMTYPE, "_eFDR.txt", sep=""), sep="\t", col.names=T, row.names=T, quote=FALSE)
      save(m_comb, file=paste(NAME, "_", PATHW, "_GSVA_", PERMTYPE, "_eFDR.Rdata", sep=""))
      print("computing eFDR for observed vs permuted score: Done")

      #**************************** make cutoff file *****************************

      # round to interval
      pw_length=data.frame(pw_newname=names(listA), geneset_length, rounded_geneset, stringsAsFactors=F)

      sim_df2=do.call(rbind, sim_df)
      sim_df2=sim_df2[!duplicated(rownames(sim_df2)),]

      geneset_l=unique(gsub("*.*_geneset_length_", "", rownames(sim_df2)))

      vals_gsva=do.call(rbind, mclapply(geneset_l, function(geneset_l){
        df=sim_df2[gsub("*.*_geneset_length_", "", rownames(sim_df2))%in%geneset_l,]
        B=sort(as.numeric(df), decreasing=T)
        l=length(B)
        c(geneset_l, B[0.05*l], B[l-0.05*l],B[0.01*l], B[l-0.01*l], B[0.001*l], B[l-0.001*l])
      }, mc.cores=CORES))

      vals_gsva=vals_gsva[order(as.numeric(vals_gsva[,1])),]

      cutoffs=do.call(rbind, lapply(pw_length$rounded_geneset, function(v){vals_gsva[vals_gsva[,1]%in%v,]}))

      df_all=cbind(pw_length, cutoffs[,2:7])

      colnames(df_all)=c("name", "geneset_length", "geneset_length_rounded", "p_0.05_pos", "p_0.05_neg", "p_0.01_pos", "p_0.01_neg", "p_0.001_pos", "p_0.001_neg")

      df_all=df_all[df_all[,1]%in%rownames(gsva_es),]

      df_all=df_all[df_all[,2]>=5,]

      write.table(df_all, paste0(NAME, "_", PATHW, "_pw_FDR_cutoffs_GSVA.txt"), sep="\t", quote=F, col.names=T, row.names=F)

    }
  }
}

# run per type
FUN=function(i){
  entname=ENT2$ENTNAME[i]
  NAME=ENT2$NAME[i]

  # Set wd and create directory
  subDir <- NAME
  dir.create(file.path(mainDir, subDir), recursive=T, showWarnings = FALSE)

  if(grepl("RData|Rdata", entname)){
    matrix <- get(load(entname))
  }else{
    matrix=data.matrix(read.delim(entname, header=T, row.names=1, stringsAsFactors=F))
  }

  print("Data matrix loaded:")
  print(NAME)

  # outputs to GSVA folder
  outDir <- paste(NAME, "/GSVA/", sep="")
  dir.create(file.path(mainDir, outDir), showWarnings = FALSE)
  setwd(file.path(mainDir, outDir))

  log=FUNPATHW(SINGLEPATHW, matrix, NAME, entname)

}

if(PANCAN){
  log=mclapply(1:length(ENT2$NAME), FUN, mc.cores=CORES)
}else{
  log=lapply(which(ENT2$NAME%in%TYPE), FUN)
}
