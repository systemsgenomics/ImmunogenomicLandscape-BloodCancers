mixtureM.gene=function(gene, d.matrix, FIX=F){
  library(mclust)
  library(parallel)
  dat=d.matrix[rownames(d.matrix)%in%gene,,drop=F]
  model=Mclust(as.numeric(dat),2)
  
  # first test if we can find clear components:
  class=model$classification
  
  class[class == 1] = -1 # background
  class[class == 2] = 1 # expressed high
  
  for (j in 1:dim(dat)[2]) { # Adjust the classes.
    if (dat[,j] < 4) {
      class[j] = -1
      next
    }
    if (model$uncertainty[j] > 0.1) {
      class[j] = 0
      next
    }
    if (dat[,j] > 10) {
      class[j] = 1
    }
  }
  
  # in some cases background distribution does not exist (gene always expressed, -1 missing)
  # or over 90% of the samples are uncertain
  # in these cases use original model component for expressed - notexpressed
  # These cases have only one distribution, try to identify if it is background distribution or highly expressed gene distribution
  if(FIX){
    
    # this means, that only one distribution was found, not two:
    if((sum(class==-1)==0|sum(class==0)>0.9*length(class))&!all(class==1)){
      
      # odds, how many samples have expression above/below 6
      odds=sum(dat[class==0]>=6)/sum(class==0)
      odds2=sum(dat[class==0]<6)/sum(class==0)
      
      # 60% of the uncertain have expression above 6
      if(odds>0.6){
        class[model$classification == 2] = 1
      }
      
      # 60% of the uncertain have expression below 6
      if(odds2>0.6){
        class[model$classification == 1] = -1
      }
    }
  }
  
  ret_class=data.matrix(t(class))
  dimnames(ret_class)=dimnames(dat)
  return(ret_class)
}

# data=t(get(load("/research/groups/sysgen/PROJECTS/HEMAP/HEMAP/dat2figs/data/data9544_with_gene_symbols.RData")))
# d.matrix=data
# profile=do.call(rbind, mclapply(rownames(data), mixtureM.gene, data, FIX=T, mc.cores=8))
# save(profile, file="mixtureM_profile.Rdata")
