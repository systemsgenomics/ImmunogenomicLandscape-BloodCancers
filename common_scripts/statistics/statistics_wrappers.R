wrapper.wilcoxtest=function(genelist, data, logicalVectors, logicalVector_normals=NULL, ALTERNATIVE="greater", adj.method="BH", CORES=1, prettynum=T){
  library(parallel)
  data=do.call(rbind, mclapply(genelist, TestGeneWilcox, data, logicalVectors, logicalVector_normals=logicalVector_normals, ALTERNATIVE, prettynum, mc.cores = CORES))
  data$adj.p=p.adjust(data$p, method = adj.method)
  
  data$adj.method=adj.method
  
  # signifcode
  data$signifCode=""
  data$signifCode[as.numeric(data$adj.p)<0.1]="."
  data$signifCode[as.numeric(data$adj.p)<0.05]="*"
  data$signifCode[as.numeric(data$adj.p)<0.01]="**"
  data$signifCode[as.numeric(data$adj.p)<0.001]="***"
  data$signifCode[as.numeric(data$adj.p)<0.0001]="****"
  
  if(prettynum){
    data$adj.p=prettyNum(signif(data$adj.p,2))
  }
  
  if(adj.method=="BH"){
    colnames(data)[colnames(data)%in%"adj.p"]="FDR"
  }
  
  data=data[,c(1:5, 7,9, 6,8)]
  
  return(data)
}

FUN_GOTHROUGH=function(i, genelist, list_cancers, logicalVector_normals, TEST_FAILS=F, HG_PVAL=0.001, MW_PVAL=0.001, ALTERNATIVE="greater"){
  name=names(list_cancers)[i]
  logicalVector=as.logical(unlist(list_cancers[[i]]))
  logicalVector_normals=logicalVector_normals[!names(logicalVector_normals)%in%name]
  
  result=FUN_HGT_MW_TEST(genelist, logicalVector,logicalVector_normals, name, TEST_FAILS=F, data = data, profile = profile,HG_PVAL = HG_PVAL, MW_PVAL=MW_PVAL, ALTERNATIVE = ALTERNATIVE, CORES = 10)
}

fisher.wrapper=function(lv.list1, lv.list2, alternative="two.sided", method="bonferroni",log10=F, prettyNumbers=F,orderdata=F, cores=1){
  
  n1=names(lv.list1)
  n2=names(lv.list2)
  
  res=do.call(rbind,mclapply(seq(lv.list1), function(i){
    do.call(rbind,lapply(seq(lv.list2), function(j){
      data.frame("featureA"=n1[i], "featureB"=n2[j],fisher.2x2(lv.list1[[i]], lv.list2[[j]]))
    }))
  }, mc.cores=cores))
  
  res$adj.p=p.adjust(res$p, method=method)
  res$signifCode=""
  res$signifCode[as.numeric(res$adj.p)<0.1]="."
  res$signifCode[as.numeric(res$adj.p)<0.05]="*"
  res$signifCode[as.numeric(res$adj.p)<0.01]="**"
  res$signifCode[as.numeric(res$adj.p)<0.001]="***"
  res$signifCode[as.numeric(res$adj.p)<0.0001]="****"
  
  if(log10){
    res$p=abs(log10(res$p))
    res$adj.p=abs(log10(res$p))
  }
  
  if(prettyNumbers){
    res$log.odds=prettyNum(signif(res$log.odds,2))
    res$p=prettyNum(signif(res$p,2))
    res$adj.p=prettyNum(signif(res$adj.p,2))
  }
  
  if(orderdata)res=res[order(res$featureA, as.numeric(res$adj.p), decreasing = F),]
  return(res)
}

fisher.matrix.pairwise=function(E, plot=F, alternative="two.sided"){
  # perform a fisher exact test for each pairs of gene and saves the oddsRatio and p-value
  res <- NULL
  for(i in seq(dim(E)[2])){
    for(j in seq(dim(E)[2])){
      if(i!=j){
        f <- fisher.test(E[,i], E[,j], alternative = alternative) 
        res <- rbind(res,cbind(geneA=colnames(E)[i],geneB=colnames(E)[j],oddsRatio=f$estimate,pvalue=f$p.value))
      }   
    }
  }
  # some formatting
  res <- as.data.frame(res)
  res$geneA <- factor(res$geneA,levels=unique(res$geneA))
  res$geneB <- factor(res$geneB,levels=unique(res$geneB))
  res$oddsRatio <- as.numeric(as.character(res$oddsRatio))
  
  # avoid inf values
  res$oddsRatio[res$oddsRatio>10]=10
  res$oddsRatio[res$oddsRatio<=0]=0.01
  
  res$pvalue <- as.numeric(as.character(res$pvalue))
  # use p.adjust to correct for multi testing using a FDR
  res <- cbind(res,fdr=p.adjust(res$pvalue,"fdr"))
  # change the FDR in labels for plotting
  res$stars <- cut(res$fdr, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))
 
  if(plot){
    # plot with ggplot 2
    require(ggplot2)
    require(cowplot) # not necessary but the plot is nicer
    p <- ggplot(res, aes(geneA, geneB)) + geom_tile(aes(fill = log2(oddsRatio+0.01)),colour = "white") + scale_fill_gradient2(midpoint=1) + geom_text(aes(label=stars), color="black", size=5)
    p
    return(p)
  }
 
  return(res)
}

fun.get.cox=function(i, logicalv, DATA, univariate=T, pretty=T, TIME, STATUS){
  logicalvector=logicalv[[i]]
  n=names(logicalv)[i]

  # remove if all are NA
  DATA=DATA[,!apply(DATA[logicalvector,,drop=F], 2, function(v)all(is.na(v)))]
  
  if(univariate){
    univ=do.call(rbind, lapply(seq(dim(DATA)[2]), function(j){
      y=Surv(TIME[logicalvector], STATUS[logicalvector])
      fit=coxph(formula = as.formula(paste("y ~ ", paste(colnames(DATA)[j], collapse="+"))), DATA[logicalvector,,drop=F])

      b=cox.zph(fit)[[1]][3] # p.value

      a=summary(fit)

      pval=a$coefficients[1,5]
      coef=a$conf.int[c(1,3,4)]
      v=data.frame(colnames(DATA)[j], t(coef), pval,"","", a$concordance[1], b<0.05, stringsAsFactors = F)
      rownames(v)=NULL
      return(v)
    }))

    # adjust P-value depending on number of genes tested:
    univ[,6]=p.adjust(univ[,5], method="BH")

  }else{
    y=Surv(TIME[logicalvector], STATUS[logicalvector])
    fit=coxph(formula = as.formula(paste("y ~ ", paste(colnames(DATA), collapse="+"))), DATA[logicalvector,])

    b=cox.zph(fit)$table[,3] # p.value
    b=b[!names(b)%in%"GLOBAL"]

    a=summary(fit)

    pval=a$coefficients[,5]
    coef=a$conf.int[,c(1,3,4)]
    univ=data.frame(rownames(a$coefficients), coef, pval, "", "", a$concordance[1], b<0.05, stringsAsFactors = F)
    rownames(univ)=NULL

    # no need to adjust P-value here:
    univ[,6]=univ[,5]
  }

  univ=univ[order(univ[,5]),]

  univ[,7][as.numeric(univ[,6 ])<0.1]="."
  univ[,7][as.numeric(univ[,6 ])<0.05]="*"
  univ[,7][as.numeric(univ[,6 ])<0.01]="**"
  univ[,7][as.numeric(univ[,6 ])<0.001]="***"
  univ[,7][as.numeric(univ[,6 ])<0.0001]="****"
  univ[,7][is.na(univ[,7])]=""

  val=data.frame(univ, n, stringsAsFactors = F)
  colnames(val)=c("Feature", "exp(coef)", "lower .95", "upper .95", "P", "Adj.P", "Signif", "concordance", "zph P < 0.05", "Name")

  if(pretty){
    for(i in c(2:6,8)){
      val[,i]=prettyNum(signif(val[,i],2))
    }
  }

  return(val)
}

fun.cox.elasticnet=function(DATA_ORG, time, status, cores=8, nfold=10,min.elnet=0.01,max.elnet=0.99, summary.km="85th_15th_percentile", percentage=0.25, REPEATS=100){
  library(glmnet)

  # run regularized regression
  a <- seq(min.elnet, max.elnet, 0.05)

  if(percentage==0){
    y=Surv(time, status)

    formdf1 <- as.formula(paste(" ~ ", paste(colnames(DATA_ORG),collapse="+")))
    x=model.matrix(formdf1, data=DATA_ORG)

    # go through alpha sequence with nfolds and repeat (100 times etc) cross validation on different sets of samples
    enet.rep=do.call(rbind, parallel::mclapply(seq(REPEATS), function(r){

      s <- do.call(rbind, lapply(a, function(i){
        cv <- cv.glmnet(x=x,y=y, family = "cox", nfold = nfold, type.measure = "deviance", alpha = i, standardize=F)
        data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.min], lambda.min = cv$lambda.min, alpha = i)
      }))

      return(s)
    }, mc.cores=cores))

    # find mean alpha and lambda using all repeats
    search_m=do.call(rbind, lapply(a, function(alpha){
      cvm=mean(enet.rep$cvm[enet.rep$alpha==alpha])
      lambda=mean(enet.rep$lambda.min[enet.rep$alpha==alpha])
      data.frame(alpha, cvm, lambda)
    }))

    # minimum cvm error:
    cv3 <- search_m[search_m$cvm == min(search_m$cvm), ]

    print(cv3$alpha)
    print(cv3$lambda)
    print(cv3$cvm)

    # fit model with tuned alpha, use lambda sequence here, not single value
    md3 <- glmnet(x=x,y=y, family = "cox", alpha = cv3$alpha, standardize=F)

    coefficients <- coef(md3, s = cv3$lambda)[-1,,drop=F]
    coefficients=coefficients[order(abs(coefficients[,1]), decreasing = T),,drop=F]
    active_coefficients <- coefficients[,1] != 0

    coefficients_all=coefficients[active_coefficients,,drop=F]

    comprisk=DATA_ORG[,match(rownames(coefficients_all), colnames(DATA_ORG))]
    risk_patient=as.numeric(as.numeric(coefficients_all) %*% data.matrix(t(comprisk)))

    df=data.frame(x=risk_patient)
    ggplot(df, aes(x=x))+
      geom_density(color="lightblue", fill="lightblue")

    fun.kapplanMeier(time, status, CONTINUOUS=risk_patient, CONTINUOUS_SUMMARY = summary.km, MONTHS=F, PVAL=1, INDIVIDUAL_GROUPS=F, NAME = "predicted risk")

    a=list(summary(coxph(y ~ PI.train, data.frame("PI.train"=risk_patient))))
    c=list(coefficients_all)

    out=c(a,c)
    names(out)=c("PI.test", "coefficients")
    return(out)
  }else{
    print("Splitting to test and train data")
    # divide data so that there are equal proportions of outcomes with training and test set!
    nr.samples=floor(percentage*dim(DATA_ORG)[1])
    nr.events=floor(table(status)*percentage)

    set.seed(0)
    x <- which(status==1)
    vf <- logical(length(status==1))
    vf[x[sample.int(length(x), nr.events[2])]] <- TRUE

    set.seed(0)
    x <- which(status==0)
    vf2 <- logical(length(status==0))
    vf2[x[sample.int(length(x), nr.events[1])]] <- TRUE

    find=vf|vf2

    # define training and testing datasets
    TRAIN=DATA_ORG[!find,]
    TEST=DATA_ORG[find,]

    ytrain=Surv(time[!find], status[!find])
    ytest=Surv(time[find], status[find])

    # this must look the same:
    log=fun.kapplanMeier(time, status, GROUPS=(rownames(DATA_ORG)%in%rownames(TRAIN)*1), MONTHS=F, PVAL=1, INDIVIDUAL_GROUPS=T, NAME = "These must be similar\nand high pval: 1 TRAIN, 0 TEST")

    formdf1 <- as.formula(paste(" ~ ", paste(colnames(TRAIN),collapse="+")))
    x=model.matrix(formdf1, data=TRAIN)

    # go through alpha sequence with nfolds and repeat (100 times etc) cross validation on different sets of samples
    enet.rep=do.call(rbind, parallel::mclapply(seq(REPEATS), function(r){

      s <- do.call(rbind, lapply(a, function(i){
        cv <- cv.glmnet(x=x,y=ytrain, family = "cox", nfold = floor((dim(DATA_ORG)[1]/10)), type.measure = "deviance", alpha = i, standardize=F)
        data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.min], lambda.min = cv$lambda.min, alpha = i)
      }))

      return(s)
    }, mc.cores=cores))

    # find mean alpha and lambda using all repeats
    search_m=do.call(rbind, lapply(a, function(alpha){
      cvm=mean(enet.rep$cvm[enet.rep$alpha==alpha])
      lambda=mean(enet.rep$lambda.min[enet.rep$alpha==alpha])
      data.frame(alpha, cvm, lambda)
    }))

    # minimum cvm error:
    cv3 <- search_m[search_m$cvm == min(search_m$cvm), ]

    print(cv3$alpha)
    print(cv3$lambda)
    print(cv3$cvm)

    # fit model with tuned alpha, use lambda sequence here, not single value
    md3 <- glmnet(x=x,y=ytrain, family = "cox", alpha = cv3$alpha, standardize=F)

    coefficients <- coef(md3, s = cv3$lambda)[-1,,drop=F]
    coefficients=coefficients[order(abs(coefficients[,1]), decreasing = T),,drop=F]
    active_coefficients <- coefficients[,1] != 0

    coefficients_all=coefficients[active_coefficients,,drop=F]
    comprisk=TRAIN[,match(rownames(coefficients_all), colnames(TRAIN))]
    risk_patient=as.numeric(as.numeric(coefficients_all) %*% data.matrix(t(comprisk)))

    comprisk=TEST[,match(rownames(coefficients_all), colnames(TEST))]
    risk_patient_test=as.numeric(as.numeric(coefficients_all) %*% data.matrix(t(comprisk)))

    df=data.frame(x=risk_patient)
    ggplot(df, aes(x=x))+
      geom_density(color="lightblue", fill="lightblue")

    df=data.frame(x=risk_patient_test)
    ggplot(df, aes(x=x))+
      geom_density(color="indianred", fill="indianred")

    fun.kapplanMeier(time[!find], status[!find], CONTINUOUS=risk_patient, CONTINUOUS_SUMMARY = summary.km, MONTHS=F, PVAL=1, INDIVIDUAL_GROUPS=F, NAME = "predicted risk")
    fun.kapplanMeier(time[find], status[find], CONTINUOUS=risk_patient_test, CONTINUOUS_SUMMARY = summary.km, MONTHS=F, PVAL=1, INDIVIDUAL_GROUPS=F, NAME = "validated risk")

    a=list(summary(coxph(ytrain ~ PI.train, data.frame("PI.train"=risk_patient))))
    b=list(summary(coxph(ytest ~ PI.test, data.frame("PI.test"=risk_patient_test))))
    c=list(coefficients_all)

    out=c(a,b,c)
    names(out)=c("PI.train", "PI.test", "coefficients")
    return(out)
  }
}
