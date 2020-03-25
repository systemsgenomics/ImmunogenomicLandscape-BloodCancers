
plot.boxplot=function(gene, logicalVectors, names.lv=NULL, clean.name=c("Cancer_", ":::::", "DUFVA_", ".*.:"), data=NULL, matrix=NULL, color.v=NULL, order.bl=F, y.range=NULL, spread=F, sample.annotation=NULL, sample.annotation.color="yellow", N=F,sample.color.continuous=NULL,sample.color.continuous.colors=c("#a52a2a", "#d4d3d1", "#fea719"), sample.color.continuous.breaks=c(0, 0.5, 1), sample.color.continuous.limits=c(0,1), sample.color.continuous.title="CIITA methylation\nbeta value", intercept.y=NULL, intercept.y.2=NULL, ylab="Expression (log2)", outlier.size=0.5) {
  library(ggplot2)
  
  if(is.null(matrix)&is.null(data))stop("No data to plot, check data/matrix")
  
  if(is.null(names.lv))names.lv=names(logicalVectors)
  
  if(!all(clean.name=="")|!is.null(clean.name)){
    
    clean.name=paste(clean.name, collapse="|")
    
    title=gsub("N:....:|:::::|DUFVA_|.*.:", "", gene)
    names.lv=gsub("Cancer_", " ", names.lv)
    
    names.lv=gsub(eval(clean.name), "",  names.lv)
    title=gsub(eval(clean.name), "", title)
    
    names.lv=gsub("_", " ",  names.lv)
    title=gsub("_", " ", title)
  }
  
  if(is.null(color.v)){
    # http://tools.medialab.sciences-po.fr/iwanthue/
    color.v=c("#d7a85b",
              "#4d56b9",
              "#acb839",
              "#5e2883",
              "#42c87f",
              "#bf4da5",
              "#75b550",
              "#8479e6",
              "#cea632",
              "#5488e3",
              "#d38725",
              "#3e397f",
              "#a4a94e",
              "#be7cde",
              "#4d7122",
              "#8460b5",
              "#62ac6a",
              "#86275d",
              "#43c8ac",
              "#cf483f",
              "#748ed7",
              "#ca6e37",
              "#de91dc",
              "#926a26",
              "#94589d",
              "#822c17",
              "#d76092",
              "#d2745b",
              "#b24258",
              "#d35760")[seq(logicalVectors)]
  }
  
  
  if(!is.null(matrix)){
    gene2=ifelse(grepl("GEXP", gene), gene, paste0("'N:GEXP:", gene, ":::::'"))
    D=as.numeric(read.delim(pipe(paste0("grep -Fw ", gene2, " ", matrix)), row.names = 1, header=F))
    names=scan(pipe(paste0("head -n1 ", matrix)), "n")
  }
  
  if(!is.null(data)){
    if("TP53"%in%colnames(data)){
      D = as.numeric(data[,colnames(data)%in%gene])
      names=colnames(data)
    }else{
      D = as.numeric(data[rownames(data)%in%gene,])
      names=colnames(data)
    }
    
  }
  
  bplot_list=lapply(logicalVectors, function(v){
    D[v]
  })
  
  if(order.bl){
    ord=sapply(bplot_list, median)
    color.v=color.v[order(ord, decreasing = T)]
    bplot_list=bplot_list[order(ord, decreasing = T)]
  }
  
  if(!is.null(sample.annotation)){
    sample.annotation=lapply(logicalVectors, function(v){
      D[v]%in%D[names%in%sample.annotation&v]
    })
    names(sample.annotation)=gsub("_", " ", names.lv)
    if(order.bl)sample.annotation=sample.annotation[order(ord, decreasing = T)]
  }
  
  if(!is.null(sample.color.continuous)){
    sample.color.continuous=lapply(logicalVectors, function(v){
      sample.color.continuous[v]
    })
    
    if(order.bl)sample.color.continuous=sample.color.continuous[order(ord, decreasing = T)]
  }
  
  if(N){
    names(bplot_list)=paste0(names(bplot_list), " (n=", sapply(logicalVectors, sum, na.rm=T), ")")
  }
  
  plots=plot.boxplot.list(bplot_list, gene, color.v, order.bl, y.range, spread, sample.annotation, sample.annotation.color, sample.color.continuous,sample.color.continuous.colors,sample.color.continuous.breaks,sample.color.continuous.limits,sample.color.continuous.title, intercept.y, intercept.y.2, ylab, outlier.size)
  return(plots)
}

plot.boxplot.list=function(bplot_list, title="boxplot", color.v, order.bl=F, y.range=NULL, spread=F, sample.annotation=NULL, sample.annotation.color="red", sample.color.continuous=NULL,sample.color.continuous.colors=c("#a52a2a", "#d4d3d1", "#fea719"), sample.color.continuous.breaks=c(0, 0.5, 1), sample.color.continuous.limits=c(0,1), sample.color.continuous.title="methylation\nbeta value", intercept.y=NULL, intercept.y.2=NULL, ylab="Expression (log2)", outlier.size=0.5){
  
  
  if(!is.null(sample.color.continuous)){
    df=data.frame(reshape2::melt(bplot_list), "cont.colorv"=unlist(sample.color.continuous))
  }else{
    df=reshape2::melt(bplot_list)
  }
  
  df$class <- factor(df[,2], levels = unique(as.character(df[,2])),ordered = TRUE)
  
  df$Expression=as.numeric(as.vector(df[,1]))
  
  p <- ggplot(data=df, aes(x=class, y=Expression))+stat_boxplot(geom = "errorbar", width = 0.5)
  
  
  if(spread){
    
    names(color.v)=levels(df$class)
    colScale=scale_colour_manual(name="class", values=color.v)
    df$spreadCol=color.v[match(df$class, names(color.v))]
    
    if(!is.null(sample.annotation)){
      df$spreadCol2=ifelse(melt(sample.annotation)[,1], sample.annotation.color, "grey20")
      
      p=p+geom_boxplot(size = 0.0001, outlier.colour = NA, outlier.size = 0, fill=unique(df$spreadCol), fatten = 10000)+
        geom_jitter(width = .3, size=outlier.size, shape=21, color="grey20", fill=df$spreadCol2, stroke = 0)
    }
    
    if(is.null(sample.annotation)&!is.null(sample.color.continuous)){
      library(ggnewscale)
      
      p=p+geom_boxplot(size = 0.0001, outlier.size = 0,outlier.colour = NA,outlier.shape = NA, fill=unique(df$spreadCol), fatten = 10000)+
        scale_fill_manual(values = unique(df$spreadCol), na.value = "#d4d3d1") + # scale for box fill (cancer color)
        new_scale_fill() +
        geom_jitter(aes(fill = cont.colorv), width = .3, size=outlier.size, shape=21, stroke = 0) +
        scale_fill_gradientn(colors = sample.color.continuous.colors, breaks = sample.color.continuous.breaks, limits = sample.color.continuous.limits, guide = "colourbar") + # scale for dot fill (methylation beta value)
        guides(fill = guide_colorbar(title = sample.color.continuous.title))
      
    }
    
    if(is.null(sample.annotation)&is.null(sample.color.continuous)){
      p=p+geom_boxplot(size = 0.0001, outlier.colour = NA, outlier.size = 0, fill=unique(df$spreadCol), fatten = 10000)+
        geom_jitter(width = .3, size=outlier.size, shape=21, color="grey20", fill=df$spreadCol)
    }
  }else{
    p=p+geom_boxplot(size = 0.0001, outlier.colour = "grey20", outlier.size = 1, fill=color.v, fatten = 10000)
  }
  
  if(!is.null(intercept.y)) p=p+geom_hline(yintercept=intercept.y, linetype="dashed", color = "red")
  if(!is.null(intercept.y.2)) p=p+geom_hline(yintercept=intercept.y.2, linetype="dashed", color = "blue")
  
  p2 <- p +
    
    #theme with white background
    theme_bw() +
    
    # titles
    theme(plot.title = element_text(color="black", size=10, hjust=0)) +
    theme(axis.title = element_text(color="black", face=NULL, size=8,angle = 90)) +
    theme(axis.title.y = element_text(size = 8, angle = 90, color="black", face=NULL)) +
    
    ylab(ylab) +
    xlab("") +
    labs(title=title) +
    
    #eliminates background, gridlines, and chart border
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    theme(panel.border= element_blank())+
    theme(plot.margin = unit(c(0.1,0.1,0.1,5), "cm"))+
    
    #draws x and y axis line
    theme(axis.line = element_line(),
          axis.line.x = element_line(color="black", size = 0.5),
          axis.line.y = element_line(color="black", size = 0.5)) +
    
    # X - axis text
    theme(axis.text.x = element_text(angle=45, hjust=1, color="black", size = 8, face=NULL),
          axis.text.y = element_text(hjust=1, color="black", size = 8, face=NULL))+
    
    # if want to limit to range
    if(!is.null(y.range))scale_y_continuous(breaks=seq(y.range[1], y.range[2], (y.range[2]-y.range[1])/5), limits = y.range)
  
  return(p2)
}

get.logical=function(annovector, a=NULL, filterv=NULL, PREFIX=NULL){
  if(sum(filterv)<1&!is.null(filterv))stop("Logical filtering vector is empty")
  if(class(annovector)!="list")stop("Annotation vector should be a list")
  
  binaryfeatures=unlist(lapply(annovector, function(annov, a, filterv, PREFIX){
    
    if(is.null(filterv))filterv=rep(T, length(annov))
    
    if(is.null(a)){
      a=unique(annov[filterv])
      a=a[!(is.na(a)|a=="na")]
      annov[!filterv]=NA
    }else{
      b=unique(annov[filterv])
      a=a[a%in%b]
    }
    FUN_BINARYFEAT=function(type, annovector){
      logv=annovector%in%type
    }
    # make binary features
    binaryfeats=lapply(a, FUN_BINARYFEAT, annov)
    if(!is.null(PREFIX))names(binaryfeats)=paste0(a ,"_", PREFIX)
    if(is.null(PREFIX))names(binaryfeats)=a
    return(binaryfeats)
  }, a, filterv, PREFIX), recursive = F)
}

logicalList2vector=function(lv.list){
  classes=names(lv.list)
  
  vector.classes=rep(NA, length(lv.list[[1]]))
  
  for(i in seq(classes)){
    vector.classes[lv.list[[i]]]=classes[i]
  }
  
  return(vector.classes)
  
}

get.data.type=function(type, feats, fm.f, ord, FDR=F, clean_name=T){
  pl=fm.f[match(feats[grepl(type, feats)], rownames(fm.f)),ord,drop=F]
  if(dim(pl)[1]==0)return(NULL)
  
  n=do.call(rbind, strsplit(rownames(pl), ":"))[,1:3]
  
  if(is.null(dim(n)))n=data.frame(t(n))
  
  rownames(pl)=apply(n, 1, paste, collapse = ":")
  
  if(clean_name){
    rownames(pl)=make.unique(gsub(".*.:", "", rownames(pl)))
    rownames(pl)=gsub("_", " ", rownames(pl))
    
  }
  
  return(pl)
}

make.hm=function(i, dat, param, cluster_rows=F, cluster_columns=F, split=NULL, text_annot=NULL, clean_name=T, WIDTH=50, use_raster=T, show_column_names=F, ...){
  pl=data.matrix(dat[[i]])
  par=param[i,]
  name=par$name
  color=gsub(" ", "", as.character(unlist(strsplit(par$colors, ","))))
  values=as.numeric(unlist(strsplit(par$values, ",")))
  type=ifelse(grepl("^B:", par[1]), "discrete", "numeric")
  labels_legend=values
  
  if(grepl("^B:", par[1]))labels_legend=c(0,1)
  if(par$scale){
    pl=t(scale(t(pl)))
    pl[pl>2]=2
    pl[pl<(-2)]=-2
    labels_legend=c(-2,-1,0,1,2)
    name=paste(name, "Z-score", sep="\n")
  }
  
  ha.r=NULL
  if(!is.null(text_annot)&!cluster_rows){
    
    if(clean_name){
      names(text_annot)=make.unique(gsub(".*.:", "",  names(text_annot)))
      names(text_annot)=gsub("_", " ",  names(text_annot))
    }
    
    text=text_annot[match(rownames(pl), names(text_annot))]
    if(any(!is.na(text))) ha.r = rowAnnotation(value = anno_text(text, gp = gpar(fontsize = 8)))
  }
  
  p=Heatmap(pl, cluster_columns = cluster_columns, cluster_rows = cluster_rows, row_names_side = "right", 
            column_names_gp = gpar(fontsize = 8), 
            row_names_gp = gpar(fontsize = 8),
            column_title_gp = gpar(fontsize = 9),
            name=par$type, col = colorRamp2(values, color), 
            show_column_names = show_column_names, width = unit(WIDTH, "mm"),
            # rect_gp = gpar(col = "white", lty = 1, lwd = 0.5),
            height = unit(3*dim(pl)[1], "mm"),
            left_annotation = ha.r,
            use_raster = use_raster,
            heatmap_legend_param = list(legend_direction = "horizontal",title=name,
                                        legend_width = unit(20, "mm"), title_position = "topcenter",color_bar = type, at = labels_legend, title_gp = gpar(fontsize = 8), legend_gp = gpar(fontsize = 8)),
            column_split = split)
  
  
  return(p)
}


plot.complexHM.fm=function(feats, fm.f, annotdf=NULL, order_columns=NULL, order_rows=NULL, clean_name=T, feats.barplot=NULL,text_annot=NULL, split.columns=T, plotting.param="/research/work/ppolonen/plotting_tools/plotting_param_fm.txt", NAME=NULL, barplot.height=3, WIDTH=50, show_column_names=F, use_raster=F){
  
  if(is.null(order_columns)){
    cluster_columns=T
    ord=colnames(fm.f)
  }else{
    cluster_columns=F
    ord=order_columns
  } 
  if(is.null(order_rows)){
    cluster_rows=T
    ord_rows=feats
  }else{
    cluster_rows=F
    ord_rows=order_rows
  }
  
  if("adj.p"%in%colnames(feats)){
    FDR=feats$adj.p
    feats=unique(feats$featureB)
  }else{
    feats=as.character(feats)
  }
  
  feats=feats[feats%in%rownames(fm.f)]
  
  param=data.table::fread(plotting.param, data.table = F)
  datatypes=param$type
  
  dat=lapply(datatypes, get.data.type, feats, fm.f, ord, FDR=F)
  names(dat)=datatypes
  dat=dat[!sapply(dat, is.null)]
  
  # remove NA columns if they exist, otherwise clustering fails:
  if(cluster_columns){
    rm=lapply(dat, function(d){colSums(is.na(d))==dim(d)[1]})
    
    filt=Reduce("|", rm)
    
    dat=lapply(dat, function(d)d[,!filt])
    
    if(sum(filt)>0)warning(paste("Removed ", sum(filt), "samples, all samples NA for clustering"))
  }
  
  if(split.columns&!is.null(annotdf)){
    split=as.character(annotdf[ord,1])
    split=factor(split, levels=unique(split))
  }else{
    split=NULL
  }
  
  param=param[match(names(dat), param[,1]),]
  
  panels=lapply(seq(dat), make.hm, dat, param, cluster_rows, cluster_columns, split, text_annot, clean_name, WIDTH,use_raster, show_column_names)
  panels=panels[!is.null(panels)]
  
  if(!is.null(annotdf)){
    annotdf=annotdf[ord,,drop=F]
    h=dim(annotdf)[2]*3
    
    if(!is.null(feats.barplot)){
      f=annotdf[,feats.barplot,drop=F]
      annotdf=annotdf[,!colnames(annotdf)%in%feats.barplot, drop=F]
      
      if(dim(annotdf)[2]==0)annotdf=data.frame("dummy"=c(rep(0, dim(f)[1])[-1], 1))
      h=dim(annotdf)[2]*3+barplot.height
      
      l.barplot=paste("HeatmapAnnotation(", paste(paste0("'", feats.barplot, "' = anno_barplot(as.numeric(f[,",seq(dim(f)[2]),"]), axis = T, gp = gpar(fill = 'grey25'), height=unit(", barplot.height, ", 'mm'))"), collapse=", "),", height = unit(dim(f)[2]*", h,", 'mm'), annotation_name_gp = gpar(fontsize = 8))")
      
      ha = eval(parse(text = l.barplot))
    }else{
      ha=HeatmapAnnotation(df = annotdf)
    }
    
    add=Heatmap(t(annotdf), cluster_columns = cluster_columns, top_annotation = ha, cluster_rows = cluster_rows, row_names_side = "right", 
                column_names_gp = gpar(fontsize = 8), 
                row_names_gp = gpar(fontsize = 8), 
                column_title_gp = gpar(fontsize = 9),
                show_column_names = show_column_names, width = unit(WIDTH, "mm"), 
                height = unit(h+10, "mm"),
                use_raster=use_raster,
                heatmap_legend_param = list(legend_direction = "horizontal",title="Annotations",
                                            legend_width = unit(20, "mm"), title_position = "topcenter", title_gp = gpar(fontsize = 8), legend_gp = gpar(fontsize = 8)), column_split = split)
    panels=append(list("Annotations"=add), panels)
    dat=append(list("Annotations"=t(annotdf)), dat)
  }
  
  heights=colSums(sapply(panels, component_height))+60
  nrows=sum(heights)
  rows_needed=heights
  ind=findInterval(seq(nrows), c(1,cumsum(rows_needed)), rightmost.closed = T, left.open=T)
  
  figure <-multipanelfigure:: multi_panel_figure(width = 150+WIDTH, height = 40+sum(heights), rows = nrows, columns = 1, panel_label_type = "lower-alpha", unit = "mm", row_spacing = unit(5, "mm"))
  
  for(i in seq(panels))figure <- multipanelfigure::fill_panel(figure,  panels[[i]], row = which(ind==i), column = 1)
  
  if(!is.null(NAME))multipanelfigure::save_multi_panel_figure(figure, filename=paste0(NAME, "_fm_complexheatmap.pdf"), limitsize = FALSE)
  return(figure)
}


plot.scatter=function(x, y,namev, lv=NULL, group=NULL, colorv=NULL, main="XY-plot", axis.plot=F, xlab="x", ylab="y", SIZE=3, add.proportions=F,add.density=F, add.legend=T,add.labels=F, rasterize=T, width=77, height=74, text.size=10){
  library(ggplot2)
  library(reshape2)
  library(RColorBrewer)
  library(ggrepel)
  library(ggrastr)
  library(cowplot)
  library(dplyr)
  
  if(is.null(colorv)){
    # http://tools.medialab.sciences-po.fr/iwanthue/
    colorv=c( "#acb839",
              "#42c87f",
              "#bf4da5",
              "#8479e6",
              "#cea632",
              "#5488e3",
              "#d76092",
              "#d38725",
              "#3e397f",
              "#a4a94e",
              "#be7cde",
              "#4d7122",
              "#8460b5",
              "#62ac6a",
              "#d7a85b",
              "#86275d",
              "#43c8ac",
              "#cf483f",
              "#748ed7",
              "#ca6e37",
              "#de91dc",
              "#926a26",
              "#94589d",
              "#5e2883",
              "#822c17",
              "#4d56b9",
              "#d2745b",
              "#b24258",
              "#75b550",
              "#d35760")
    # image(seq(colorv), 1, as.matrix(1:length(colorv)), col=colorv,
    #       xlab="", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
  }
  
  dat=data.frame("name"=namev, "x"=x,"y"=y)
  n=paste0("N = ", dim(dat)[1])
  
  if(is.null(lv))lv=rep(F, length(x))
  dat$Significant=lv
  
  if(is.null(group)){
    dat$group=""
    colorv=colorv[1]
  }else{
    dat$group=factor(group, levels=unique(group))
    dat=dat[order(dat$group),]
  }
  
  myColors <- colorv
  
  names(myColors) <- levels(dat$group)
  colScale <- scale_colour_manual(name = "group",values = myColors)
  
  if(add.legend){
    legend <- cowplot::get_legend(ggplot(dat, aes(x, y, colour = group)) + theme_bw() +
                                    geom_point(size=0.6) + colScale +theme(legend.position = "bottom", legend.direction = "horizontal", legend.title = element_blank(), 
                                                                           legend.key=element_blank(), legend.box.background=element_blank(),
                                                                           legend.text=element_text(size = text.size, face = "bold", family="Helvetica"))+
                                    theme(plot.margin=unit(c(1,1,1,1), "mm"))+
                                    guides(colour = guide_legend(override.aes = list(size=text.size/3))))
    nrow=4
  }else{
    nrow=3
  }
  
  if(rasterize){
    p=ggplot(dat, aes(x, y, colour = group)) +
      geom_point_rast(size=SIZE, raster.width = 3, raster.height = 3) + colScale +
      geom_point_rast(data = subset(dat, dat$Significant), size=SIZE*1.5, raster.width = 3, raster.height = 3)
  }else{
    p=ggplot(dat, aes(x, y, colour = group)) +
      geom_point(size=SIZE) + colScale +
      geom_point(data = subset(dat, dat$Significant), size=SIZE*1.5)
  }
  
  p=p +
    theme_classic(base_size = 12) +
    ggtitle(paste(main, n, sep="\n")) +
    labs(x = xlab, y = ylab, fill = "") +
    scale_y_continuous(breaks = pretty(dat$y, n = 6))+
    theme(legend.position = "none", legend.direction = "horizontal",
          axis.line = element_line(size=1, colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          # panel.border = element_blank(), panel.background = element_blank(),
          plot.title = element_text(size = text.size, face = "bold", family="Helvetica"),
          text=element_text()
    )+
    theme(
      axis.text.x = element_text(colour="grey20",size=text.size,face="plain", family="Helvetica"),
      axis.text.y = element_text(colour="grey20",size=text.size,face="plain", family="Helvetica"),
      axis.title.x = element_text(colour="grey20",size=text.size,face="plain", family="Helvetica"),
      axis.title.y = element_text(colour="grey20",size=text.size,face="plain", family="Helvetica")
    ) +
    theme(plot.margin=unit(c(1,1,1,1), "mm"))
  
  if(!axis.plot){
    p=p+theme(axis.title = element_blank(), # remove if needed
    axis.ticks = element_blank(), # remove if needed
    axis.line = element_blank(),
    axis.text = element_blank())
  }
  
  # add lv labels (useful if plotting genes):
  if(add.labels){
    p=p+geom_text_repel(
      data = subset(dat, dat$Significant),
      aes(label = as.character(dat$name[dat$Significant])),
      size = SIZE*1.25,
      color="black",
      box.padding = unit(0.2, "lines"),
      point.padding = unit(0.2, "lines"),
      family="Helvetica"
    )
  }
  
  theme0 <- function(...) theme( legend.position = "none",
                                 panel.background = element_blank(),
                                 panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank(),
                                 panel.spacing = unit(0,"null"),
                                 axis.ticks = element_blank(),
                                 axis.text.x = element_blank(),
                                 axis.text.y = element_blank(),
                                 axis.title.x = element_blank(),
                                 axis.title.y = element_blank(),
                                 axis.ticks.length = unit(0,"null"),
                                 panel.border=element_rect(color=NA),...)
  
  # 74mm * 99mm per final panel
  figure <-multipanelfigure:: multi_panel_figure(width = width, height = height, rows = nrow, columns = 5, panel_label_type = "none", unit = "mm", row_spacing = unit(2, "mm"), column_spacing = unit(0, "mm"))
  figure <- multipanelfigure::fill_panel(figure,  p, row = 1:3, column = 2:5)
  
  if(add.density){
    
    library(ggplot2)
    library(gridExtra)
    
    DF=dat
    
    p2 <- ggplot(DF,aes(x=x,colour=factor(group))) + colScale +
      geom_density(alpha=0.5) + 
      scale_x_continuous(breaks=NULL,expand=c(0.02,0)) +
      scale_y_continuous(breaks=NULL,expand=c(0.02,0)) +
      theme_bw() +
      theme0(plot.margin = unit(c(1,0,0,2.2),"lines")) 
    
    p3 <- ggplot(DF,aes(x=y,colour=factor(group))) + colScale +
      geom_density(alpha=0.5) + 
      coord_flip()  + 
      scale_x_continuous(labels = NULL,breaks=NULL,expand=c(0.02,0)) +
      scale_y_continuous(labels = NULL,breaks=NULL,expand=c(0.02,0)) +
      theme_bw() +
      theme0(plot.margin = unit(c(0,1,1.2,0),"lines"))
    
    p=grid.arrange(arrangeGrob(p2,ncol=2,widths=c(3,1)),
                   arrangeGrob(p,p3,ncol=2,widths=c(3,1)),
                   heights=c(1,3))
    
  }
  
  
  if(add.proportions){
    prop.dat=data.frame(main,as.data.frame(table(group)/length(group), stringsAsFactors=F), stringsAsFactors = F)
    ord=order(prop.dat$Freq, decreasing = F)
    prop.dat=prop.dat[ord,]
    prop.dat$group=factor(prop.dat$group, levels=prop.dat$group)
    myColors=colorv
    names(myColors) <- levels(dat$group)
    colScale2 <- scale_fill_manual(values = myColors)
    p2=ggplot(data=prop.dat,
              aes(x=main, y=Freq, fill=group)) + colScale2 +
      geom_bar(position="fill", stat = "identity", inherit.aes = T) +
      ggtitle(paste(" ", " ", sep="\n")) +
      theme(legend.position = "none", legend.direction = "horizontal",
            axis.line.y = element_line(size=1, colour = "black"),
            axis.line.x = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_text(colour="grey20",size=text.size,face="plain", family="Helvetica"),  
            axis.ticks.length = unit(2,"mm"),
            panel.border = element_blank(), panel.background = element_blank(),
            plot.title = element_text(size = text.size, face = "bold", family="Helvetica"),
            text=element_text())+
      theme(plot.margin=unit(c(0,0,0,0), "mm"))
    
    
    # 74mm * 99mm per panel
    figure <- multipanelfigure::fill_panel(figure,  p2, row = 1:3, column = 1)
    
  }  
  
  if(add.legend)figure <- multipanelfigure::fill_panel(figure,  legend, row = 4, column = 2:4)
  
  
  return(figure)
}


get.only.legend=function(colorv, group, position="right", direction="vertical", text.size=10){
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
                                                                         legend.text=element_text(size = text.size, face = "bold", family="Helvetica"))+
                                  theme(plot.margin=unit(c(1,1,1,1), "mm"))+
                                  guides(colour = guide_legend(override.aes = list(size=text.size/3))))
  return(legend)
}




plot.DotPlot=function (object, assay = NULL, features, cols = c("lightgrey", "blue"), col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6, group.by = NULL, split.by = NULL, scale.by = "radius", scale.min = NA, scale.max = NA){
  library(cowplot)
  
  # if seurat input:
  library(Seurat)
  
  assay <- DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  scale.func <- switch(EXPR = scale.by, size = scale_size, 
                       radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
  data.features <- FetchData(object = object, vars = features)
  data.features$id <- if (is.null(x = group.by)) {
    Idents(object = object)
  }else {
    object[[group.by, drop = TRUE]]
  }
  if (!is.factor(x = data.features$id)) {
    data.features$id <- factor(x = data.features$id)
  }
  id.levels <- levels(x = data.features$id)
  data.features$id <- as.vector(x = data.features$id)
  if (!is.null(x = split.by)) {
    splits <- object[[split.by, drop = TRUE]]
    if (length(x = unique(x = splits)) > length(x = cols)) {
      stop("Not enought colors for the number of groups")
    }
    cols <- cols[1:length(x = unique(x = splits))]
    names(x = cols) <- unique(x = splits)
    data.features$id <- paste(data.features$id, splits, sep = "_")
    unique.splits <- unique(x = splits)
    id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)), 
                        "_", rep(x = unique(x = splits), times = length(x = id.levels)))
  }
  data.plot <- lapply(X = unique(x = data.features$id), FUN = function(ident) {
    data.use <- data.features[data.features$id == ident, 
                              1:(ncol(x = data.features) - 1), drop = FALSE]
    avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
      return(mean(x = expm1(x = x)))
    })
    pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove, 
                     threshold = 0)
    return(list(avg.exp = avg.exp, pct.exp = pct.exp))
  })
  names(x = data.plot) <- unique(x = data.features$id)
  data.plot <- lapply(X = names(x = data.plot), FUN = function(x) {
    data.use <- as.data.frame(x = data.plot[[x]])
    data.use$features.plot <- rownames(x = data.use)
    data.use$id <- x
    return(data.use)
  })
  data.plot <- do.call(what = "rbind", args = data.plot)
  if (!is.null(x = id.levels)) {
    data.plot$id <- factor(x = data.plot$id, levels = id.levels)
  }
  
  # could continue here with other kinds of data too
  
  avg.exp.scaled <- sapply(X = unique(x = data.plot$features.plot), 
                           FUN = function(x) {
                             data.use <- data.plot[data.plot$features.plot == 
                                                     x, "avg.exp"]
                             data.use <- scale(x = data.use)
                             data.use <- MinMax(data = data.use, min = col.min, 
                                                max = col.max)
                             return(data.use)
                           })
  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
  if (!is.null(x = split.by)) {
    avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled, 
                                         breaks = 20))
  }
  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$features.plot <- factor(x = data.plot$features.plot, 
                                    levels = rev(x = features))
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  if (!is.null(x = split.by)) {
    splits.use <- vapply(X = strsplit(x = as.character(x = data.plot$id), 
                                      split = "_"), FUN = "[[", FUN.VALUE = character(length = 1L), 
                         2)
    data.plot$colors <- mapply(FUN = function(color, value) {
      return(colorRampPalette(colors = c("grey", color))(20)[value])
    }, color = cols[splits.use], value = avg.exp.scaled)
  }
  color.by <- ifelse(test = is.null(x = split.by), yes = "avg.exp.scaled", 
                     no = "colors")
  if (!is.na(x = scale.min)) {
    data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
  }
  if (!is.na(x = scale.max)) {
    data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
  }
  plot <- ggplot(data = data.plot, mapping = aes_string(x = "id", y = "features.plot")) + 
    geom_point(mapping = aes_string(size = "pct.exp", color = color.by)) + 
    scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) + 
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + 
    guides(size = guide_legend(title = "Percent Expressed")) + 
    labs(x = "", y = ifelse(test = is.null(x = split.by), yes = "", no = "Split Identity")) + 
    theme_cowplot() +
    theme(axis.text.x=element_text(angle=45, hjust=1)) + 
    theme(axis.text.y = element_text(colour="grey20",size=10,face="plain", family="Helvetica"), axis.text.x = element_text(colour="grey20",size=10,face="plain", family="Helvetica")) +
    theme(legend.key=element_blank(), legend.box.background=element_blank(),
           legend.text=element_text(size = 10, family="Helvetica"))+
    theme(plot.margin=unit(c(1,1,1,1), "mm"))
          
  
  if (!is.null(x = split.by)) {
    plot <- plot + scale_color_identity()
  }
  else if (length(x = cols) == 1) {
    plot <- plot + scale_color_distiller(palette = cols)
  }
  else {
    plot <- plot + scale_color_gradient(low = cols[1], high = cols[2])
  }
  if (is.null(x = split.by)) {
    plot <- plot + guides(color = guide_colorbar(title = "Average Expression"))
  }
  return(plot)
}

PercentAbove <- function(x, threshold) {
  return(length(x = x[x > threshold]) / length(x = x))
}

RandomName <- function(length = 5L, ...) {
  CheckDots(..., fxns = 'sample')
  return(paste(sample(x = letters, size = length, ...), collapse = ''))
}

MinMax=function (data, min, max) {
  data2 <- data
  data2[data2 > max] <- max
  data2[data2 < min] <- min
  return(data2)
}

scale.func <- function (name = waiver(), breaks = waiver(), labels = waiver(), 
                        limits = NULL, range = c(1, 6), trans = "identity", guide = "legend") 
{
  continuous_scale("size", "radius", scales::rescale_pal(range), name = name, 
                   breaks = breaks, labels = labels, limits = limits, trans = trans, 
                   guide = guide)
}

#' Contains columns: variable.1, variable.2, features (y name), id (x name)
#' @examples
#' #'\dontrun{
#'df=data.frame("variable.1"=c(-2, 4), "variable.2"=c(10,40), "features"=c("feature1", "feature2"), "id"=c("group1", "group2"))
#'plot.DotPlot.df(df)
#'# example, modified:
#'plot.DotPlot.df(df, name.variable.1 = "Fold-Change", name.variable.2 = "FDR (-log10)", cols = c("blue","white", "red"), col.min = -2, col.max = 2, scale.min = 0, scale.max = 20, dot.scale = 4)
#'}
plot.DotPlot.df=function(data.plot, scale.min=NULL, scale.max=NULL, col.min=NULL,col.max=NULL, dot.scale = 6, fontsize=8, cols = c("lightgrey", "red"), name.variable.1="variable.1", name.variable.2="variable.2", number.legend.points=NULL){
  library(cowplot)

  # set scales for data:
  if(is.null(col.min))col.min=min(data.plot$variable.1)
  if(is.null(col.max))col.max=max(data.plot$variable.1)
  if(is.null(scale.min))scale.min=min(data.plot$variable.2)
  if(is.null(scale.max))scale.max=max(data.plot$variable.2)
  
  if(!is.null(number.legend.points)){
    breaks=seq(scale.min, scale.max, length.out=number.legend.points)
  }else{
    breaks=waiver()
  }
  
  
  # set scales:
  data.plot$variable.1 <- MinMax(data = data.plot$variable.1, min = col.min, max = col.max)
  data.plot$variable.2 <- MinMax(data = data.plot$variable.2, min = scale.min, max = scale.max)
  
  # rearrenge, feature1 on top:
  data.plot$features=factor(data.plot$features, levels=rev(unique(data.plot$features)))
  
  plot <- ggplot(data = data.plot, mapping = aes_string(x = "id", y = "features")) + 
    geom_point(mapping = aes_string(size = "variable.2", color = "variable.1")) + 
    scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max), breaks = breaks) + 
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + 
    guides(size = guide_legend(title = name.variable.2)) + 
    guides(colour=guide_legend(title = name.variable.1)) + 
    labs(x = "", y = "") + 
    theme_cowplot() +
    theme(axis.text.x=element_text(angle=45, hjust=1)) + 
    theme(axis.text.y = element_text(colour="grey20",size=fontsize,face="plain", family="Helvetica"), axis.text.x = element_text(colour="grey20",size=fontsize,face="plain", family="Helvetica")) +
    theme(legend.key=element_blank(), legend.box.background=element_blank(),
          legend.text=element_text(size=fontsize, family="Helvetica"))+
    theme(plot.margin=unit(c(1,1,1,1), "mm"))

  # 2 color gradient
  if(length(cols)==2)plot <- plot + scale_colour_gradient(low = cols[1], high = cols[2])
  
  # 3 color gradient
  # if(length(cols)==3)plot <- plot + scale_colour_gradient2(low = cols[1],mid = cols[2], high = cols[3], midpoint = 0)
  
  # forces to range (no midpoint, but middle color will be midpoint):
  if(length(cols)>2)plot <- plot + scale_colour_gradientn(colours = cols, limits = c(col.min, col.max))
  
  return(plot)
}