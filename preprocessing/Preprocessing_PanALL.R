# PECAN ALL data processing:
# https://pecan.stjude.cloud/proteinpaint/study/PanALL
files=list.files("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/Pecan_ALL", "HTSeq$", full.names = T)

m=do.call(cbind, parallel::mclapply(files, read.delim, header=F, row.names=1, mc.cores=8))

colnames(m)=gsub(".HTSeq|_.*.", "", list.files("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/Pecan_ALL", "HTSeq$", full.names = F))

save(m, file="PECAN_ALL_counts.Rdata")


# load data:
load("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/Pecan_ALL/PECAN_ALL_counts.Rdata")

# convert ENSAMBL to symbol:

# filter:
keep <- rowSums(edgeR::cpm(m) > 1) >= ceiling(dim(m)[2]*0.025)

m=m[keep,]

# map IDs:
library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl", host="useast.ensembl.org"))
genes <- rownames(m)
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)

G_list=G_list[!(G_list$hgnc_symbol%in%""|duplicated(G_list$hgnc_symbol)),]
m=m[match(G_list$ensembl_gene_id, rownames(m)),]

rownames(m)=G_list$hgnc_symbol

# load existing subtype coords:
annot=data.table::fread("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/Pecan_ALL/clinical_pecan.txt", data.table = F,dec = ",")
annot=annot[match(gsub("_.*.", "", colnames(gexp)), annot$patient),]

# normalize library size
DGE <- edgeR::DGEList(m)
DGE <- edgeR::calcNormFactors(DGE,method =c("TMM"))

# voom:
DGE.voom=limma::voom(DGE, plot=T)$E

# run combat:
# library(sva)
# annot$batch = annot$`RNA-seq library`
# modcombat = model.matrix(~1, data=annot)
# gexp = ComBat(dat=as.matrix(DGE.voom), batch=annot$batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

gexp=limma::removeBatchEffect(x = DGE.voom, batch = annot$`RNA-seq library`)

# write data out
coordinates.subtype=CancerMap(data = t(as.matrix(gexp)), name = "Pecan_pre-B-ALL", VAR = 10, BW = 1.75, perplexity = 30, PATH_OUTPUT = "/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/")
coordinates.subtype$subtype=annot$`primary subtype`

save(list = c("gexp", "coordinates.subtype", "annot"), file="/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/PecanALL_subtypes.Rdata")


coord_pecan=CancerMap(data = t(as.matrix(gexp)), name = "Pecan_pre-B-ALL", VAR = 10, BW = 1.75, perplexity = 30, PATH_OUTPUT = "/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/")
plot.scatter(x=coord_pecan$x, y = coord_pecan$y, group =  coordinates.subtype$`primary subtype`, namev = coordinates.subtype$`primary subtype`, main = "Pecan ALL", rasterize = F, width = 70*2, height = 74*2, SIZE = 0.5)

plot.scatter(x=coord_pecan$x, y = coord_pecan$y, group =  annot$`RNA-seq library`, namev = annot$`RNA-seq library`, main = "Pecan ALL", rasterize = F, width = 70*2, height = 74*2, SIZE = 0.5)
plot.scatter(x=coord_pecan$x, y = coord_pecan$y, group =  annot$institute, namev = annot$institute, main = "Pecan ALL", rasterize = F, width = 70*2, height = 74*2, SIZE = 0.5)

mod=names(table(annot$protocol))[table(annot$protocol)<3]
annot$protocol[annot$protocol%in%mod]="other"
plot.scatter(x=coord_pecan$x, y = coord_pecan$y, group =  annot$protocol, namev = annot$protocol, main = "Pecan ALL", rasterize = F, width = 70*2, height = 74*2, SIZE = 0.5)

mod=names(table(annot$protocol))[table(annot$protocol)<3]
annot$protocol[annot$protocol%in%mod]="other"

plot.scatter(x=coord_pecan$x, y = coord_pecan$y, group =  batch, namev =batch, main = "Pecan ALL", rasterize = F, width = 70*2, height = 74*2, SIZE = 0.5)