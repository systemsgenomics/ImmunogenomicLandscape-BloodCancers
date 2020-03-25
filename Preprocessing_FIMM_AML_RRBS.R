# sbatch
# https://www.sciencedirect.com/science/article/pii/S0168165617315936
# RnBeads chosen as it uses limma style statistics suitable if you have replicates.

# BiocManager::install("RnBeads")
# BiocManager::install("RnBeads.hg38")

library(RnBeads)
logger.start(fname=NA)
parallel.setup(8)
parallel.isEnabled()

# set dirs:
data.dir <- file.path(getwd(), "input_data")
analysis.dir=file.path(getwd(), "analysis")
report.dir <-file.path(getwd(), "report_ttest")

# Set some analysis options
rnb.options(
  filtering.sex.chromosomes.removal=T,
  identifiers.column="SampleID",
  replicate.id.column="treatment",
  import.bed.style="bismarkCov",
  differential.site.test.method = "ttest",
  differential.enrichment.go=TRUE,
  differential.enrichment.lola=TRUE,
  filtering.coverage.threshold=10,
  assembly="hg38",
  differential.report.sites=FALSE,
  filtering.sex.chromosomes.removal = TRUE,
  import.table.separator="\t"
)

# add CpGisland shore and shelves to analysis:
d=data.table::fread("cpgIslandExt.hg38.bed", data.table = F, fill = T)

CpG.expanded <- data.frame(
  Chromosome = d$V1,
  Start = as.integer(d$V2-4000),
  End = as.integer(d$V3+4000),
  row.names = paste(make.unique(d$V4), "island.shore.shelf", sep="@"), stringsAsFactors = F)

shelf.left <- data.frame(
  Chromosome = d$V1,
  Start = as.integer(d$V2-2000),
  End = as.integer(d$V2),
  row.names = paste(make.unique(d$V4), "shelf.l", sep="@"), stringsAsFactors = F)

shelf.right <- data.frame(
  Chromosome = d$V1,
  Start = as.integer(d$V3),
  End = as.integer(d$V3+2000),
  row.names = paste(make.unique(d$V4), "shelf.r", sep="@"), stringsAsFactors = F)

shore.left <- data.frame(
  Chromosome = d$V1,
  Start = as.integer(d$V2-4000),
  End = as.integer(d$V2-2000),
  row.names = paste(make.unique(d$V4), "shore.l", sep="@"), stringsAsFactors = F)

shore.right <- data.frame(
  Chromosome = d$V1,
  Start = as.integer(d$V3+2000),
  End = as.integer(d$V3+4000),
  row.names = paste(make.unique(d$V4), "shore.r", sep="@"), stringsAsFactors = F)

cpg <- data.frame(
  Chromosome = d$V1,
  Start = as.integer(d$V2),
  End = as.integer(d$V3),
  row.names = paste(make.unique(d$V4), "island", sep="@"), stringsAsFactors = F)


CpG_shelf_shore=rbind(cpg, shelf.left, shelf.right, shore.left, shore.right)
CpG_shelf_shore=CpG_shelf_shore[!(CpG_shelf_shore[,2]<0|CpG_shelf_shore[,3]<0),]

txt <- "CpG islands + CpG shelves + CpG shores"
rnb.set.annotation(type = "CpGregion", regions = CpG_shelf_shore, description = txt, assembly = "hg38")

rnb.save.annotation("CpGregion_annotations.Rdata", "CpGregion", assembly = "hg38")

write.table(cbind(CpG_shelf_shore[,c(1:3)], rownames(CpG_shelf_shore)), "CpG_shelf_shore.bed", quote = F, row.names = F, col.names = F, sep="\t")


rnb.set.annotation(type = "CpG.expanded", regions = CpG.expanded, description = txt, assembly = "hg38")

rnb.save.annotation("CpGexpanded_annotations.Rdata", "CpG.expanded", assembly = "hg38")

# CIITA CpG
CIITA_CpG=CpG_shelf_shore[grep("CpG:41.91|CpG:21.172|CpG:20.179",rownames(CpG_shelf_shore)),]
rnb.set.annotation(type = "CpGregion_CIITA", regions = CIITA_CpG, description = txt, assembly = "hg38")
rnb.save.annotation("CpGregion_CIITA_annotations.Rdata", "CpGregion_CIITA", assembly = "hg38")

# options(fftempdir=file.path(getwd(), "temp"))
# getOption("fftempdir")

fileList=list.files(data.dir, pattern = ".cov", full.names = T)
sample.id=gsub("_FRB.*.-1a.bismark.cov|M13_","", basename(fileList))
treatment=gsub("M13_|_1$|_2$|_3$", "", sample.id)
fileList=list.files(data.dir, pattern = ".cov", full.names = F)

# compare DC treated to DMSO/IFN treated
DAC_REST=ifelse(grepl("DC", treatment), "DC/DCIFN", "DMSO/IFN")

# CIITA 
pData = data.frame(
  files=fileList,
  SampleID = sample.id,
  treatment = treatment,
  row.names = sample.id,
  DAC_REST,
  stringsAsFactors = FALSE)

sample.annotation <- file.path(getwd(), "pData.txt")

# write pData
write.table(pData, sample.annotation, quote = F, sep = "\t", row.names = F, col.names = T)

data.source <- c(data.dir, sample.annotation)

# initialize
rnb.initialize.reports(report.dir)

result=rnb.run.import(data.source,data.type = "bs.bed.dir", dir.reports=report.dir, init.configuration = !file.exists(file.path(report.dir, "configuration")), close.report = TRUE, show.report = FALSE)
rnb.set <- result$rnb.set

## Quality Control
rnb.run.qc(rnb.set, report.dir)

rnb.set <- rnb.run.preprocessing(rnb.set, dir.reports=report.dir)$rnb.set

save.dir <- file.path(report.dir, "analysis")
save.rnb.set(rnb.set, path=save.dir, archive=TRUE)

rnb.options("differential.variability"=TRUE)
dm <- rnb.execute.computeDiffMeth(rnb.set,pheno.cols=c("DAC_REST"))

save.rnb.diffmeth(dm, save.dir)

## Data export
rnb.run.tnt(rnb.set, report.dir)



# OLD 
# sample.id=gsub("_FRB.*.-1a.bismark.cov|M13_","", basename(fileList))
# treatment=gsub("M13_|_1$|_2$|_3$", "", sample.id)
# fileList=list.files(data.dir, pattern = ".cov", full.names = F)
# 
# DAC_DMSO=treatment
# DAC_DMSO[!DAC_DMSO%in%c("Ctrl_DC", "Ctrl_DMSO")]=NA
# 
# IFNG_DMSO=treatment
# IFNG_DMSO[!IFNG_DMSO%in%c("Ctrl_IFN", "Ctrl_DMSO")]=NA
# 
# DAC_IFNG_DMSO=treatment
# DAC_IFNG_DMSO[!DAC_IFNG_DMSO%in%c("Ctrl_DCIFN", "Ctrl_DMSO")]=NA
# 
# # compare DAC + IFNG stimulus when sgCIITA and wtCIITA
# sgCIITA_DMSO=treatment
# sgCIITA_DMSO[!sgCIITA_DMSO%in%c("CIITA_DMSO", "Ctrl_DMSO")]=NA
# 
# sgCIITA_DAC_DMSO=treatment
# sgCIITA_DAC_DMSO[!sgCIITA_DAC_DMSO%in%c("CIITA_DC", "Ctrl_DC")]=NA
# 
# sgCIITA_IFNG_DMSO=treatment
# sgCIITA_IFNG_DMSO[!sgCIITA_IFNG_DMSO%in%c("CIITA_IFN", "Ctrl_IFN")]=NA
# 
# sgCIITA_DAC_IFNG_DMSO=treatment
# sgCIITA_DAC_IFNG_DMSO[!sgCIITA_DAC_IFNG_DMSO%in%c("CIITA_DCIFN", "Ctrl_DCIFN")]=NA
# 
# # compare DAC + IFNG stimulus when sgTET2 and wtTET2
# gTET2_DMSO=treatment
# gTET2_DMSO[!gTET2_DMSO%in%c("TET2_DMSO", "Ctrl_DMSO")]=NA
# 
# sgTET2_DAC_DMSO=treatment
# sgTET2_DAC_DMSO[!sgTET2_DAC_DMSO%in%c("TET2_DC", "Ctrl_DC")]=NA
# 
# sgTET2_IFNG_DMSO=treatment
# sgTET2_IFNG_DMSO[!sgTET2_IFNG_DMSO%in%c("TET2_IFN", "Ctrl_IFN")]=NA
# 
# sgTET2_DAC_IFNG_DMSO=treatment
# sgTET2_DAC_IFNG_DMSO[!sgTET2_DAC_IFNG_DMSO%in%c("TET2_DCIFN", "Ctrl_DCIFN")]=NA
# 
# DAC_REST=ifelse(grepl("DC", treatment), "DC", "REST")
# 
# DACIFN_DMSOIFN=ifelse(grepl("DCIFN", treatment), "DCIFN", "DMSO_IFN")
# DACIFN_DMSOIFN[grepl("_DC$", treatment)]=NA
