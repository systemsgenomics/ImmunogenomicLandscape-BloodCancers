GSE="GSE116256"

library(GEOquery)
gse <- getGEO(GSE) # retrieves a GEO list set for your SRA id.

## see what is in there:
show(gse)

##  what you want is table a with SRR to download and some sample information:
## lets see what the first set contains:
df <- as.data.frame(gse[[1]])
head(df)

df.filt=df[grepl("D0$", df$title),]

file.list=gsub(".*.=", "", as.character(df.filt$relation.1))

data.table::fwrite(list(file.list), file="accList.txt")

#**********************************
# run in shell, use screen or nohup, version must be prefetch.2.9.0
# This command is important to change the default download folder... change the path accordingly
# echo "/repository/user/main/public/root = \"/research/groups/sysgen/raw_data/petri/scRNA/Galen_AML_GSE116256\"" > $HOME/.ncbi/user-settings.mkfg
# prefetch $(<accList.txt) --max-size 25000000
# fastq-dump sra/*.sra
#**********************************

# make input pdata.txt file

old_name=df.filt$title
new_name=df.filt$title
project="Hemaimmunology"
assay="scRNAseq"
gseid=GSE
gsmid=df.filt$geo_accession
celltype=df.filt$source_name_ch1
treatment="na"
time="na"
replicate=1
genome_fromversion="hg19"
folder="align_scRNA"


data.table::fwrite(data.frame(old_name, new_name, project, assay, gseid, celltype, treatment, time,replicate,genome_fromversion,folder), file="pdata.txt", sep="\t")



file.list[!file.list%in%gsub(".sra", "", list.files("/research/groups/sysgen/raw_data/petri/scRNA/Galen_AML_GSE116256/sra", ".sra"))]




# library(SRAdb)
# sqlfile = "/research/work/ppolonen/SRAmetadb.sqlite" # obtained by getSRAdbFile()
# sra_con <- dbConnect(SQLite(),sqlfile)
