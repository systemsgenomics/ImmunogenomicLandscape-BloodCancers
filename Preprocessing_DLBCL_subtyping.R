GIT_HOME="/research/users/ppolonen/git_home/common_scripts"
source(file.path(GIT_HOME, "visualisation/plotting_functions.R"))
source(file.path(GIT_HOME, "statistics/functions_statistics.R"))
source(file.path(GIT_HOME, "pathway_analysis/functions.GSEA.R"))

setwd("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/")

annot = get(load("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/Hemap_immunology_Annotations.Rdata"))

# gexp data
data=t(get(load("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/data9544_with_gene_symbols.RData")))
data=data[,colnames(data)%in%annot$GSM.identifier..sample.]


gexp=data[,annot$tbLY%in%c("Lymphoma_BCL_DLBCL_ABC", "Lymphoma_BCL_DLBCL_GCB")]

coordinates.subtype=gsub("Lymphoma_BCL_DLBCL_", "", annot$tbLY[annot$tbLY%in%c("Lymphoma_BCL_DLBCL_ABC", "Lymphoma_BCL_DLBCL_GCB")])

save(list = c("gexp", "coordinates.subtype"), file = "Hemap_DLBCL_subtypes.Rdata")
