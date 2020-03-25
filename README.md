# ImmunogenomicLandscape-BloodCancers

## Description:
Scripts related to results presented in Dufva and Pölönen et al.

- Modify the **working directory** to the scripts
- setwd("path/data"), download data from 

- Modify the **GIT_HOME** variable in R script (folder where the git folder is cloned). GIT_HOME points to **common_scripts** that contains various tools that were used in the analysis

## Scripts to reproduce the analysis, plots and tables

#### Figure 1, Figure S1:
```
Fig1_plots.R (Figure 1F-H)
Statistical_analysis_Cytolytic_Score_development.R (FigureS1B, G-I)
```
#### Figure 2, Figure S2 and Table S2:
```
Fig2_microenvironment_analysis.R (TableS2, Fig2A-B, FigS2C)
FigS2AB_microenvironment_analysis_GSEA_GSVA.R (FigS2A-B)
FigS2DE_microenvironment_validation_CLL_AML.R (FigS2D-E)
```
#### Figure 3, Figure S3 and Table S3:
```
Fig3_plots_scRNA.R  (Fig3F-H and K, FigS3K and O)
Fig3_FigS3_MDSsignature.R (Fig3E, FigS3D)
FigS3L.R (TableS3 tab, FigS3L)
Statistical_analysis_DE_analysis_MDSlike.R (TableS3 tab)
Statistical_analysis_scRNA_MDSlike_analysis.R  (TableS3 tabs)
Statistical_analysis_Szabo_TCell_analysis.R (Fig3I, FigS3M)
Statistical_analysis_Yang_NKCell_analysis.R (Fig3J, FigS3N)
```
#### Figure 4, Figure S4 and Table S4:
```
Statistical_analysis_HLAII_Score_development.R (FigS4A-B)
Fig4D_TCGA_AML_complexheatmap_CIITA.R (Fig4D)
FigS4EF_CIITA_methylation_validation_ERRBS.R (FigS4E-F)
FigS4K_CCLE_CIITA_methylation.R (FigS4K)
FigS4J_CIITA_methylation_validation_GSE49031.R (FigS4J)
FigS4D_TCGA_AML_global_hypermethylation.R (FigS4D)
Statistical_analysis_FIMM_AML_RRBS.R (Fig4G and H, FigS4L)
```
#### Figure 5, Figure S5 and Table S5:
```
Fig5_DE_analysis_costim.R  (Fig5B, FigS5B, Table S5 tabs)
FigS5C.R (FigS5C)
```
#### Figure 6, Figure S6 and Table S6:
```
Fig6B_S6F_CGA_tSNEplot_hemap.R (Fig6B, FigS6F)
Fig6C_S6B_CGA_Hemap.R (Fig6C, FigS6B)
Fig6H_FigS6H_CGA_heatmap_GSE98588.R (Fig6H, FigS6H, FigS6I)
FigS6G_CGA_GSEA_hemapMM.R (FigS6G)
Fig6D_CCLE_CGA_methylation.R (Fig6D)
Statistical_analysis_CGA_discovery_Hemap.R
```
#### Figure 7, Figure S7 and Table S7:
```
Fig7_Univariate_Coxph_survival.R (TableS7)
Fig7_univariate_survival_forestplot.R (Fig7A, Fig7B-D, FigS7A-B, )
Fig7_multivariable_regression_eNet_survival.R (Fig7E-F, FigS7C-G)
```
#### Statistical association analysis, Supplementary tables 1-6
```
Statistical_analysis_correlations_BeatAML.R (TableS3-5)
Statistical_analysis_correlations_CoMMpass.R (TableS4-6)
Statistical_analysis_correlations_GSE98588_DLBCL.R (TableS3-6)
Statistical_analysis_correlations_Reddy.R (TableS3-5)
Statistical_analysis_correlations_TCGA_AML.R (TableS3-5)
Statistical_analysis_correlations_TCGA_DLBCL.R
```
## Scripts related to data preprocessing (under folder preprocessing)
These scripts are for reference only. Raw data and input data would have to be downloaded and processed for these scripts. Check the publication for raw data accession codes.

#### Data Preprocessing
```
Preprocessing_CIBERSORT_MCPcounter.R
Preprocessing_normalize_hguarray_GSE98588.R
```
#### Generating featurematrix for each cohort
```
Preprocessing_Hemap_featurematrix_generation.R (TableS1)
Preprocessing_REDDY_DLBCL_featurematrix_generation.R
Preprocessing_TCGA_AML_featurematrix_generation.R
Preprocessing_TCGA_DLBCL_featurematrix_generation.R
Preprocessing_coMMpass_featurematrix_generation.R
Preprocessing_GSE98588_DLBCL_featurematrix_generation.R
Preprocessing_PanALL.R
```
#### Subtype analysis for each cohort
```
Preprocessing_MM_subtyping.R
Preprocessing_DLBCL_subtyping.R
Preprocessing_ALL_subtyping.R
Preprocessing_AML_subtyping.R
```
#### Methylation data processing
```
Preprocessing_TCGA_AML_add_meth_probes.R
Preprocessing_TCGA_meth_data_genelist.R
Preprocessing_FIMM_AML_RRBS.R
Preprocessing_AML_RRBS_meth_de_analysis.R
```
#### scRNA data preprocessing for statistical analysis
```
Preprocessing_scRNA_CLL_GSE111014.R
Preprocessing_scRNA_FIMM_AML.R
Preprocessing_scRNA_Galen_AML.R
Preprocessing_scRNA_HCA.R
Preprocessing_scRNA_PB_Citeseq.R
Preprocessing_scRNA_Szabo_Tcells_dataprocessing.R
Preprocessing_scRNA_Yang_NK_dataprocessing.R
Preprocessing_scRNA_integrate_Tcells_NKcells.R
```


## SessionInfo (R)
```
R version 3.6.0 (2019-04-26)
Platform: x86_64-redhat-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

Matrix products: default
BLAS/LAPACK: /usr/lib64/R/lib/libRblas.so

locale:
 [1] LC_CTYPE=C                 LC_NUMERIC=C
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
 [9] LC_ADDRESS=C               LC_TELEPHONE=C
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

attached base packages:
 [1] grid      parallel  stats4    stats     graphics  grDevices utils
 [8] datasets  methods   base

other attached packages:
 [1] tidyr_1.0.0
 [2] dplyr_0.8.3
 [3] cowplot_1.0.0
 [4] tibble_2.1.3
 [5] stringr_1.4.0
 [6] survminer_0.4.6
 [7] survcomp_1.34.0
 [8] prodlim_2018.04.18
 [9] Seurat_3.1.1.9021
[10] RnBeads_2.2.0
[11] plyr_1.8.4
[12] methylumi_2.30.0
[13] minfi_1.30.0
[14] bumphunter_1.26.0
[15] locfit_1.5-9.1
[16] iterators_1.0.12
[17] foreach_1.4.7
[18] Biostrings_2.52.0
[19] XVector_0.24.0
[20] SummarizedExperiment_1.14.1
[21] DelayedArray_0.10.0
[22] BiocParallel_1.18.1
[23] FDb.InfiniumMethylation.hg19_2.2.0
[24] org.Hs.eg.db_3.8.2
[25] TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2
[26] GenomicFeatures_1.36.4
[27] reshape2_1.4.3
[28] scales_1.0.0
[29] illuminaio_0.26.0
[30] matrixStats_0.54.0
[31] gridExtra_2.3
[32] gplots_3.0.1.1
[33] fields_10.0
[34] maps_3.3.0
[35] spam_2.4-0
[36] dotCall64_1.0-0
[37] ff_2.2-14
[38] bit_1.1-14
[39] cluster_2.1.0
[40] MASS_7.3-51.4
[41] GenomicRanges_1.36.0
[42] GenomeInfoDb_1.20.0
[43] RColorBrewer_1.1-2
[44] openxlsx_4.1.4
[45] multipanelfigure_2.0.2
[46] mclust_5.4.5
[47] Matrix_1.2-17
[48] Hmisc_4.2-0
[49] Formula_1.2-3
[50] survival_2.44-1.1
[51] lattice_0.20-38
[52] GSVA_1.32.0
[53] ggridges_0.5.1
[54] ggrastr_0.1.7
[55] ggpubr_0.2.3
[56] future_1.14.0
[57] forestplot_1.9
[58] checkmate_1.9.4
[59] magrittr_1.5
[60] EnhancedVolcano_1.3.5
[61] ggrepel_0.8.1
[62] ggplot2_3.2.1
[63] edgeR_3.26.8
[64] limma_3.40.6
[65] data.table_1.12.4
[66] ComplexHeatmap_2.0.0
[67] circlize_0.4.7
[68] caTools_1.17.1.2
[69] AnnotationDbi_1.46.1
[70] IRanges_2.18.2
[71] S4Vectors_0.22.0
[72] Biobase_2.44.0
[73] BiocGenerics_0.30.0
[74] methylSig_0.1
```
