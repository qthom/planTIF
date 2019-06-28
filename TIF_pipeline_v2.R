
#In this script all functions are called for producing the plots and analysis of the paper

#the version used in the analysis are as follow:  [1] CAGEfightR_1.2.0                         
 [2] itsadug_2.3                              
 [3] plotfunctions_1.3                        
 [4] mgcv_1.8-28                              
 [5] nlme_3.1-140                             
 [6] LSD_4.0-0                                
 [7] VennDiagram_1.6.20                       
 [8] futile.logger_1.4.3                      
 [9] tibble_2.1.3                             
[10] SummarizedExperiment_1.12.0              
[11] DelayedArray_0.8.0                       
[12] BiocParallel_1.16.6                      
[13] matrixStats_0.54.0                       
[14] scales_1.0.0                             
[15] reshape_0.8.8                            
[16] BSgenome.Athaliana.TAIR.TAIR9_1.3.1000   
[17] BSgenome_1.50.0                          
[18] Biostrings_2.50.2                        
[19] XVector_0.22.0                           
[20] biomaRt_2.38.0                           
[21] ggplot2_3.1.1                            
[22] rtracklayer_1.42.2                       
[23] TxDb.Athaliana.BioMart.plantsmart28_3.2.2
[24] GenomicFeatures_1.34.8                   
[25] AnnotationDbi_1.44.0                     
[26] Biobase_2.42.0                           
[27] GenomicRanges_1.34.0                     
[28] GenomeInfoDb_1.18.2                      
[29] IRanges_2.16.0                           
[30] S4Vectors_0.20.1                         
[31] BiocGenerics_0.28.0                      
[32] RevoUtils_11.0.1                         
[33] RevoUtilsMath_11.0.0                     

loaded via a namespace (and not attached):
 [1] ProtGenerics_1.14.0       bitops_1.0-6             
 [3] bit64_0.9-7               RColorBrewer_1.1-2       
 [5] progress_1.2.2            httr_1.4.0               
 [7] backports_1.1.4           tools_3.5.1              
 [9] R6_2.4.0                  rpart_4.1-15             
[11] Hmisc_4.2-0               DBI_1.0.0                
[13] lazyeval_0.2.2            Gviz_1.26.5              
[15] colorspace_1.4-1          nnet_7.3-12              
[17] withr_2.1.2               gridExtra_2.3            
[19] tidyselect_0.2.5          prettyunits_1.0.2        
[21] curl_3.3                  bit_1.1-14               
[23] compiler_3.5.1            htmlTable_1.13.1         
[25] formatR_1.6               checkmate_1.9.3          
[27] stringr_1.4.0             digest_0.6.19            
[29] Rsamtools_1.34.1          foreign_0.8-71           
[31] dichromat_2.0-0           htmltools_0.3.6          
[33] base64enc_0.1-3           pkgconfig_2.0.2          
[35] ensembldb_2.6.8           htmlwidgets_1.3          
[37] rlang_0.3.4               rstudioapi_0.10          
[39] RSQLite_2.1.1             acepack_1.4.1            
[41] dplyr_0.8.1               VariantAnnotation_1.28.13
[43] RCurl_1.95-4.12           magrittr_1.5             
[45] GenomeInfoDbData_1.2.0    Formula_1.2-3            
[47] Matrix_1.2-17             Rcpp_1.0.1               
[49] munsell_0.5.0             stringi_1.4.3            
[51] zlibbioc_1.28.0           plyr_1.8.4               
[53] blob_1.1.1                crayon_1.3.4             
[55] lattice_0.20-38           splines_3.5.1            
[57] hms_0.4.2                 knitr_1.23               
[59] pillar_1.4.1              futile.options_1.0.1     
[61] XML_3.98-1.20             glue_1.3.1               
[63] biovizBase_1.30.1         latticeExtra_0.6-28      
[65] data.table_1.12.2         lambda.r_1.2.3           
[67] gtable_0.3.0              purrr_0.3.2              
[69] assertthat_0.2.1          xfun_0.7                 
[71] AnnotationFilter_1.6.0    survival_2.44-1.1        
[73] GenomicAlignments_1.18.1  memoise_1.1.0            
[75] cluster_2.0.9  

#the script take approximatly 20 minutes to run and require that you downloaded the processed bed files from the GEO submission in the right directory 

list.of.packages <- c("BiocParallel","CAGEfightR","itsadug","LSD","VennDiagram","tibble","SummarizedExperiment","scales","GenomicRanges","ggplot2","TxDb.Athaliana.BioMart.plantsmart28","rtracklayer","biomaRt","BSgenome.Athaliana.TAIR.TAIR9","reshape")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(GenomicRanges)
library(TxDb.Athaliana.BioMart.plantsmart28)
library(rtracklayer)
library(ggplot2)
library(biomaRt)
library(BSgenome.Athaliana.TAIR.TAIR9)
library(reshape)
library(scales)
library(SummarizedExperiment)
library(tibble)
library(VennDiagram)
library(LSD)
library(itsadug)
#library(CAGEfightR)						
#library(BiocParallel)


txdb <- TxDb.Athaliana.BioMart.plantsmart28

seqlevels(Athaliana) <- c("1", "2", "3", "4", "5", "Mt", "Pt")

#load functions called in the scripts


source("assignTxType_custom.R")
source("001-batchReadTrackData.R")
source("metageneMatrix.R")
source("drawMetagenePlot.R")
source("extractSummits.R")
source("metagenePipeline.R")
source("002-getOverlappingScores.R")
source("findMatchedControl_v2.R")

#Generate the objects needed for the pipeline (i.e genomic coordinates)

source("Process_TIFs.R")

#The following scripts plot the main figure

#loading dplyr now because of a function conflict with biomart (select())
library(dplyr)

#Calculates controls sets of genes

source("subset_equaldistribution_v4.R")


#histscatterplots for TUs - the plotting takes a while for that script

source("histscatter_all.R")
source("histscatter_FACT.R")
source("histscatter_small.R")


#Metagenes plots

source("Plot_metagenes_sppRNA.R")

#histograms

source("histogram.R")

#metagenePlot FP vs Mock

source("Overlapped_mock_FP.R")

#VennDiagram

source("VennDiagram.R")
