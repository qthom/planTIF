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
