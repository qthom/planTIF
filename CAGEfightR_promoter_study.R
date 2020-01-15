# Load the required libraries

library(CAGEfightR)						# version 0.99.0
library(TxDb.Athaliana.BioMart.plantsmart28)
txdb <- TxDb.Athaliana.BioMart.plantsmart28
library(BiocParallel)
register(MulticoreParam(4), default=T)
library(BSgenome.Athaliana.TAIR.TAIR9)
seqlevels(Athaliana) <- c("1", "2", "3", "4", "5", "Mt", "Pt")
library(tibble)
library(dplyr)
library(edgeR)
library(DESeq2)

source("/home/bcm215/Desktop/Quentin/scripts/Nielsen_et_al_2018-master/assignTxType_custom.R")

# Load the TSS-Seq BigWig files (see the 01-Alignment_of_5Cap-Seq_data.sh pipeline):

setwd("/home/bcm215/Desktop/sequencing_files/CAP-Seq/bedgraph_capseq/expanded_cov3_bw/")
bw_plus_filenames <- list.files(".", pattern="fw_cov3_expanded.bw$")
bw_minus_filenames <- list.files(".", pattern="rev_cov3_expanded.bw$")
bw_plus <- BigWigFileList(bw_plus_filenames)
bw_minus <- BigWigFileList(bw_minus_filenames)
sample_names <- sub('fw_cov3_expanded.bw', '', bw_plus_filenames)
names(bw_plus) <- sample_names
names(bw_minus) <- sample_names

# Make the design matrix:
design_matrix <- data.frame("Name"=sample_names, "BigWigPlus"=bw_plus_filenames, "BigWigMinus"=bw_minus_filenames,
row.names=sample_names)


# Quantify all tag clusters (TCs):
ctss <- quantifyCTSSs(plusStrand=bw_plus, minusStrand=bw_minus, design=design_matrix, genome=seqinfo(Athaliana))

# Call candidate TSS:
tss <- quickTSSs(ctss)


# Annotate TSS and enhancers by genomic features (observe that a custom assignTxType() function is used):
rowRanges(tss)$txType <- suppressWarnings(assignTxType_custom(rowRanges(tss), txdb=txdb, asFactor=TRUE))

# Combine candidate TSS and enhancers into a single RangedSummarizedExperiment object:
rowRanges(tss)$clusterType <- "TSS"

rse <- combineClusters(tss, tss, removeIfOverlapping="object1")

# Remove low expressed TCs:
rse <- subsetBySupport(rse, inputAssay = "counts", outputColumn = "support", unexpressed = 0, minSamples = 1) # n = 96232

# Annotate TCs by gene IDs:
rse <- suppressWarnings(assignGeneID(rse, geneModels=txdb))

# Annotate TCs by gene names:
tair_ann <- import.gff3("/home/bcm215/groupdirs/SCIENCE-PLEN-Marquardt_lab/Quentin/sequencing_files/TAIR10/Arabidopsis_thaliana.TAIR10.26.gff3")
mapping <- data.frame("geneID" = sub("gene:", "", tair_ann$ID), "name"=tair_ann$external_name)
tmp <- left_join(x=tibble(geneID=rowRanges(rse)$geneID), y=mapping, by="geneID", na_matches = "never")
mcols(rse) <- DataFrame(mcols(rse), tmp[,-1])
rm(tmp)