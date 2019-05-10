# Load the required custom code:
setwd("/home/bcm215/Desktop/Quentin/scripts")

#pNET-Seq metagene plot

source("/home/bcm215/Desktop/Quentin/scripts/TIF-Seq_2019/metagenePipeline.R")
seqlevels(Athaliana) <- c("1", "2", "3", "4", "5", "Mt", "Pt")
seqlevels(TSS.both2) <- c("1", "2", "3", "4", "5")
seqlevels(TSS.ctrlboth2) <- c("1", "2", "3", "4", "5")
seqlevels(TTS_hen2) <- c("1", "2", "3", "4", "5")

nub <- length(TSS.both2)

na1 <- paste("Control genes without detectable sppRNAs,n=",nub)
grl_pnet <- list("genes without sppRNAs (n=1153)" = TSS.ctrlboth2, "genes with sppRNAs,n=1155" = TSS.both2)
grl_PAS <- list("PAS of genes with sppRNAs (n=1155)" = TTS_hen2)


yaxlab <- "Average COL.0 pNET-Seq signal (with 95% CI)"
pNET_dir <- "/home/bcm215/groupdirs/SCIENCE-PLEN-Marquardt_lab/The Holy Grail/Arabidopsis/Zhu 2018 (PMID 30374093) - pNET-Seq/pNET-Seq/Norm1M/f/"
pNET_filenames <- list.files(pNET_dir, pattern=".bedgraph.gz$")
pNET_data <- suppressWarnings(batchReadTrackData(pNET_filenames, pNET_dir, format="bedGraph", seqinfo=seqinfo(Athaliana)))
r2 <- suppressWarnings(metagenePipeline(signal=pNET_data, intervals=grl_pnet, out_dir="./plots", expand=500))
source("/home/bcm215/Desktop/Quentin/scripts/Nielsen_et_al_2018-master/metagenePipeline_v3.R")
r2 <- suppressWarnings(metagenePipeline(signal=pNET_data, intervals=grl_PAS, out_dir="./plots/PAS", expand=500))

rm(pNET_data)


