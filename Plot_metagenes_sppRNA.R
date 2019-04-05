

#pNET-Seq metagene plot

source("metagenePipeline.R")
seqlevels(Athaliana) <- c("1", "2", "3", "4", "5", "Mt", "Pt")
seqlevels(TSS.both2) <- c("1", "2", "3", "4", "5")
seqlevels(TSS.ctrlboth2) <- c("1", "2", "3", "4", "5")
seqlevels(TTS_hen2) <- c("1", "2", "3", "4", "5")

nub <- length(TSS.both2)

na1 <- paste("Control genes without detectable sppRNAs,n=",nub)
grl_pnet <- list("genes without sppRNAs (n=602)" = TSS.ctrlboth2, "genes with sppRNAs,n=602" = TSS.both2)
grl_PAS <- list("PAS of genes with sppRNAs (n=602)" = TTS_hen2)


yaxlab <- "Average COL.0 pNET-Seq signal (with 95% CI)"
pNET_dir <- "bedgraph for unphosphorilated pNET-Seq data from zhu et al"
pNET_filenames <- list.files(pNET_dir, pattern=".bedgraph.gz$")
pNET_data <- suppressWarnings(batchReadTrackData(pNET_filenames, pNET_dir, format="bedGraph", seqinfo=seqinfo(Athaliana)))
r2 <- suppressWarnings(metagenePipeline(signal=pNET_data, intervals=grl_pnet, out_dir="./metagene", expand=500))
source("metagenePipeline_v3.R")
r2 <- suppressWarnings(metagenePipeline(signal=pNET_data, intervals=grl_PAS, out_dir="./metagene/PAS", expand=500))

rm(pNET_data)


