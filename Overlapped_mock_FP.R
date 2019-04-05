

source("001-batchReadTrackData.R")
source("003-metageneMatrix.R")
source("randomPositions.R")
source("005-drawMetagenePlot.R")
source("extractSummits.R")
source("metagenePipeline.R")

pNET_dir <- "directory for the pNET-Seq FP treated and Mock"
pNET_filenames <- list.files(pNET_dir, pattern=".bedgraph.gz$")
pNET_data <- suppressWarnings(batchReadTrackData(pNET_filenames, pNET_dir, format="bedGraph", seqinfo=seqinfo(Athaliana)))

TSS <- resize(TSS.both2, width = 1000 , fix = "center")
TSS <- TSS[width(TSS) == 1000]

TSSc <- resize(TSS.ctrlboth2, width = 1000 , fix = "center")
TSSc <- TSS[width(TSS) == 1000]

matlist_sppRNA <- lapply(pNET_data,metageneMatrix, interval = TSS)
matlist_ctrl <- lapply(pNET_data,metageneMatrix, interval = list)


drawMetagenePlot(matlist_sppRNA,vline = 500, title = "Mock_vs_FP_treatment",xlabel = "Distance from the TSS (bp)",ylabel = "Average pNET-Seq signal")

drawMetagenePlot(matlist_ctrl,vline = 500, title = "Mock_vs_FP_treatment",xlabel = "Distance from the TSS (bp)",ylabel = "Average pNET-Seq signal")



