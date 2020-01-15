

source("001-batchReadTrackData.R")
source("003-metageneMatrix.R")
source("randomPositions.R")
source("005-drawMetagenePlot.R")
source("extractSummits.R")
source("metagenePipeline.R")

pNET_dir <- "DATA from Zhu 2018 (PMID 30374093)"
pNET_filenames <- list.files(pNET_dir, pattern=".bedgraph.gz$")
pNET_data <- suppressWarnings(batchReadTrackData(pNET_filenames, pNET_dir, format="bedGraph", seqinfo=seqinfo(Athaliana)))

TSS <- resize(TSS.both2, width = 1000 , fix = "center")
TSS <- TSS[width(TSS) == 1000]

matlist <- lapply(pNET_data,metageneMatrix, interval = TSS)

pdf("metagene_FPvsMock_sppRNAgenes.pdf")
drawMetagenePlot(matlist,vline = 500, title = "Mock_vs_FP_treatment",xlabel = "Distance from the TSS (bp)",ylabel = "Average pNET-Seq signal")
dev.off()

TSSc <- resize(TSS.ctrlboth2, width = 1000 , fix = "center")
TSSc <- TSSc[width(TSSc) == 1000]

matlist <- lapply(pNET_data,metageneMatrix, interval = TSSc)

pdf("metagene_FPvsMock_ctrlgenes.pdf")
drawMetagenePlot(matlist,vline = 500, title = "Mock_vs_FP_treatment",xlabel = "Distance from the TSS (bp)",ylabel = "Average pNET-Seq signal")
dev.off()




library(ggsignif)

pNET_FP <- import("FP TREATED DATA FROM Zhu 2018 (PMID 30374093)")
pNET_Mock <- import("MOCK TREATED DATA FROM Zhu 2018 (PMID 30374093)")


genes <- genes(txdb)

genes$match <- match(genes$gene_id,TSS.both2$gene_id,nomatch=0)
genes_sppRNA <- genes[genes$match > 0,]
shrinked_genes_sppRNA <- extend(genes_sppRNA,upstream=-250,downstream=-250)

genes$match <- NULL
genes$match <- match(genes$gene_id,TSS.ctrlboth2$gene_id,nomatch=0)
genes_nosppRNA <- genes[genes$match > 0,]
shrinked_genes_nosppRNA <- extend(genes_nosppRNA,upstream=-250,downstream=-250)



data <- list("FP" = pNET_FP, "Mock" = pNET_Mock)
Norm_nosppRNA <- getOverlappingScores_v2(shrinked_genes_nosppRNA,data)
Norm_sppRNA <- getOverlappingScores_v2(shrinked_genes_sppRNA,data)
Norm_sppRNA$FP_norm <- Norm_sppRNA$FP/width(Norm_sppRNA)
Norm_sppRNA$Mock_norm <- Norm_sppRNA$Mock/width(Norm_sppRNA)
Norm_nosppRNA$FP_norm <- Norm_nosppRNA$FP/width(Norm_nosppRNA)
Norm_nosppRNA$Mock_norm <- Norm_nosppRNA$Mock/width(Norm_nosppRNA)


sppRNAFP <- data.frame("pNET" = Norm_sppRNA$FP_norm, "geno" = as.factor("sppRNA genes FP"))
sppRNAMock <- data.frame("pNET" = Norm_sppRNA$Mock_norm, "geno" = as.factor("sppRNA genes Mock"))
ctrlFP <- data.frame("pNET" = Norm_nosppRNA$FP_norm, "geno" = as.factor("no sppRNA genes FP"))
ctrlMock <- data.frame("pNET" = Norm_nosppRNA$Mock_norm, "geno" = as.factor("no sppRNA genes Mock"))

dat <- rbind(sppRNAMock,sppRNAFP,ctrlMock,ctrlFP)

anno <- t.test(sppRNAFP$pNET,sppRNAMock$pNET)$p.value

ggplot(dat, aes(y=pNET,x=geno)) + labs(x = "Genes Sets") +  
    geom_boxplot() + theme_classic() + geom_hline(yintercept = 0, size = 1) + 
    geom_signif(comparisons = list(c("no sppRNA genes Mock", "sppRNA genes Mock"),c("sppRNA genes Mock","sppRNA genes FP"),c("no sppRNA genes Mock","no sppRNA genes FP")), map_signif_level=TRUE) 

#...................................................................
