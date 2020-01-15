library(groupdata2)

coding_genes_n0big <- coding_genes_nO[width(coding_genes_nO) > 1000, ]
coding_genes_n0big$widthg <- width(coding_genes_n0big)

TSS <- resize(coding_genes_n0big, width = 1 , fix = "start")
extended_TSS <- extend(TSS,upstream=149,downstream=150)
shrunk_genes <- extend(coding_genes_n0big,upstream=-300,downstream=-300)

pNET <- import('/home/bcm215/groupdirs/SCIENCE-PLEN-Marquardt_lab/The Holy Grail/Arabidopsis/Zhu 2018 (PMID 30374093) - pNET-Seq/pNET-Seq/Norm1M/f/s03_s04_unph_fw_rev_norm1M.bedgraph.gz', format = 'bedgraph')

seqlevels(pnet) <- c("1", "2", "3", "4", "5")

data <- list("pNET" = pNET)
Norm <- getOverlappingScores_v2(shrunk_genes,data)
Norm$normpNET <- Norm$pNET/Norm$widthg
NormTSS <- getOverlappingScores_v2(extended_TSS,data)

Norm$TSS_cov <- NormTSS$pNET
Norm$pausing_density <- Norm$TSS_cov/Norm$normpNET

df.3 <- data.frame("geneID" = Norm$gene_id ,"genebody"=Norm$pNET , "PD" = Norm$pausing_density,"size" = Norm$widthg)
df.3 <- df.3[order(-df.3$genebody),]
nrow(df.3)
nb <- nrow(df.3)
nb.down <- round(0.25*nb)
df.3 <- df.3[0:(nb-nb.down),]
nrow(df.3)

df.3 <- df.3[complete.cases(df.3$PD), ]
df.3 <- df.3[df.3$PD != "Inf", ]
df.3$matchsppRNA <- match(df.3$geneID,all,nomatch="0")
df.3 <- df.3[order(-df.3$PD),]
a <- df.3$matchsppRNA != 0
df.3$matchsppRNA[a] <- 1

nb = 150
ls <- group(df.3, n = nb, method = "greedy")
list <- ls %>% group_split(.groups) 
rec <- lapply(list, `[[`, 5)
recs <- lapply(rec, sum)
recsu <- unlist(recs)
recsn <- recsu/nb

ddf <- data.frame("Pause_Index" = seq(1,length(recsn),by=1), "y" = recsn)

pdf("pausing_sppRNAgenes.pdf")
ggplot(ddf, aes(x=Pause_Index, y=y)) + geom_point() + geom_smooth() + theme_classic() + labs(title="Distribution of sppRNA in relation with 5' Pause Index values",y = "Proportion of sppRNA genes",x=" Number of bins of genes ordered from High to Low by Pausing Index value")

dev.off()