#in this script the NET-Seq signal for sppRNA genes promoters and their control genes without sppRNA was plotted and the significance calculated between the two population
#with the wilcoxon test

pNET <- import('/home/bcm215/groupdirs/SCIENCE-PLEN-Marquardt_lab/The Holy Grail/Arabidopsis/Zhu 2018 (PMID 30374093) - pNET-Seq/pNET-Seq/Norm1M/f/s03_s04_unph_fw_rev_norm1M.bedgraph.gz', format = 'bedgraph')

pnet <- dropSeqlevels(pNET,c("Pt","Mt"),pruning.mode ="coarse")

seqlevels(pnet) <- c("1", "2", "3", "4", "5")

data <- list("pNET" = pnet)



sppRNA.ex <- extend(TSS.both2,upstream=-50,downstream=200)
nosppRNA.ex <- extend(TSS.ctrlboth2,upstream=-50,downstream=200)

Norm.1 <- getOverlappingScores_v2(sppRNA.ex,data)
Norm.2 <- getOverlappingScores_v2(nosppRNA.ex,data)
dat1 <- data.frame("sppRNA"="sppRNA", "pNET" = Norm.1$pNET)
dat2 <- data.frame("sppRNA"="nosppRNA", "pNET" = Norm.2$pNET)
dat <- rbind(dat1,dat2)
pdf("violin_significance_metagene.pdf",6,6)
ggplot(dat, aes(x=sppRNA, y=pNET, fill=sppRNA)) + geom_violin() + geom_boxplot(width=0.2) +theme_classic() + scale_fill_manual(values=c("#0072B2","darkorange2")) +  stat_compare_means(method = "wilcox.test") +  labs(title="Distribution of pNET-Seq between sppRNA genes and genes without sppRNA",y = "TPM at TSS",x="1155 genes")
dev.off()
