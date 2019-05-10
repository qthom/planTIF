Norm$match <- match(Norm$gene_id,all, nomatch = 0)

pnet_sppRNA <- Norm[Norm$match > 0,]
pnet_all <- Norm[Norm$match == 0,]

dat1 <- data.frame("sppRNA"="sppRNA", "pnet"=pnet_sppRNA$normpNET)
dat2 <- data.frame("sppRNA"="all_no_sppRNA", "pnet"=pnet_all$normpNET)

dat <- rbind(dat1,dat2)

pdf("violin_expression_sppRNA.pdf",6,6)
ggplot(dat, aes(x=sppRNA, y=log((10000*pnet),10), fill=sppRNA)) + geom_violin() + geom_boxplot(width=0.2) +theme_classic() +  labs(title="Distribution of pNET-Seq between sppRNA genes and genes without sppRNA",y = "log(TPM) at TSS",x="1639 genes")
dev.off()