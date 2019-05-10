

coding_genes_n0big <- coding_genes_nO[width(coding_genes_nO) > 1000, ]

coding_genes_n0big <- dropSeqlevels(coding_genes_n0big,c("Mt","Pt"),pruning.mode ="coarse")

TSS <- resize(coding_genes_n0big, width = 1 , fix = "start")
extended_TSS <- extend(TSS,upstream=149,downstream=150)
shrunk_genes <- extend(coding_genes_n0big,upstream=-300,downstream=-300)



pNET <- import('/home/bcm215/groupdirs/SCIENCE-PLEN-Marquardt_lab/The Holy Grail/Arabidopsis/Zhu 2018 (PMID 30374093) - pNET-Seq/pNET-Seq/Norm1M/f/s03_s04_unph_fw_rev_norm1M.bedgraph.gz', format = 'bedgraph')

seqlevels(pnet) <- c("1", "2", "3", "4", "5")

data <- list("pNET" = pnet)
Norm <- getOverlappingScores_v2(shrunk_genes,data)
Norm$normpNET <- Norm$pNET/width(Norm)
NormTSS <- getOverlappingScores_v2(extended_TSS,data)

Norm$TSS_cov <- NormTSS$pNET
Norm$pausing_density <- Norm$TSS_cov/Norm$normpNET

df.3 <- data.frame("genes" = Norm$gene_id ,"genebody"=Norm$pNET , "PD" = Norm$pausing_density)
df.3 <- df.3[order(-df.3$genebody),]
nrow(df.3)
nb <- nrow(df.3)
nb.down <- round(0.25*nb)
df.3 <- df.3[0:(nb-nb.down),]
nrow(df.3)

df.3 <- df.3[complete.cases(df.3$PD), ]
df.3 <- df.3[df.3$PD != "Inf", ]
df.3$matchsppRNA <- match(df.3$genes,all_hen2,nomatch="0")

df.3 <- df.3[order(-df.3$PD),]



df_1 <- df.3[0:200,]
df_2 <- df.3[200:400,]
df_3 <- df.3[400:600,]
df_4 <- df.3[600:800,]
df_5 <- df.3[800:1000,]
df_6 <- df.3[1000:1200,]
df_7 <- df.3[1200:1400,]
df_8 <- df.3[1400:1600,]
df_9 <- df.3[1600:1800,]
df_10 <- df.3[1800:2000,]
df_11 <- df.3[2000:2200,]
df_12 <- df.3[2200:2400,]
df_13 <- df.3[2400:2600,]
df_14 <- df.3[2600:2800,]
df_15 <- df.3[2800:3000,]
df_16 <- df.3[3000:3200,]
df_17 <- df.3[3200:3400,]
df_18 <- df.3[3400:3600,]
df_19 <- df.3[3600:3800,]
df_19b <- df.3[3800:4000,]
df_20 <- df.3[4000:5000,]
df_21 <- df.3[5000:6000,]
df_22 <- df.3[6000:7000,]
df_23 <- df.3[7000:8000,]
df_24 <- df.3[8000:10000,]
df_25 <- df.3[10000:12000,]
df_26 <- df.3[10000:14126,]


n1 <- nrow(df_1[df_1$matchsppRNA >0,])/200
n2 <- nrow(df_2[df_2$matchsppRNA >0,])/200
n3 <- nrow(df_3[df_3$matchsppRNA >0,])/200
n4 <- nrow(df_4[df_4$matchsppRNA >0,])/200
n5 <- nrow(df_5[df_5$matchsppRNA >0,])/200
n6 <- nrow(df_6[df_6$matchsppRNA >0,])/200
n7 <- nrow(df_7[df_7$matchsppRNA >0,])/200
n8 <- nrow(df_8[df_8$matchsppRNA >0,])/200
n9 <- nrow(df_9[df_9$matchsppRNA >0,])/200
n10 <- nrow(df_10[df_10$matchsppRNA >0,])/200
n11 <- nrow(df_11[df_11$matchsppRNA >0,])/200
n12 <- nrow(df_12[df_12$matchsppRNA >0,])/200
n13 <- nrow(df_13[df_13$matchsppRNA >0,])/200
n14 <- nrow(df_14[df_14$matchsppRNA >0,])/200
n15 <- nrow(df_15[df_15$matchsppRNA >0,])/200
n16 <- nrow(df_16[df_16$matchsppRNA >0,])/200
n17 <- nrow(df_17[df_17$matchsppRNA >0,])/200
n18<- nrow(df_18[df_18$matchsppRNA >0,])/200
n19 <- nrow(df_19[df_19$matchsppRNA >0,])/200
n19b <- nrow(df_19b[df_19b$matchsppRNA >0,])/200
n20 <- nrow(df_20[df_20$matchsppRNA >0,])/1000
n21 <- nrow(df_21[df_21$matchsppRNA >0,])/1000
n22 <- nrow(df_22[df_22$matchsppRNA >0,])/1000
n23 <- nrow(df_23[df_23$matchsppRNA >0,])/1000
n24 <- nrow(df_24[df_24$matchsppRNA >0,])/2000
n25 <- nrow(df_25[df_25$matchsppRNA >0,])/2000
n26 <- nrow(df_25[df_25$matchsppRNA >0,])/2126


ddf <- data.frame("Pause_Index" = c(200,400,600,800,1000,1200,1400,1600,1800,2000,2200,2400,2600,2800,3000,3200,3400,3600,3800,4000,5000,6000,7000,8000,10000,12000,14126), "y" = c(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14,n15,n16,n17,n18,n19,n19b,n20,n21,n22,n23,n24,n25,n26))

pdf("pausing_sppRNAgenes.pdf")
ggplot(ddf, aes(x=Pause_Index, y=y)) + geom_point() + geom_smooth() + theme_classic() + labs(title="Distribution of sppRNA in relation with 5' Pause Index values",y = "Proportion of sppRNA genes",x=" Number of genes ordered from High to Low by Pausing Index value")
dev.off()