#This script calculates the proper WT control for each PAS-Seq dataset.
#mut refers to the mutant dataset (bedgraph)
#WT to the wild-type dataset (bedgraph)

mut <- import("file.bedgraph")
WT <- import("file.bedgraph")

mut <- dropSeqlevels(mut,c("Pt","Mt"),pruning.mode ="coarse")
mut <- resize(mut, width = 1 , fix = "center")
seqlevels(mut) <- c("1", "2", "3", "4", "5")
WT <- dropSeqlevels(WT,c("Pt","Mt"),pruning.mode ="coarse")
WT <- resize(WT, width = 1 , fix = "center")
seqlevels(WT) <- c("1", "2", "3", "4", "5")


TTS.2 <- resize(coding_genes_n0big, width = 1 , fix = "end")
TTS.e <- extend(TTS.2,upstream=25,downstream=25)
seqinfo(mut) <- seqinfo(TTS.e)
seqinfo(WT) <- seqinfo(TTS.e)
data <- list("mut" = mut,"WT"= WT)
Norm <- getOverlappingScores_v2(TTS.e,data)
names(Norm) <- Norm$gene_id
#Provide Names to geneID of sppRNA genes and without sppRNAs

df$type <- "0"
df$type[all] <- "hen2"
df$type[no_small] <- "None"

MyData <- data.frame('geneID' = Norm$gene_id, 'V3' = Norm$WT)
df2 <- merge(df, MyData, by='geneID', all=T)

#Calculate the control set for TSS of sppRNA vs control set of genes

df2 <- df2[complete.cases(df2$type,df2$V3), ]
sorteddf2.2 <- df2[order(df2$V3),]
nb <- nrow(df2)
nb.down <- round(0.20*nb)
nb.up <- round(0.10*nb)
df2.2 <- sorteddf2.2[nb.down:(nb-nb.up),]
none2 <- data.frame("geneID" = df2.2$geneID[df2.2$type == "None"], "exp" = df2.2$V3[df2.2$type == "None"])
both2 <- data.frame("geneID" = df2.2$geneID[df2.2$type == "hen2"], "exp" = df2.2$V3[df2.2$type == "hen2"])
none2 <- none2[complete.cases(none2$exp), ]
rowgenes_both2 <- findMatchedControl_v2(both2$exp,none2$exp)
ctrlboth2 <- none2[rowgenes_both2,]


ctrl <- Norm[ctrlboth2$geneID,]
sppRNA <- Norm[both2$geneID,]