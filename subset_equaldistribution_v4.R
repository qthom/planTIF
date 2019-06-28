#this script finds a subset of genes with the same expression distribution betwween hen2 small transcript genes and col0 genes without any small transcripts

#create data frame with gene's pNET-Seq expression

#at the end of the script the files of sppRNA per samples are created

pNET <- import('/home/bcm215/groupdirs/SCIENCE-PLEN-Marquardt_lab/The Holy Grail/Arabidopsis/Zhu 2018 (PMID 30374093) - pNET-Seq/pNET-Seq/Norm1M/f/s03_s04_unph_fw_rev_norm1M.bedgraph.gz', format = 'bedgraph')

pnet <- dropSeqlevels(pNET,c("Pt","Mt"),pruning.mode ="coarse")

seqlevels(pnet) <- c("1", "2", "3", "4", "5")

data <- list("pNET" = pnet)
Norm <- getOverlappingScores_v2(shrinked_genes,data)
txdb2 <- dropSeqlevels(txdb,c("6","7"),pruning.mode ="coarse")
Norm$normpNET <- Norm$pNET/width(Norm)


#download the sppRNA clusters from each data set

seqlevels(coding_genes_nO) <- seqlevels(gr_wt)

#Subsets the sppRNA from the total cluster files

gr <- gr_wt
source("sepclusters_forsmallclusters.R")
small_wt <- small

gr <- gr_hen2
source("sepclusters_forsmallclusters.R")
small_hen2 <- small

gr <- gr_wtcold
source("sepclusters_forsmallclusters.R")
small_wtcold <- small

gr <- gr_hen2cold
source("sepclusters_forsmallclusters.R")
small_hen2cold <- small

gr <- gr_ssrp1
source("sepclusters_forsmallclusters.R")
small_ssrp1 <- small

gr <- gr_spt16
source("sepclusters_forsmallclusters.R")
small_spt16 <- small

#Annotate the sppRNAs

gr <- small_hen2
source("annotate.R")
small_hen2 <- gr

gr <- small_wt
source("annotate.R")
small_wt <- gr

gr <- small_hen2cold
source("annotate.R")
small_hen2cold <- gr

gr <- small_wtcold
source("annotate.R")
small_wtcold <- gr

gr <- small_ssrp1
source("annotate.R")
small_ssrp1 <- gr

gr <- small_spt16
source("annotate.R")
small_spt16 <- gr

#select for "long" genes > 1000 bp for finding the control genes set

genes <- coding_genes[width(coding_genes) > 1000]

#finds genes with sppRNA for each data set
#22C

df <- data.frame("geneID" = genes$gene_id)
wt_genes <- small_wt$geneID
hen2_genes <- small_hen2$geneID
df$small_wt <- match(df$geneID, wt_genes,nomatch=0)
df$small_hen2 <- match(df$geneID, hen2_genes,nomatch=0)

#cold

wtcold_genes <- small_wtcold$geneID
hen2cold_genes <- small_hen2cold$geneID
df$small_wtcold <- match(df$geneID, wtcold_genes,nomatch=0)
df$small_hen2cold <- match(df$geneID, hen2cold_genes,nomatch=0)

#FACT

spt16_genes <- small_spt16$geneID
ssrp1_genes <- small_ssrp1$geneID
df$small_spt16 <- match(df$geneID, spt16_genes,nomatch=0)
df$small_ssrp1 <- match(df$geneID, ssrp1_genes,nomatch=0)

#finds genes with no sppRNA and genes with and assign categories

a <- df$small_wt != 0
df$small_wt[a] <- 1
b <- df$small_hen2 != 0
df$small_hen2[b] <- 1
d <- df$small_wtcold != 0
df$small_wtcold[d] <- 1
e <- df$small_hen2cold != 0
df$small_hen2cold[e] <- 1
f <- df$small_spt16 != 0
df$small_spt16[f] <- 1
g <- df$small_ssrp1 != 0
df$small_ssrp1[g] <- 1

all_hen2 <- df$geneID[df$small_hen2 == 1 | df$small_hen2cold == 1 & df$small_wt == 0 & df$small_wtcold == 0 & df$small_ssrp1 == 0 & df$small_spt16 == 0]
all_wt <- df$geneID[df$small_hen2 == 0 & df$small_hen2cold == 0 & df$small_wt == 1 | df$small_wtcold == 1]
all_warm <- df$geneID[df$small_hen2 == 1 | df$small_wt == 1 | df$small_ssrp1 == 1 | df$small_spt16 == 1]
all <- df$geneID[df$small_hen2 == 1 | df$small_hen2cold == 1 | df$small_wt == 1 | df$small_wtcold == 1 | df$small_ssrp1 == 1 | df$small_spt16 == 1]
no_small <- df$geneID[df$small_hen2 == 0 & df$small_wt == 0 & df$small_hen2cold == 0 & df$small_wtcold == 0 & df$small_ssrp1 == 0 & df$small_spt16 == 0]

#Provide Names to geneID of sppRNA genes and without sppRNAs

tair_ann <- import.gff3("/home/bcm215/groupdirs/SCIENCE-PLEN-Marquardt_lab/Quentin/sequencing_files/TAIR10/Arabidopsis_thaliana.TAIR10.26.gff3")
mapping <- data.frame("geneID" = sub("gene:", "", tair_ann$ID), "name"=tair_ann$external_name)
tmp <- left_join(x=tibble(geneID=df$geneID), y=mapping, by="geneID", na_matches = "never")
tmp <- as.data.frame(tmp)
df$geneName <- tmp$name

df$type <- "0"
df$type[all_hen2] <- "hen2"
df$type[no_small] <- "None"

MyData <- data.frame('geneID' = Norm$gene_id, 'V3' = Norm$normpNET)
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
TSS.both2 <- genes[both2$geneID]
TSS.both2 <- resize(TSS.both2, width = 1 , fix = "start")
TSS.ctrlboth2 <- genes[ctrlboth2$geneID]
TSS.ctrlboth2 <- resize(TSS.ctrlboth2, width = 1 , fix = "start")


grl_pnet <- list("B" = TSS.both2, "Control of B" = TSS.ctrlboth2)


#Termination site for sppRNAs

small_all_hen2 <- do.call("c",(list(small_hen2,small_hen2cold)))
small_all_hen2$match <- match(small_all_hen2$geneID,TSS.both2$gene_id,nomatch = 0)
small_all_hen2 <- small_all_hen2[small_all_hen2$match > 0,]
small_all_hen2$dup <- duplicated(small_all_hen2$geneID)
small_all_hen2 <- small_all_hen2[small_all_hen2$dup == FALSE,]

TTS_hen2 <- resize(small_all_hen2, width = 1, fix = "end")
grl_PAS <- list("PAS of sppRNA in Hen2" = TTS_hen2)


#End of genes

genes <- genes(txdb)

genes$match <- match(genes$gene_id,TSS.both2$gene_id,nomatch=0)
genes_sppRNA <- genes[genes$match > 0,]
genes_sppRNA_PAS <- resize(genes_sppRNA, width = 1, fix = "end")
genes$matchnospp <- match(genes$gene_id,TSS.ctrlboth2$gene_id,nomatch=0)
genes_nosppRNA <- genes[genes$matchnospp > 0,]
genes_nosppRNA_PAS <- resize(genes_nosppRNA, width = 1, fix = "end")

grl_pnetend <- list("B" = genes_sppRNA_PAS, "Control of B" = genes_nosppRNA_PAS)



filename = "hen2_sppRNA.csv"
rite.table(as.data.frame(small_hen2), file=filename, quote=F, sep="\t", row.names=F, col.names=T)

filename = "wt_sppRNA.csv"
write.table(as.data.frame(small_wt), file=filename, quote=F, sep="\t", row.names=F, col.names=T)

filename = "hen2cold_sppRNA.csv"
write.table(as.data.frame(small_hen2cold), file=filename, quote=F, sep="\t", row.names=F, col.names=T)

filename = "spt16_sppRNA.csv"
write.table(as.data.frame(small_spt16), file=filename, quote=F, sep="\t", row.names=F, col.names=T)

filename = "ssrp1_sppRNA.csv"
write.table(as.data.frame(small_ssrp1), file=filename, quote=F, sep="\t", row.names=F, col.names=T)

filename = "wtcold_sppRNA.csv"
write.table(as.data.frame(small_wtcold), file=filename, quote=F, sep="\t", row.names=F, col.names=T)
