

df <- data.frame("geneID" = genes$gene_id)
wt_genes <- small_wt$geneID
hen2_genes <- small_hen2$geneID
df$small_wt <- match(df$geneID, wt_genes,nomatch=0)
df$small_hen2 <- match(df$geneID, hen2_genes,nomatch=0)
wtcold_genes <- small_wtcold$geneID
hen2cold_genes <- small_hen2cold$geneID
df$small_wtcold <- match(df$geneID, wtcold_genes,nomatch=0)
df$small_hen2cold <- match(df$geneID, hen2cold_genes,nomatch=0)
ssrp1_genes <- small_ssrp1$geneID
spt16_genes <- small_spt16$geneID
df$small_ssrp1 <- match(df$geneID, ssrp1_genes,nomatch=0)
df$small_spt16 <- match(df$geneID, spt16_genes,nomatch=0)


a <- df$small_wt != 0
df$small_wt[a] <- 1
b <- df$small_hen2 != 0
df$small_hen2[b] <- 1
d <- df$small_wtcold != 0
df$small_wtcold[d] <- 1
e <- df$small_hen2cold != 0
df$small_hen2cold[e] <- 1
f <- df$small_ssrp1 != 0
df$small_ssrp1[f] <- 1
g <- df$small_spt16 != 0
df$small_spt16[g] <- 1

l1 <- df$geneID[df$small_hen2 == 1]
l2 <- df$geneID[df$small_hen2cold == 1]
l3 <- df$geneID[df$small_wt == 1]
l4 <- df$geneID[df$small_wtcold == 1]

l5 <- df$geneID[df$small_hen2 == 1 & df$small_hen2cold == 1 ]
l6 <- df$geneID[df$small_hen2 == 1 & df$small_wt == 1 ]
l7 <- df$geneID[df$small_hen2 == 1 & df$small_wtcold == 1]
l7b <- df$geneID[df$small_hen2cold == 1 & df$small_wt == 1]
l8 <- df$geneID[df$small_hen2cold == 1 & df$small_wtcold == 1]
l9 <- df$geneID[df$small_wt == 1 & df$small_wtcold == 1]

l10 <- df$geneID[df$small_hen2 == 1 & df$small_hen2cold == 1 & df$small_wt == 1 ]
l11 <- df$geneID[df$small_hen2 == 1 & df$small_hen2cold == 1 &  df$small_wtcold == 1]
l12 <- df$geneID[df$small_hen2 == 1 &  df$small_wt == 1 & df$small_wtcold == 1]
l13 <- df$geneID[df$small_hen2cold == 1 & df$small_wt == 1 & df$small_wtcold == 1]

lf <- df$geneID[df$small_hen2 == 1 & df$small_hen2cold == 1 & df$small_wt == 1 & df$small_wtcold == 1]


#plot

ll1 <- length(l1)
ll2 <- length(l2)
ll3 <- length(l3)
ll4 <- length(l4)
ll5 <- length(l5)
ll6 <- length(l6)
ll7 <- length(l7)
ll7b <- length(l7b)
ll8 <- length(l8)
ll9 <- length(l9)
ll10 <- length(l10)
ll11 <- length(l11)
ll12 <- length(l12)
ll13 <- length(l13)
llf <- length(lf)


pdf("Overlapping_sppRNA genes in data sets.pdf",7,8)

draw.quad.venn(ll1, ll2, ll3, ll4, ll5, ll6, ll7, ll7b, ll8,ll9, ll10, ll11, ll12,ll13,llf, lty = "blank", fill = c("skyblue", "pink1","mediumorchid", "orange"), category = c("hen2-2","hen2-2 Cold", "wild tyoe", "wild type Cold"))

dev.off()


#Strict categories

s1 <- df$geneID[df$small_hen2 == 1 & df$small_hen2cold == 0 & df$small_wt == 0 & df$small_wtcold == 0]
s2 <- df$geneID[df$small_hen2 == 0 & df$small_hen2cold == 1 & df$small_wt == 0 & df$small_wtcold == 0]
s3 <- df$geneID[df$small_hen2 == 0 & df$small_hen2cold == 0 & df$small_wt == 1 & df$small_wtcold == 0]
s4 <- df$geneID[df$small_hen2 == 0 & df$small_hen2cold == 0 & df$small_wt == 0 & df$small_wtcold == 1]

s12 <- df$geneID[df$small_hen2 == 1 & df$small_hen2cold == 1 & df$small_wt == 0 & df$small_wtcold == 0 ]
s13 <- df$geneID[df$small_hen2 == 1 & df$small_hen2cold == 0 & df$small_wt == 1 & df$small_wtcold == 0]
s14 <- df$geneID[df$small_hen2 == 1 & df$small_hen2cold == 0 & df$small_wt == 0  & df$small_wtcold == 1]
s23 <- df$geneID[df$small_hen2 == 0 & df$small_hen2cold == 1 & df$small_wt == 1 & df$small_wtcold == 0]
s24 <- df$geneID[df$small_hen2 == 0 & df$small_hen2cold == 1 & df$small_wt == 0 & df$small_wtcold == 1]
s34 <- df$geneID[df$small_hen2 == 0 & df$small_hen2cold == 0 & df$small_wt == 1 & df$small_wtcold == 1]

s123 <- df$geneID[df$small_hen2 == 1 & df$small_hen2cold == 1 & df$small_wt == 1 & df$small_wtcold == 0]
s124 <- df$geneID[df$small_hen2 == 1 & df$small_hen2cold == 1 & df$small_wt == 0 &  df$small_wtcold == 1]
s134 <- df$geneID[df$small_hen2 == 1 & df$small_hen2cold == 0 &  df$small_wt == 1 & df$small_wtcold == 1]
s234 <- df$geneID[df$small_hen2 == 0 & df$small_hen2cold == 1 & df$small_wt == 1 & df$small_wtcold == 1]

s1234 <- df$geneID[df$small_hen2 == 1 & df$small_hen2cold == 1 & df$small_wt == 1 & df$small_wtcold == 1]

df$type <- "0"
df$type[s1] <- "s1"
df$type[s2] <- "s2"
df$type[s3] <- "s3"
df$type[s4] <- "s4"
df$type[s12] <- "s12"
df$type[s13] <- "s13"
df$type[s14] <- "s14"
df$type[s23] <- "s23"
df$type[s24] <- "s24"
df$type[s34] <- "s34"
df$type[s123] <- "s123"
df$type[s124] <- "s124"
df$type[s134] <- "s134"
df$type[s234] <- "s234"
df$type[s1234] <- "s1234"

library(tibble)

tair_ann <- import.gff3("Arabidopsis_thaliana.TAIR10.26.gff3")
mapping <- data.frame("geneID" = sub("gene:", "", tair_ann$ID), "name"=tair_ann$external_name)
tmp <- left_join(x=tibble(geneID=df$geneID), y=mapping, by="geneID", na_matches = "never")
tmp <- as.data.frame(tmp)
df$geneName <- tmp$name

df <- df[(df$type != "0"),]

write.table(df, file="filtered_sppRNAs_categories_VennDiagram.csv", quote=F, sep="\t", row.names=F, col.names=T)
