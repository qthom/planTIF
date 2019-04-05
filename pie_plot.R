
#this summarize the cluster position distribution over genes and calculate the Fold Change between the datasets and plot it

genes <- genes(txdb)

gr <- gr_hen2
source("sepclusters_orf_simplified.R")
hen2_list <- c(l_over_orf,l_Ntrunca,l_Ctrunca,l_within,l_antisens,l_intergenic,l_2ormore)
gr <- gr_hen2
source("sepclusters_overlapping_orf.R")
hen2orf_list <- c(l_o1,l_o2,l_o3,l_o4,l_o5,l_o6,l_o7,l_o8,l_o9)

gr <- gr_wt
source("sepclusters_orf_simplified.R")
Col0_list <- c(l_over_orf,l_Ntrunca,l_Ctrunca,l_within,l_antisens,l_intergenic,l_2ormore)
gr <- gr_wt
source("sepclusters_overlapping_orf.R")
Col0orf_list <- c(l_o1,l_o2,l_o3,l_o4,l_o5,l_o6,l_o7,l_o8,l_o9)

gr <- gr_hen2cold
source("sepclusters_orf_simplified.R")
hen2cold_list <- c(l_over_orf,l_Ntrunca,l_Ctrunca,l_within,l_antisens,l_intergenic,l_2ormore)
gr <- gr_hen2cold
source("sepclusters_overlapping_orf.R")
hen2coldorf_list <- c(l_o1,l_o2,l_o3,l_o4,l_o5,l_o6,l_o7,l_o8,l_o9)

gr <- gr_wtcold
source("sepclusters_orf_simplified.R")
Col0cold_list <- c(l_over_orf,l_Ntrunca,l_Ctrunca,l_within,l_antisens,l_intergenic,l_2ormore)
gr <- gr_wtcold
source("sepclusters_overlapping_orf.R")
Col0coldorf_list <- c(l_o1,l_o2,l_o3,l_o4,l_o5,l_o6,l_o7,l_o8,l_o9)

gr <- gr_ssrp1
source("sepclusters_orf_simplified.R")
ssrp1_list <- c(l_over_orf,l_Ntrunca,l_Ctrunca,l_within,l_antisens,l_intergenic,l_2ormore)
gr <- gr_ssrp1
source("sepclusters_overlapping_orf.R")
ssrp1orf_list <- c(l_o1,l_o2,l_o3,l_o4,l_o5,l_o6,l_o7,l_o8,l_o9)

gr <- gr_spt16
source("sepclusters_orf_simplified.R")
spt16_list <- c(l_over_orf,l_Ntrunca,l_Ctrunca,l_within,l_antisens,l_intergenic,l_2ormore)
gr <- gr_spt16
source("sepclusters_overlapping_orf.R")
spt16orf_list <- c(l_o1,l_o2,l_o3,l_o4,l_o5,l_o6,l_o7,l_o8,l_o9)


df <- data.frame(group = c("Overlapping whole ORF", "Internal Start","Internal End", "Internal Start and End","Antisense","Intergenic","Overlapping 2 or more"), Col0 = Col0_list, hen2 = hen2_list, Col0_cold = Col0cold_list, Hen2_cold = hen2cold_list, SSRP1 = ssrp1_list, SPT16 =spt16_list )
rownames(df) <- df$group
dforf <- data.frame(group = c("Exact TU", "Within 3'UTR only","Extended 3' UTR", "Within both UTRs","Extended 3' and 5' UTR","Within 3' UTR","extended 5' UTR","Over TSS and within 3utr","Over PAS and within 5utr"), Col0 = Col0orf_list, hen2 = hen2orf_list, Col0_cold = Col0coldorf_list, Hen2_cold = hen2coldorf_list, SSRP1 = ssrp1orf_list, SPT16 =spt16orf_list )
rownames(dforf) <- dforf$group


df$log2FC_COL0_Hen2 <- log(df$hen2,2) - log(df$Col0,2)
df$log2FC_COL0cold_Hen2cold <- log(df$Hen2_cold,2) - log(df$Col0_cold,2)
df$log2FC_COL0_SSRP1 <- log(df$SSRP1,2) - log(df$Col0,2)
df$log2FC_COL0_SPT16 <- log(df$SPT16,2) - log(df$Col0,2)

dforf$log2FC_COL0_Hen2 <- log(dforf$hen2,2) - log(dforf$Col0,2)
dforf$log2FC_COL0cold_Hen2cold <- log(dforf$Hen2_cold,2) - log(dforf$Col0_cold,2)
dforf$log2FC_COL0_SSRP1 <- log(dforf$SSRP1,2) - log(dforf$Col0,2)
dforf$log2FC_COL0_SPT16 <- log(dforf$SPT16,2) - log(dforf$Col0,2)

write.table(df,file="FoldChange_TUs.csv",row.names=TRUE, col.names=TRUE)

