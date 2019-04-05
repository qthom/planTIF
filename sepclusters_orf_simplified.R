#this script is made to measure the diversity of clusters that overlap or not known genes

genes_nO <- coding_genes_nO

gr$score <- NULL
gr$name  <- as.numeric(gr$name)
gr$cluster_number <- as.character(seq_len(nrow(mcols(gr))))

#finds clusters overlapping coding genes all genes and ncRNA genes

gr$overlap_allgenes <- countOverlaps(gr,genes,type=c("any"),ignore.strand=TRUE)
gr$overlap_allcodinggenes <- countOverlaps(gr,coding_genes,type=c("any"),ignore.strand=TRUE)
gr$overlap_anyncRNA <- countOverlaps(gr,ncRNA,type=c("any"),ignore.strand=TRUE)
gr$overlap_nO <- countOverlaps(gr,genes_nO,type=c("any"),ignore.strand=TRUE)

intragenic_all <- gr[gr$overlap_allgenes != 0]
intragenic_nO <- gr[gr$overlap_nO != 0]
intergenic <- gr[gr$overlap_allcodinggenes == 0]
new_intergenic <- gr[gr$overlap_allgenes == 0 & gr$overlap_anyncRNA == 0 & gr$overlap_nO == 0 & gr$overlap_allcodinggenes == 0]

#identifies clusters overlaping one or several orf in the sense and antisense direction

intragenic_nO$overlap_orf_anti <- countOverlaps(intragenic_nO,genes_nO,type=c("any"),ignore.strand=TRUE)
intragenic_nO$overlap_orf_sens <- countOverlaps(intragenic_nO,genes_nO,type=c("any"),ignore.strand=FALSE)

sens_intragenic <- intragenic_nO[intragenic_nO$overlap_orf_sens != 0]
antisens_intragenic <- intragenic_nO[intragenic_nO$overlap_orf_sens == 0 & intragenic_nO$overlap_orf_anti == 1]

#identifies clusters overlaping one or several orf
overlapping_one <- sens_intragenic[sens_intragenic$overlap_orf_sens == 1 ]
overlapping_two <- sens_intragenic[sens_intragenic$overlap_orf_sens == 2 ]
overlapping_more <- sens_intragenic[sens_intragenic$overlap_orf_sens > 2 ]

#identify where overlapping clusters of one orf locate

overlapping_one$overlap_3 <- countOverlaps(overlapping_one,exTTS,type=c("any"),ignore.strand=FALSE)
overlapping_one$overlap_5 <- countOverlaps(overlapping_one,exTSS,type=c("any"),ignore.strand=FALSE)
overlapping_one$overlap_3orf <- countOverlaps(overlapping_one,To,type=c("any"),ignore.strand=FALSE)
overlapping_one$overlap_5orf <- countOverlaps(overlapping_one,So,type=c("any"),ignore.strand=FALSE)

#defines the clusters overlapping just the start, the end  or within the annotated ORF. And the one covering the orf

over_orf <- overlapping_one[overlapping_one$overlap_5orf == 1 & overlapping_one$overlap_3orf == 1]

Ntrunca4 <- overlapping_one[overlapping_one$overlap_3 == 1 | overlapping_one$overlap_3orf == 1]
Ntrunca4 <- Ntrunca4[Ntrunca4$overlap_5 == 0 & Ntrunca4$overlap_5orf == 0]

Ctrunca4 <- overlapping_one[overlapping_one$overlap_5 == 1 | overlapping_one$overlap_5orf == 1]
Ctrunca4 <- Ctrunca4[Ctrunca4$overlap_3 == 0 & Ctrunca4$overlap_3orf == 0]

within <- overlapping_one[overlapping_one$overlap_3 == 0 & overlapping_one$overlap_5 == 0 & overlapping_one$overlap_3orf == 0 & overlapping_one$overlap_5orf == 0]

l_clusters =length(gr)
l_ovclusters = length(intragenic_nO)
l_ov1clusters = length(overlapping_one)
l_over_orf = (length(over_orf))/l_clusters
l_Ntrunca = (length(Ntrunca4))/l_clusters
l_Ctrunca = (length(Ctrunca4))/l_clusters
l_within = (length(within))/l_clusters
l_antisens = (length(antisens_intragenic))/l_clusters
l_intergenic = (length(intergenic))/l_clusters
l_2ormore = (length(overlapping_two)+length(overlapping_more))/l_clusters


all = (l_over_orf + l_Ntrunca + l_Ctrunca + l_within + l_antisens + l_intergenic + l_2ormore) 

l_over_orf = l_over_orf/all
l_Ntrunca = l_Ntrunca/all
l_Ctrunca = l_Ctrunca/all
l_within = l_within/all
l_antisens = l_antisens/all
l_intergenic = l_intergenic/all
l_2ormore = l_2ormore/all

