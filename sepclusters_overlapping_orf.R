#this script is made to measure the diversity of clusters that overlap or not known genes

genes_nO <- coding_genes_nO
genes <- genes(txdb)
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

overlapping_one$overlap_3orf <- countOverlaps(overlapping_one,To,type=c("any"),ignore.strand=FALSE)
overlapping_one$overlap_5orf <- countOverlaps(overlapping_one,So,type=c("any"),ignore.strand=FALSE)


#defines the clusters overlapping just the start, the end  or within the annotated ORF. And the one covering the orf
over_orf_all <- overlapping_one[overlapping_one$overlap_3orf == 1 & overlapping_one$overlap_5orf == 1]

over_orf_allTSS <- resize(over_orf_all, width = 1 , fix = "start")
over_orf_allTTS <- resize(over_orf_all, width = 1 , fix = "end")

over_orf_all$overlap_3 <- countOverlaps(over_orf_all,exTTS,type=c("any"),ignore.strand=FALSE)
over_orf_all$overlap_5 <- countOverlaps(over_orf_all,exTSS,type=c("any"),ignore.strand=FALSE)
over_orf_all$overlap_TTS <- countOverlaps(over_orf_allTTS,exTTS,type=c("any"),ignore.strand=FALSE)
over_orf_all$overlap_TSS <- countOverlaps(over_orf_allTSS,exTSS,type=c("any"),ignore.strand=FALSE)
over_orf_all$overlap_3orf <- countOverlaps(over_orf_all,To,type=c("any"),ignore.strand=FALSE)
over_orf_all$overlap_5orf <- countOverlaps(over_orf_all,So,type=c("any"),ignore.strand=FALSE)

#exact overlapping
o1 <- over_orf_all[over_orf_all$overlap_3 == 1 & over_orf_all$overlap_3orf == 1 & over_orf_all$overlap_5 == 1 & over_orf_all$overlap_5orf == 1 &  over_orf_all$overlap_TSS == 1 &  over_orf_all$overlap_TTS == 1 ]
#Terminating within 3'UTR before extended PAS and at TSS
o2 <- over_orf_all[over_orf_all$overlap_3 == 0 & over_orf_all$overlap_3orf == 1 & over_orf_all$overlap_5 == 1 & over_orf_all$overlap_5orf == 1 &  over_orf_all$overlap_TSS == 1 &  over_orf_all$overlap_TTS == 0 ]
#extended 3 UTR - over the PAS - at TSS
o3 <- over_orf_all[over_orf_all$overlap_3 == 1 & over_orf_all$overlap_3orf == 1 & over_orf_all$overlap_5 == 1 & over_orf_all$overlap_5orf == 1 &  over_orf_all$overlap_TSS == 1 &  over_orf_all$overlap_TTS == 0 ]
#within both UTRs strictly
o4 <- over_orf_all[over_orf_all$overlap_3 == 0 & over_orf_all$overlap_3orf == 1 & over_orf_all$overlap_5 == 0 & over_orf_all$overlap_5orf == 1 &  over_orf_all$overlap_TSS == 0 &  over_orf_all$overlap_TTS == 0 ]
#over both TSS and PAS
o5 <- over_orf_all[over_orf_all$overlap_3 == 1 & over_orf_all$overlap_3orf == 1 & over_orf_all$overlap_5 == 1 & over_orf_all$overlap_5orf == 1 &  over_orf_all$overlap_TSS == 0 &  over_orf_all$overlap_TTS == 0 ]
#Within 5' UTR strictly and on 3'UTR
o6 <- over_orf_all[over_orf_all$overlap_3 == 1 & over_orf_all$overlap_3orf == 1 & over_orf_all$overlap_5 == 0 & over_orf_all$overlap_5orf == 1 &  over_orf_all$overlap_TSS == 0 &  over_orf_all$overlap_TTS == 1 ]
#Over the 5'UTR and at the 3' UTR
o7 <- over_orf_all[over_orf_all$overlap_3 == 1 & over_orf_all$overlap_3orf == 1 & over_orf_all$overlap_5 == 1 & over_orf_all$overlap_5orf == 1 &  over_orf_all$overlap_TSS == 0 &  over_orf_all$overlap_TTS == 1 ]
#within 3' UTR and at TSS
o8 <- over_orf_all[over_orf_all$overlap_3 == 0 & over_orf_all$overlap_3orf == 1 & over_orf_all$overlap_5 == 1 & over_orf_all$overlap_5orf == 1 &  over_orf_all$overlap_TSS == 0 &  over_orf_all$overlap_TTS == 0 ]
#over PAS and within 5'UTR
o9 <- over_orf_all[over_orf_all$overlap_3 == 1 & over_orf_all$overlap_3orf == 1 & over_orf_all$overlap_5 == 0 & over_orf_all$overlap_5orf == 1 &  over_orf_all$overlap_TSS == 0 &  over_orf_all$overlap_TTS == 0 ]


l_o1 = (length(o1))
l_o2 = (length(o2))
l_o3 = (length(o3))
l_o4 = (length(o4))
l_o5 = (length(o5))
l_o6 = (length(o6))
l_o7 = (length(o7))
l_o8 = (length(o8))
l_o9 = (length(o9))

size <- sum(l_o1,l_o2,l_o3,l_o4,l_o5,l_o6,l_o7,l_o8,l_o9)

l_o1 = (length(o1))/size
l_o2 = (length(o2))/size
l_o3 = (length(o3))/size
l_o4 = (length(o4))/size
l_o5 = (length(o5))/size
l_o6 = (length(o6))/size
l_o7 = (length(o7))/size
l_o8 = (length(o8))/size
l_o9 = (length(o9))/size