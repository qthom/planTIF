gr$score <- NULL
gr$wi <- width(gr)
gr <- gr[gr$wi < 350]

#this get the intergenic and intragenic genomic ranges

gr$overlap_any <- countOverlaps(gr,coding_genes_nO,type=c("any"),ignore.strand=TRUE)
intragenic_any <- gr[gr$overlap_any != 0]

#identifies clusters overlaping one or several orf

intragenic_any$overlap_orf_sens <- countOverlaps(intragenic_any,coding_genes_nO,type=c("any"),ignore.strand=FALSE)
sens_intragenic <- intragenic_any[intragenic_any$overlap_orf_sens != 0]
antisens_intragenic <- intragenic_any[intragenic_any$overlap_orf_sens == 0 & intragenic_any$overlap_orf_anti != 0]

#identifies clusters overlaping one or several orf
overlapping_one <- sens_intragenic[sens_intragenic$overlap_orf_sens == 1 ]

#identify where overlapping clusters of one orf locate

overlapping_one$overlap_3 <- countOverlaps(overlapping_one,TTS,type=c("any"),ignore.strand=FALSE)
overlapping_one$overlap_5ext <- countOverlaps(overlapping_one,extended_TSS,type=c("any"),ignore.strand=FALSE)
overlapping_one$overlap_5 <- countOverlaps(overlapping_one,TSS,type=c("any"),ignore.strand=FALSE)
overlapping_one$overlap_3orf <- countOverlaps(overlapping_one,To,type=c("any"),ignore.strand=FALSE)
overlapping_one$overlap_5orf <- countOverlaps(overlapping_one,So,type=c("any"),ignore.strand=FALSE)


small <- overlapping_one[overlapping_one$overlap_5ext == 1 | overlapping_one$overlap_5orf == 1 | overlapping_one$overlap_5 == 1]
small <- small[small$overlap_3 == 0 & small$overlap_3orf == 0 ]
