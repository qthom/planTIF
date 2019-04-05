#this scrit take genomic ranges and associate annotation


coding_genes_nO$genes_number <- as.character(seq_len(nrow(mcols(coding_genes_nO))))

map <- data.frame("geneID" = coding_genes_nO$names, "gene_number" = coding_genes_nO$genes_number)

ov <- as.data.frame(findOverlaps(coding_genes_nO,gr, ignore.strand = FALSE))

colnames(ov)[1] <- 'gene_number'

m <- merge(ov,map,by='gene_number',all=F)

gr$subjectHits <- as.character(seq_len(nrow(mcols(gr))))

df <- as.data.frame(gr)

m2 <- merge(m,df,by='subjectHits',all=F)

gr <- GRanges(seqnames = m2$seqnames, ranges = IRanges(m2$start, m2$end), strand = m2$strand, geneID = m2$geneID, score = m2$name)
