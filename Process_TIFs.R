#load the clusters

gr_wt <- import('/home/bcm215/groupdirs/SCIENCE-PLEN-Marquardt_lab/Quentin/TIF-Seq/Data_processing/Quentin_secondBAM/EndTodEnd_align_polyXfiltered/clusters/commonTSS/COL0.Aligned.sortedByCoord.out.dedup.uniq.merged.polyN.commonTSS.bed')
gr_hen2 <- import('/home/bcm215/groupdirs/SCIENCE-PLEN-Marquardt_lab/Quentin/TIF-Seq/Data_processing/Quentin_secondBAM/EndTodEnd_align_polyXfiltered/clusters/commonTSS/HEN2-2.Aligned.sortedByCoord.out.dedup.uniq.merged.polyN.commonTSS.bed')
gr_hen2cold <- import('/home/bcm215/groupdirs/SCIENCE-PLEN-Marquardt_lab/Quentin/TIF-Seq/Data_processing/Quentin_secondBAM/EndTodEnd_align_polyXfiltered/clusters/commonTSS/HEN2-2cold.Aligned.sortedByCoord.out.dedup.uniq.merged.polyN.commonTSS.bed')
gr_wtcold <- import('/home/bcm215/groupdirs/SCIENCE-PLEN-Marquardt_lab/Quentin/TIF-Seq/Data_processing/Quentin_secondBAM/EndTodEnd_align_polyXfiltered/clusters/commonTSS/COL0cold.Aligned.sortedByCoord.out.dedup.uniq.merged.polyN.commonTSS.bed')
gr_spt16 <- import('/home/bcm215/groupdirs/SCIENCE-PLEN-Marquardt_lab/Quentin/TIF-Seq/Data_processing/Quentin_secondBAM/EndTodEnd_align_polyXfiltered/clusters/commonTSS/SPT16.Aligned.sortedByCoord.out.dedup.uniq.merged.polyN.commonTSS.bed')
gr_ssrp1 <- import('/home/bcm215/groupdirs/SCIENCE-PLEN-Marquardt_lab/Quentin/TIF-Seq/Data_processing/Quentin_secondBAM/EndTodEnd_align_polyXfiltered/clusters/commonTSS/SSRP1.Aligned.sortedByCoord.out.dedup.uniq.merged.polyN.commonTSS.bed')

#import the ENSEMBL database

mart = useMart(biomart="plants_mart", host="plants.ensembl.org", dataset="athaliana_eg_gene")
listg <- select(mart,keys=c("protein_coding"),columns=c("chromosome_name","start_position","end_position","strand","gene_biotype","ensembl_gene_id"), keytype="biotype")
lista <- select(mart,keys=c("protein_coding","antisense_RNA","ncRNA","lncRNA","snoRNA"),c("chromosome_name","start_position","end_position","strand","gene_biotype","ensembl_gene_id"), keytype="biotype")
lncRNA <- select(mart,keys=c("antisense_RNA","ncRNA","snRNA","rRNA","tRNA","lncRNA","snoRNA","miRNA"),c("chromosome_name","start_position","end_position","strand","gene_biotype","ensembl_gene_id"), keytype="biotype")

#Find TAIR10 coding genes
genes <- genes(txdb)
keys <- genes$gene_id
genes$txtype <- select(txdb,keys=keys,keytype="GENEID",columns="TXTYPE")$TXTYPE
coding_genes <- genes[genes$txtype == "protein_coding",]

#all genes including ncRNAs

all_genes <- GRanges(names = lista$ensembl_gene_id,seqnames = lista$chromosome_name, ranges = IRanges(lista$start_position, lista$end_position), strand = lista$strand)
names(all_genes) <- lista$ensembl_gene_id
all_genes$names <- lista$ensembl_gene_id
all_genes$biotype <- lista$gene_biotype

#Find ENSEMBL non coding genes

ncRNA <- GRanges(names = lncRNA$ensembl_gene_id,seqnames = lncRNA$chromosome_name, ranges = IRanges( lncRNA$start_position,  lncRNA$end_position), strand =  lncRNA$strand)
names(ncRNA) <- lncRNA$ensembl_gene_id
ncRNA$names <- lncRNA$ensembl_gene_id
ncRNA$biotype <- lncRNA$gene_biotype

#Find non-overlapping genes for the downstream analysis

coding_genes$overlapping <- countOverlaps(coding_genes,coding_genes,type=c("any"),ignore.strand=TRUE)
coding_genes_nO <- coding_genes[coding_genes$overlapping == 1]
coding_genes_nO <- dropSeqlevels(coding_genes_nO,c("Mt","Pt"),pruning.mode ="coarse")

#Gets the TAIR10 annotation of genes

TSS <- resize(coding_genes_nO, width = 1 , fix = "start")
TTS <- resize(coding_genes_nO, width = 1, fix = "end")

exons <- cdsBy(txdb,"gene")
orf <- unlist(range(exons))

five <- fiveUTRsByTranscript(txdb)
futr <- unlist(range(five))
futr <- dropSeqlevels(futr,c("Mt","Pt"),pruning.mode ="coarse")
three <- threeUTRsByTranscript(txdb)
tutr <- unlist(range(three))
tutr <- dropSeqlevels(tutr,c("Mt","Pt"),pruning.mode ="coarse")

orf$overlapping_orf <- countOverlaps(orf,coding_genes_nO,type=c("any"),ignore.strand=FALSE)
orfs <- orf[orf$overlapping_orf == 1]

So <- resize(orfs, width = 1 , fix = "start")
To <- resize(orfs, width = 1, fix = "end")

extend <- function(x, upstream=0, downstream=0)     
{
    if (any(strand(x) == "*"))
        warning("'*' ranges were treated as '+'")
    on_plus <- strand(x) == "+" | strand(x) == "*"
    new_start <- start(x) - ifelse(on_plus, upstream, downstream)
    new_end <- end(x) + ifelse(on_plus, downstream, upstream)
    ranges(x) <- IRanges(new_start, new_end)
    trim(x)
}

extended_TSS <- extend(TSS,upstream=40,downstream=40)

coding_genes_n0big <- coding_genes_nO[width(coding_genes_nO) > 500, ]
shrinked_genes <- extend(coding_genes_n0big,upstream=-200,downstream=-200)
shrinked_genes <- dropSeqlevels(shrinked_genes,c("Mt","Pt"),pruning.mode ="coarse")

exTSS <- extend(TSS, upstream=10,downstream=10)
exTTS <- extend(TTS, upstream=10,downstream=10)


gr_wt_PN <- import('/home/bcm215/groupdirs/SCIENCE-PLEN-Marquardt_lab/Quentin/TIF-Seq/Data_processing/Quentin_secondBAM/EndTodEnd_align_polyXfiltered/clusters/polyN/COL0.Aligned.sortedByCoord.out.dedup.uniq.merged.polyN.bed')
gr_hen2_PN <- import('/home/bcm215/groupdirs/SCIENCE-PLEN-Marquardt_lab/Quentin/TIF-Seq/Data_processing/Quentin_secondBAM/EndTodEnd_align_polyXfiltered/clusters/polyN/HEN2-2.Aligned.sortedByCoord.out.dedup.uniq.merged.polyN.bed')
gr_hen2cold_PN <- import('/home/bcm215/groupdirs/SCIENCE-PLEN-Marquardt_lab/Quentin/TIF-Seq/Data_processing/Quentin_secondBAM/EndTodEnd_align_polyXfiltered/clusters/polyN/HEN2-2cold.Aligned.sortedByCoord.out.dedup.uniq.merged.polyN.bed')
gr_wtcold_PN <- import('/home/bcm215/groupdirs/SCIENCE-PLEN-Marquardt_lab/Quentin/TIF-Seq/Data_processing/Quentin_secondBAM/EndTodEnd_align_polyXfiltered/clusters/polyN/COL0cold.Aligned.sortedByCoord.out.dedup.uniq.merged.polyN.bed')
gr_spt16_PN <- import('/home/bcm215/groupdirs/SCIENCE-PLEN-Marquardt_lab/Quentin/TIF-Seq/Data_processing/Quentin_secondBAM/EndTodEnd_align_polyXfiltered/clusters/polyN/SPT16.Aligned.sortedByCoord.out.dedup.uniq.merged.polyN.bed')
gr_ssrp1_PN <- import('/home/bcm215/groupdirs/SCIENCE-PLEN-Marquardt_lab/Quentin/TIF-Seq/Data_processing/Quentin_secondBAM/EndTodEnd_align_polyXfiltered/clusters/polyN/SSRP1.Aligned.sortedByCoord.out.dedup.uniq.merged.polyN.bed')


gr_wt_r <- import('/home/bcm215/groupdirs/SCIENCE-PLEN-Marquardt_lab/Quentin/TIF-Seq/Data_processing/Quentin_secondBAM/EndTodEnd_align_polyXfiltered/clusters/raw/COL0.Aligned.sortedByCoord.out.dedup.uniq.merged.bed')
gr_hen2_r <- import('/home/bcm215/groupdirs/SCIENCE-PLEN-Marquardt_lab/Quentin/TIF-Seq/Data_processing/Quentin_secondBAM/EndTodEnd_align_polyXfiltered/clusters/raw/HEN2-2.Aligned.sortedByCoord.out.dedup.uniq.merged.bed')
gr_hen2cold_r <- import('/home/bcm215/groupdirs/SCIENCE-PLEN-Marquardt_lab/Quentin/TIF-Seq/Data_processing/Quentin_secondBAM/EndTodEnd_align_polyXfiltered/clusters/raw/HEN2-2cold.Aligned.sortedByCoord.out.dedup.uniq.merged.bed')
gr_wtcold_r <- import('/home/bcm215/groupdirs/SCIENCE-PLEN-Marquardt_lab/Quentin/TIF-Seq/Data_processing/Quentin_secondBAM/EndTodEnd_align_polyXfiltered/clusters/raw/COL0cold.Aligned.sortedByCoord.out.dedup.uniq.merged.bed')
gr_spt16_r <- import('/home/bcm215/groupdirs/SCIENCE-PLEN-Marquardt_lab/Quentin/TIF-Seq/Data_processing/Quentin_secondBAM/EndTodEnd_align_polyXfiltered/clusters/raw/SPT16.Aligned.sortedByCoord.out.dedup.uniq.merged.bed')
gr_ssrp1_r <- import('/home/bcm215/groupdirs/SCIENCE-PLEN-Marquardt_lab/Quentin/TIF-Seq/Data_processing/Quentin_secondBAM/EndTodEnd_align_polyXfiltered/clusters/raw/SSRP1.Aligned.sortedByCoord.out.dedup.uniq.merged.bed')
