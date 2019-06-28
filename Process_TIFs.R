#load the clusters

gr_wt <- import('COL0.bed')
gr_hen2 <- import('HEN2-2.bed')
gr_hen2cold <- import('/HEN2-2cold.bed')
gr_wtcold <- import('COL0cold.bed')
gr_spt16 <- import('SPT16.bed')
gr_ssrp1 <- import('SSRP1.bed')

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
