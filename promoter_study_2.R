set.seed(42) # for reproducibility when random
library(tidyverse) # data wrangling package
library(magrittr) # pipes
library(reshape2) # data wrangling package
library(ggplot2) # plotting package
library(RColorBrewer) # colors
library(BiocParallel) # parallelization 
register(MulticoreParam(workers=8)) # put your number of CPU cores here
library(TxDb.Athaliana.BioMart.plantsmart28) # Arabidopsis thaliana transcript database: TAIR10 from Ensemble, release 28
library(BSgenome.Athaliana.TAIR.TAIR9) # Arabidopsis thaliana full genome sequence (TAIR9 is same sequence as TAIR10: so ok)
library(stringr) # string manipulation
library(seqPattern) # to look up patterns
library(JASPAR2016) # JASPAR library
library(GenomicRanges) # to work with genome coordinates (GRanges)
options(scipen=10) # disable scientific notation
library(TFBSTools) # Transcription Factor Binding Sites Toolkit
library(motifStack) # to look up motifs

genome <- BSgenome.Athaliana.TAIR.TAIR9 # easier to work with this name
genome@seqinfo@genome[] <- "TAIR10" # little hack so that genome versions are compatible because oriignal BSgenome is from TAIR9 and might give an error otherwise


remove_out_of_bound <- function(GR) {idx=GenomicRanges:::get_out_of_bound_index(GR) # function to remove out-of-chr-boundaries regions from a GR object
                                     if(length(idx) != 0) { o <- GR[-idx]}
                                     else {o <- GR}
                                     o}


#make a gr object of canonical TSS positions based on TSS-Seq data
source("CAGEfightR_promoter_study.R")

gr <- GRanges(seqnames = seqnames(rse), ranges=mcols(rse)$thick, strand = strand(rse), geneID = mcols(rse)$geneID, txType = mcols(rse)$txType, score = mcols(rse)$score)
gr <- dropSeqlevels(gr,c("Mt","Pt"),pruning.mode ="coarse")
seqlevels(gr) <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")
gr <- gr[(gr$txType == "promoter")|(gr$txType == "proximal")|(gr$txType == "fiveUTR"),]
df <- as.data.frame(gr)
df <- df[order(df$geneID, -df$score),]
df$duplicate <- duplicated(df$geneID)
df <- df[df$duplicate == FALSE,]
gr <- GRanges(seqnames = df$seqnames, ranges = IRanges(start = df$start,end = df$end), strand = df$strand, geneID = df$geneID)
names(gr) <- gr$geneID

#make gr object for sppRNA and control genes canonical TSSs

df$matchsppRNA <- match(df$geneID,both2$geneID, nomatch =0)
df$matchnosppRNA <- match(df$geneID,ctrlboth2$geneID, nomatch =0)
sppRNA <- df[df$matchsppRNA > 0,]
ctrl <- df[df$matchnosppRNA > 0,]
gr_sppRNA <- GRanges(seqnames = sppRNA$seqnames, ranges = IRanges(start = sppRNA$start,end = sppRNA$end), strand = sppRNA$strand, geneID = sppRNA$geneID,sample = 'sppRNA')
names(gr_sppRNA) <- gr_sppRNA$geneID
gr_ctrl <- GRanges(seqnames = ctrl$seqnames, ranges = IRanges(start = ctrl$start,end = ctrl$end), strand = ctrl$strand, geneID = ctrl$geneID,sample = 'ctrl')
names(gr_ctrl) <- gr_ctrl$geneID


df_gr <- as.data.frame(gr)
#here we join the pausing index value to the TSS position calculated by TSS-Seq and CAGE-fightR
leftJoinDf <- left_join(df_gr,df.3,by="geneID")                                                                                                                                           
df_PIndex <- leftJoinDf[complete.cases(leftJoinDf), ]
gr_PI <- GRanges(seqnames = df_PIndex$seqnames, ranges = IRanges(start = df_PIndex$start,end = df_PIndex$end),strand = df_PIndex$strand, geneID = df_PIndex$geneID,PausingIndex = df_PIndex$PD, size = df_PIndex$size,match = df_PIndex$matchsppRNA)


upstream=200 # bp
downstream=150

seqPI <- gr_PI %>% # this is a GR object of 1-bp TSSs
  promoters(upstream=upstream, downstream=downstream) %>% # make promoters GR regions
  remove_out_of_bound() %>% # remove out of Chr (if any)
  getSeq(genome, .) # fetch genome sequence as a DNAStringSet


p<-plotPatternDensityMap(regionsSeq = seqPI,
             patterns = c("RGCCCAW"),
             outFile = "motifs_cordered_by_sppRNA",
             flankUp = upstream, flankDown = downstream,
             cexLabel = 8,
             seqOrder = order(-gr_PI$match),
             bandWidth = c(1,18), color = c("red"), addPatternLabel = FALSE)
        
p<-plotPatternDensityMap(regionsSeq = seqPI,
             patterns = c("TATAWA"),
             outFile = "motifs_cordered_by_sppRNA",
             flankUp = upstream, flankDown = downstream,
             cexLabel = 8,
             seqOrder = order(-gr_PI$match),
             bandWidth = c(1,18), color = c("red"), addPatternLabel = FALSE)
        
p<-plotPatternDensityMap(regionsSeq = seqPI,
             patterns = c("GAGAR"),
             outFile = "motifs_cordered_by_sppRNA",
             flankUp = upstream, flankDown = downstream,
             cexLabel = 8,
             seqOrder = order(-gr_PI$match),
             bandWidth = c(1,18), color = c("red"), addPatternLabel = FALSE)
        


upstream=200 # bp
downstream=150


list <- list(gr_ctrl,gr_sppRNA)
list <- do.call("c",list)
names(list)<-list$sample 


seqlist <- list %>% # this is a GR object of 1-bp TSSs
  promoters(upstream=upstream, downstream=downstream) %>% # make promoters GR regions
  remove_out_of_bound() %>% # remove out of Chr (if any)
  getSeq(genome, .) # fetch genome sequence as a DNAStringSet


TATA <-plotPatternDensityMap(regionsSeq = seqlist,
             patterns = c("TATAWA"),
             outFile = "motifs_twosets_cordered_by_genesets",
             flankUp = upstream, flankDown = downstream,
             cexLabel = 8,
             seqOrder = order(list$sample),
             bandWidth = c(1,14), color = c("red"), addPatternLabel = FALSE)
      
TCP1 <-plotPatternDensityMap(regionsSeq = seqlist,
             patterns = c("RGCCCAW"),
             outFile = "motifs_twosets_cordered_by_genesets",
             flankUp = upstream, flankDown = downstream,
             cexLabel = 8,
             seqOrder = order(list$sample),
             bandWidth = c(1,14), color = c("red"), addPatternLabel = FALSE)
      

GAGA <-plotPatternDensityMap(regionsSeq = seqlist,
             patterns = c("GAGAR"),
             outFile = "motifs_twosets_cordered_by_genesets",
             flankUp = upstream, flankDown = downstream,
             cexLabel = 8,
             seqOrder = order(list$sample),
             bandWidth = c(1,14), color = c("red"), addPatternLabel = FALSE)
      


seqs <- gr_sppRNA %>% # this is a GR object of 1-bp TSSs
  promoters(upstream=upstream, downstream=downstream) %>% # make promoters GR regions
  remove_out_of_bound() %>% # remove out of Chr (if any)
  getSeq(genome, .) # fetch genome sequence as a DNAStringSet
seqc <- gr_ctrl %>% # this is a GR object of 1-bp TSSs
  promoters(upstream=upstream, downstream=downstream) %>% 
  remove_out_of_bound() %>% 
  getSeq(genome, .) 



pdf("occurrence in sppRNA genes.pdf")
plotPatternOccurrenceAverage(regionsSeq = seqs,
             patterns = c("TATAWA","RGCCCAW", "GAGAR"), flankUp = upstream, flankDown = downstream,xLabel = "Distance from the canonical TSS (bp)",
             smoothingWindow = 8, color = c("red3", "blue3","green3","black","purple"), cex.axis = 0.9)
dev.off()

pdf("motif occurrence in ctrl genes.pdf")
plotPatternOccurrenceAverage(regionsSeq = seqc,
             patterns = c("TATAWA","RGCCCAW", "GAGAR"), flankUp = upstream, flankDown = downstream,xLabel = "Distance from the canonical TSS (bp)",
             smoothingWindow = 8, color = c("red3", "blue3","green3","black","purple"), cex.axis = 0.9)
dev.off()
