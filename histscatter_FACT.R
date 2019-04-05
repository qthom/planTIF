
fact_TSS <- "import the FACT TSS"
gr_factTSS <- GRanges(seqnames = fact_TSS$V1, ranges = IRanges(start = fact_TSS$V2, end = fact_TSS$V3), strand = fact_TSS$V4, geneID = fact_TSS$V9)

ssrp1_TSS <- resize(gr_ssrp1_PN, width = 1 , fix = "start")
gr_ssrp1_PN$overlapping_fact <- countOverlaps(ssrp1_TSS,gr_factTSS,type=c("any"),ignore.strand=FALSE)
gr_fact_ssrp1 <- gr_ssrp1_PN[gr_ssrp1_PN$overlapping_fact > 0,]


spt16_TSS <- resize(gr_spt16_PN, width = 1 , fix = "start")
gr_spt16_PN$overlapping_fact <- countOverlaps(spt16_TSS,gr_factTSS,type=c("any"),ignore.strand=FALSE)
gr_fact_spt16 <- gr_spt16_PN[gr_spt16_PN$overlapping_fact > 0,]

lii <- list(gr_fact_ssrp1,gr_fact_spt16)
fact_gr <- do.call("c", lii)



df.genes <- as.data.frame(coding_genes)
df.genes.p <- df.genes[which(df.genes$strand == "+"), ]
df.genes.n <- df.genes[which(df.genes$strand == "-"), ]

gr <- fact_gr
source("annotate.R")
spt16_an <- gr
df.spt16 <- as.data.frame(spt16_an)
df.spt16$geneID <- as.character(df.spt16$geneID)


start = NULL
end = NULL
line <- df.spt16

for (i in 1:nrow(line)) {
    sign <- line[i, "strand"]
    startc <- line[i, "start"]
    endc <- line[i,"end"]
    genes <- line[i,"geneID"]
    if (sign == "+") {
        gg <- df.genes[which(df.genes.p$names == genes ), ]$names
        startg <- df.genes.p[which(df.genes.p$names == genes ), ]$start
        endg <- df.genes.p[which(df.genes.p$names == genes ), ]$end
        width <- df.genes.p[which(df.genes.p$names == genes ), ]$width
        x <- ((startg - startc)/width)
        y <- ((endc-endg)/width)
        start = rbind(start, data.frame(x))
        end = rbind(end, data.frame(y))
    }

     if (sign == "-") {
        gg <- df.genes.n[which(df.genes.n$names == genes ), ]$names
        startg <- df.genes.n[which(df.genes.n$names == genes ), ]$start
        endg <- df.genes.n[which(df.genes.n$names == genes ), ]$end
        width <- df.genes.n[which(df.genes.n$names == genes ), ]$width
        x <- ((endc - endg)/width)
        y <- ((startg-startc)/width)
        start = rbind(start, data.frame(x))
        end = rbind(end, data.frame(y))
    }
}


spt16 <- data.frame("distance_Start" = start,"distance_End" = end)


pdf('scatterhist_fact.pdf',10,8)

comparisonplot(spt16$y,spt16$x,pch=20, xlab = "Distance from PAS (normalized by gene length)", ylab = "Distance from TSS (normalized by gene length)",histbreak = 150, main = "TUs distance from annotated TSS and PAS" ) 
gradientLegend(c(0,1), color = c("#BEBEBE","#00008B","#FF0000","#FFA500","#FFD700"))

dev.off()
