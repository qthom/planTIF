
#plot distribution of clusters along the start and end in a histogram/scatter

library("LSD")
library("itsadug")

setwd("/home/bcm215/Desktop/Quentin/scripts/plots/")

df.genes <- as.data.frame(coding_genes)
df.genes.p <- df.genes[which(df.genes$strand == "+"), ]
df.genes.n <- df.genes[which(df.genes$strand == "-"), ]

gr <- gr_hen2
source("/home/bcm215/Desktop/Quentin/scripts/TIF-Seq_2019/annotate.R")
hen2_an <- gr
df.hen2 <- as.data.frame(hen2_an)
df.hen2$geneID <- as.character(df.hen2$geneID)
df.hen2 <- df.hen2[which(df.hen2$width < 200), ]

gr <- gr_wt
source("/home/bcm215/Desktop/Quentin/scripts/TIF-Seq_2019/annotate.R")
wt_an <- gr
df.wt <- as.data.frame(wt_an)
df.wt$geneID <- as.character(df.wt$geneID)
df.wt <- df.wt[which(df.wt$width < 200), ]


start = NULL
end = NULL
line <- df.wt


for (i in 1:nrow(line)) {
    sign <- line[i, "strand"]
    startc <- line[i, "start"]
    endc <- line[i,"end"]
    genes <- line[i,"geneID"]
    if (sign == "+") {
        gg <- df.genes[which(df.genes.p$gene_id == genes ), ]$gene_id
        startg <- df.genes.p[which(df.genes.p$gene_id == genes ), ]$start
        endg <- df.genes.p[which(df.genes.p$gene_id == genes ), ]$end
        width <- df.genes.p[which(df.genes.p$gene_id == genes ), ]$width
        x <- ((startg - startc)/width)
        y <- ((endc-endg)/width)
        start = rbind(start, data.frame(x))
        end = rbind(end, data.frame(y))
    }

     if (sign == "-") {
        gg <- df.genes.n[which(df.genes.n$gene_id == genes ), ]$gene_id
        startg <- df.genes.n[which(df.genes.n$gene_id == genes ), ]$start
        endg <- df.genes.n[which(df.genes.n$gene_id == genes ), ]$end
        width <- df.genes.n[which(df.genes.n$gene_id == genes ), ]$width
        x <- ((endc - endg)/width)
        y <- ((startg-startc)/width)
        start = rbind(start, data.frame(x))
        end = rbind(end, data.frame(y))
    }
}


wt <- data.frame("distance_Start" = start,"width" = end)

start = NULL
end = NULL
line <- df.hen2

for (i in 1:nrow(line)) {
    sign <- line[i, "strand"]
    startc <- line[i, "start"]
    endc <- line[i,"end"]
    genes <- line[i,"geneID"]
    if (sign == "+") {
        gg <- df.genes[which(df.genes.p$gene_id == genes ), ]$gene_id
        startg <- df.genes.p[which(df.genes.p$gene_id == genes ), ]$start
        endg <- df.genes.p[which(df.genes.p$gene_id == genes ), ]$end
        width <- df.genes.p[which(df.genes.p$gene_id == genes ), ]$width
        x <- ((startg - startc)/width)
        y <- ((endc-endg)/width)
        start = rbind(start, data.frame(x))
        end = rbind(end, data.frame(y))
    }

     if (sign == "-") {
        gg <- df.genes.n[which(df.genes.n$gene_id == genes ), ]$gene_id
        startg <- df.genes.n[which(df.genes.n$gene_id == genes ), ]$start
        endg <- df.genes.n[which(df.genes.n$gene_id == genes ), ]$end
        width <- df.genes.n[which(df.genes.n$gene_id == genes ), ]$width
        x <- ((endc - endg)/width)
        y <- ((startg-startc)/width)
        start = rbind(start, data.frame(x))
        end = rbind(end, data.frame(y))
    }
}


hen2 <- data.frame("distance_Start" = start,"distance_End" = end)


pdf('scatterhist_wt_smalltranscripts.pdf',10,8)

comparisonplot(wt$y,wt$x,pch=20, xlab = "Distance from PAS (normalized by gene length)", ylab = "Distance from TSS (normalized by gene length)",histbreak = 150, main = "TUs distance from annotated TSS and PAS" ) 
gradientLegend(c(0,1), color = c("#BEBEBE","#00008B","#FF0000","#FFA500","#FFD700"))
dev.off()

pdf('scatterhist_hen2_smalltranscripts.pdf',10,8)
comparisonplot(hen2$y,hen2$x,pch=20, xlab = "Distance from PAS (normalized by gene length)", ylab = "Distance from TSS (normalized by gene length)",histbreak = 150, main = "TUs distance from annotated TSS and PAS" ) 
gradientLegend(c(0,1), color = c("#BEBEBE","#00008B","#FF0000","#FFA500","#FFD700"))
dev.off()




gr <- gr_hen2cold
source("/home/bcm215/Desktop/Quentin/scripts/TIF-Seq_2019/annotate.R")
hen2_ancold <- gr
df.hen2cold <- as.data.frame(hen2_ancold)
df.hen2cold$geneID <- as.character(df.hen2cold$geneID)
df.hen2cold <- df.hen2cold[which(df.hen2cold$width < 180), ]
gr <- gr_wtcold
source("/home/bcm215/Desktop/Quentin/scripts/TIF-Seq_2019/annotate.R")
wt_ancold <- gr
df.wtcold <- as.data.frame(wt_ancold)
df.wtcold$geneID <- as.character(df.wtcold$geneID)
df.wtcold <- df.wtcold[which(df.wtcold$width < 180), ]


start = NULL
end = NULL
line <- df.wtcold


for (i in 1:nrow(line)) {
    sign <- line[i, "strand"]
    startc <- line[i, "start"]
    endc <- line[i,"end"]
    genes <- line[i,"geneID"]
    if (sign == "+") {
        gg <- df.genes[which(df.genes.p$gene_id == genes ), ]$gene_id
        startg <- df.genes.p[which(df.genes.p$gene_id == genes ), ]$start
        endg <- df.genes.p[which(df.genes.p$gene_id == genes ), ]$end
        width <- df.genes.p[which(df.genes.p$gene_id == genes ), ]$width
        x <- ((startg - startc)/width)
        y <- ((endc-endg)/width)
        start = rbind(start, data.frame(x))
        end = rbind(end, data.frame(y))
    }

     if (sign == "-") {
        gg <- df.genes.n[which(df.genes.n$gene_id == genes ), ]$gene_id
        startg <- df.genes.n[which(df.genes.n$gene_id == genes ), ]$start
        endg <- df.genes.n[which(df.genes.n$gene_id == genes ), ]$end
        width <- df.genes.n[which(df.genes.n$gene_id == genes ), ]$width
        x <- ((endc - endg)/width)
        y <- ((startg-startc)/width)
        start = rbind(start, data.frame(x))
        end = rbind(end, data.frame(y))
    }
}

wtcold <- data.frame("distance_Start" = start,"width" = end)

start = NULL
end = NULL
line <- df.hen2cold

for (i in 1:nrow(line)) {
    sign <- line[i, "strand"]
    startc <- line[i, "start"]
    endc <- line[i,"end"]
    genes <- line[i,"geneID"]
    if (sign == "+") {
        gg <- df.genes[which(df.genes.p$gene_id == genes ), ]$gene_id
        startg <- df.genes.p[which(df.genes.p$gene_id == genes ), ]$start
        endg <- df.genes.p[which(df.genes.p$gene_id == genes ), ]$end
        width <- df.genes.p[which(df.genes.p$gene_id == genes ), ]$width
        x <- ((startg - startc)/width)
        y <- ((endc-endg)/width)
        start = rbind(start, data.frame(x))
        end = rbind(end, data.frame(y))
    }

     if (sign == "-") {
        gg <- df.genes.n[which(df.genes.n$gene_id == genes ), ]$gene_id
        startg <- df.genes.n[which(df.genes.n$gene_id == genes ), ]$start
        endg <- df.genes.n[which(df.genes.n$gene_id == genes ), ]$end
        width <- df.genes.n[which(df.genes.n$gene_id == genes ), ]$width
        x <- ((endc - endg)/width)
        y <- ((startg-startc)/width)
        start = rbind(start, data.frame(x))
        end = rbind(end, data.frame(y))
    }
}


hen2cold <- data.frame("distance_Start" = start,"distance_End" = end)


pdf('scatterhist_wt_smalltranscripts_4c.pdf',10,8)

comparisonplot(wtcold$y,wtcold$x,pch=20, xlab = "Distance from PAS (normalized by gene length)", ylab = "Distance from TSS (normalized by gene length)",histbreak = 150, main = "TUs distance from annotated TSS and PAS" ) 
gradientLegend(c(0,1), color = c("#BEBEBE","#00008B","#FF0000","#FFA500","#FFD700"))

dev.off()

pdf('scatterhist_hen2_smalltranscripts_4c.pdf',10,8)
comparisonplot(hen2cold$y,hen2cold$x,pch=20, xlab = "Distance from PAS (normalized by gene length)", ylab = "Distance from TSS (normalized by gene length)",histbreak = 150, main = "TUs distance from annotated TSS and PAS" ) 
gradientLegend(c(0,1), color = c("#BEBEBE","#00008B","#FF0000","#FFA500","#FFD700"))
dev.off()
