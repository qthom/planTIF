



df.genes <- as.data.frame(coding_genes)
df.genes.p <- df.genes[which(df.genes$strand == "+"), ]
df.genes.n <- df.genes[which(df.genes$strand == "-"), ]

gr <- gr_hen2
source("annotate.R")
hen2_an <- gr
df.hen2 <- as.data.frame(hen2_an)
df.hen2$geneID <- as.character(df.hen2$geneID)

gr <- gr_wt
source("annotate.R")
wt_an <- gr
df.wt <- as.data.frame(wt_an)
df.wt$geneID <- as.character(df.wt$geneID)


start = NULL
end = NULL
line <- df.wt


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



hen2 <- data.frame("distance_Start" = start,"distance_End" = end)


pdf('scatterhist_wt.pdf',10,8)

comparisonplot(wt$y,wt$x,pch=20, xlim=c(-1.0,0.2),ylim=c(-1.05,0.20), xlab = "Distance from PAS (normalized by gene length)", ylab = "Distance from TSS (normalized by gene length)",histbreak = 150, main = "TUs distance from annotated TSS and PAS" ) 
gradientLegend(c(0,1), color = c("#BEBEBE","#00008B","#FF0000","#FFA500","#FFD700"))
dev.off()

pdf('scatterhist_hen2.pdf',10,8)
comparisonplot(hen2$y,hen2$x,pch=20, xlim=c(-1.0,0.2),ylim=c(-1.05,0.20), xlab = "Distance from PAS (normalized by gene length)", ylab = "Distance from TSS (normalized by gene length)",histbreak = 150, main = "TUs distance from annotated TSS and PAS" ) 
gradientLegend(c(0,1), color = c("#BEBEBE","#00008B","#FF0000","#FFA500","#FFD700"))
dev.off()




gr <- gr_hen2cold
source("/home/bcm215/Desktop/Quentin/scripts/annotate.R")
hen2_ancold <- gr
df.hen2cold <- as.data.frame(hen2_ancold)
df.hen2cold$geneID <- as.character(df.hen2cold$geneID)

gr <- gr_wtcold
source("/home/bcm215/Desktop/Quentin/scripts/annotate.R")
wt_ancold <- gr
df.wtcold <- as.data.frame(wt_ancold)
df.wtcold$geneID <- as.character(df.wtcold$geneID)


start = NULL
end = NULL
line <- df.wtcold



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


hen2cold <- data.frame("distance_Start" = start,"distance_End" = end)



pdf('scatterhist_hen2_cold_4C.pdf',10,8)
comparisonplot(hen2cold$y,hen2cold$x,pch=20,xlim=c(-1.0,0.2),ylim=c(-1.05,0.20), xlab = "Distance from PAS (normalized by gene length)", ylab = "Distance from TSS (normalized by gene length)",histbreak = 150, main = "TUs distance from annotated TSS and PAS" ) 
gradientLegend(c(0,1), color = c("#BEBEBE","#00008B","#FF0000","#FFA500","#FFD700"))
dev.off()

pdf('scatterhist_WT_cold_4C.pdf',10,8)
comparisonplot(wtcold$y,wtcold$x,pch=20,xlim=c(-1.0,0.2),ylim=c(-1.05,0.20), xlab = "Distance from PAS (normalized by gene length)", ylab = "Distance from TSS (normalized by gene length)",histbreak = 150, main = "TUs distance from annotated TSS and PAS" ) 
gradientLegend(c(0,1), color = c("#BEBEBE","#00008B","#FF0000","#FFA500","#FFD700"))
dev.off()


gr <- gr_ssrp1
source("/home/bcm215/Desktop/Quentin/scripts/annotate.R")
ssrp1_an <- gr
df.ssrp1 <- as.data.frame(ssrp1_an)
df.ssrp1$geneID <- as.character(df.ssrp1$geneID)

gr <- gr_spt16
source("/home/bcm215/Desktop/Quentin/scripts/annotate.R")
wt_spt16 <- gr
df.spt16 <- as.data.frame(wt_spt16)
df.spt16$geneID <- as.character(df.spt16$geneID)


start = NULL
end = NULL
line <- df.ssrp1



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


ssrp1 <- data.frame("distance_Start" = start,"width" = end)

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



pdf('scatterhist_ssrp1.pdf',10,8)
comparisonplot(ssrp1$y,ssrp1$x,pch=20,xlim=c(-1.0,0.2),ylim=c(-1.05,0.20), xlab = "Distance from PAS (normalized by gene length)", ylab = "Distance from TSS (normalized by gene length)",histbreak = 150, main = "TUs distance from annotated TSS and PAS" ) 
gradientLegend(c(0,1), color = c("#BEBEBE","#00008B","#FF0000","#FFA500","#FFD700"))
dev.off()

pdf('scatterhist_spt16.pdf',10,8)
comparisonplot(spt16$y,spt16$x,pch=20,xlim=c(-1.0,0.2),ylim=c(-1.05,0.20), xlab = "Distance from PAS (normalized by gene length)", ylab = "Distance from TSS (normalized by gene length)",histbreak = 150, main = "TUs distance from annotated TSS and PAS" ) 
gradientLegend(c(0,1), color = c("#BEBEBE","#00008B","#FF0000","#FFA500","#FFD700"))
dev.off()

