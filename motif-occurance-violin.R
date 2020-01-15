

#This script was used to extract the relevant enrichment information from motif density plots

gaga <- as.data.frame(GAGA)
gaga <- t(gaga[-1])

write.table(gaga, file="gaga.csv", quote=F, sep="\t", row.names=T, col.names=T)

#from here the unrelevant columns were removed manually from the file

data <- read.csv("gaga.csv",sep=",",header=TRUE) 

pdf("gaga-mean-motif.pdf")
ggplot(data, aes(x=genes, y=mean, fill=genes)) + 
    geom_violin() + 
    geom_boxplot(width=0.05) +
    theme_classic() + 
    geom_signif(comparisons = list(c("ctrl","sppRNA")))

dev.off()

tcp1 <- as.data.frame(TCP1)
tcp1 <- t(tcp1[-1])

write.table(tcp1, file="tcp1.csv", quote=F, sep="\t", row.names=T, col.names=T)

#from here the unrelevant columns were removed manually from the file

data <- read.csv("tcp1.csv",sep=",",header=TRUE) 
pdf("tcp1-mean-motif.pdf")


ggplot(data, aes(x=genes, y=mean, fill=genes)) + 
    geom_violin() + 
    geom_boxplot(width=0.05) +
    theme_classic() + 
    geom_signif(comparisons = list(c("ctrl","sppRNA")))


dev.off()
