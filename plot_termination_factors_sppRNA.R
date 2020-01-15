#this script plots the distribution of PAS peaks reads of sppRNA genes in relation with a control set
#all control set of genes were calculated with the folowing script with impot mut (mutant) and WT
#source("/science/groupdirs/bcm215/SCIENCE-PLEN-Marquardt_lab/Quentin/scripts/TIF-Seq_2019/new/equa_termination_factors.R")

#import the controls sets for each data set
list <- readRDS("/home/bcm215/Desktop/plots termination analysis/list_controlsppRNA.rda")

#CTSF77

dat1 <- data.frame("geno"="cstf77 ctrl", "count" = list$cstf77ctrl$mut)
dat2 <- data.frame("geno"="WT ctrl", "count" = list$cstf77ctrl$WT)
dat3 <- data.frame("geno"="cstf77 sppRNA","count" = list$cstf77sppRNA$mut)
dat4 <- data.frame("geno"="WT sppRNA","count" = list$cstf77sppRNA$WT)
dat <- rbind(dat1,dat3,dat2,dat4)

pdf("cstf77_sppRNA.pdf")
ggplot(dat, aes(x=geno, y=count, fill=geno)) + 
                        geom_boxplot(width=0.5) + 
                        geom_signif(comparisons = list(c("cstf77 ctrl", "cstf77 sppRNA")), map_signif_level=TRUE) +
                        geom_signif(comparisons = list(c("WT ctrl", "WT sppRNA")), map_signif_level=TRUE) +
                        theme_classic() +  
                        labs(y = "TPM at canonical PAS")

dev.off()

#treated CSTF77

dat1 <- data.frame("geno"="cstf77 ctrl", "count" = list$Tcstf77ctrl$mut)
dat2 <- data.frame("geno"="WT ctrl", "count" = list$Tcstf77ctrl$WT)
dat3 <- data.frame("geno"="cstf77 sppRNA","count" = list$Tcstf77sppRNA$mut)
dat4 <- data.frame("geno"="WT sppRNA","count" = list$Tcstf77sppRNA$WT)
dat <- rbind(dat1,dat3,dat2,dat4)

pdf("Tcstf77_sppRNA.pdf")
ggplot(dat, aes(x=geno, y=count, fill=geno)) + 
                        geom_boxplot(width=0.5) + 
                        geom_signif(comparisons = list(c("cstf77 ctrl", "cstf77 sppRNA")), map_signif_level=TRUE) +
                        geom_signif(comparisons = list(c("WT ctrl", "WT sppRNA")), map_signif_level=TRUE) +
                        theme_classic() +  
                        labs(y = "TPM at canonical PAS")

dev.off()


#CPSF100

dat1 <- data.frame("geno"="CPSF100 ctrl", "count" = list$cpsf100ctrl$mut)
dat2 <- data.frame("geno"="WT ctrl", "count" = list$cpsf100ctrl$WT)
dat3 <- data.frame("geno"="CPSF100 sppRNA","count" = list$cpsf100sppRNA$mut)
dat4 <- data.frame("geno"="WT sppRNA","count" = list$cpsf100sppRNA$WT)
dat <- rbind(dat1,dat3,dat2,dat4)

pdf("cpsf100_sppRNA.pdf")
ggplot(dat, aes(x=geno, y=count, fill=geno)) + 
                        geom_boxplot(width=0.5) + 
                        geom_signif(comparisons = list(c("CPSF100 ctrl", "CPSF100 sppRNA")), map_signif_level=TRUE) +
                        geom_signif(comparisons = list(c("WT ctrl", "WT sppRNA")), map_signif_level=TRUE) +
                        theme_classic() +  
                        labs(y = "TPM at canonical PAS")
dev.off()



#FY2


dat1 <- data.frame("geno"="fy2 ctrl", "count" = list$fy2ctrl$mut)
dat2 <- data.frame("geno"="WT ctrl", "count" = list$fy2ctrl$WT)
dat3 <- data.frame("geno"="fy2 sppRNA","count" = list$fy2sppRNA$mut)
dat4 <- data.frame("geno"="WT sppRNA","count" = list$fy2sppRNA$WT)
dat <- rbind(dat1,dat3,dat2,dat4)

pdf("fy2_sppRNA.pdf")
ggplot(dat, aes(x=geno, y=count, fill=geno)) + 
                        geom_boxplot(width=0.5) + 
                        geom_signif(comparisons = list(c("fy2 ctrl", "fy2 sppRNA")), map_signif_level=TRUE) +
                        geom_signif(comparisons = list(c("WT ctrl", "WT sppRNA")), map_signif_level=TRUE) +
                        theme_classic() +  
                        labs(y = "TPM at canonical PAS")
dev.off()

#FY6

dat1 <- data.frame("geno"="fy6 ctrl", "count" = list$fy6ctrl$mut)
dat2 <- data.frame("geno"="WT ctrl", "count" = list$fy6ctrl$WT)
dat3 <- data.frame("geno"="fy6 sppRNA","count" = list$fy6sppRNA$mut)
dat4 <- data.frame("geno"="WT sppRNA","count" = list$fy6sppRNA$WT)
dat <- rbind(dat1,dat3,dat2,dat4)

pdf("fy6_sppRNA.pdf")
ggplot(dat, aes(x=geno, y=count, fill=geno)) + 
                        geom_boxplot(width=0.5) + 
                        geom_signif(comparisons = list(c("fy6 ctrl", "fy6 sppRNA")), map_signif_level=TRUE) +
                        geom_signif(comparisons = list(c("WT ctrl", "WT sppRNA")), map_signif_level=TRUE) +
                        theme_classic() +  
                        labs(y = "TPM at canonical PAS")
dev.off()


#FY3

dat1 <- data.frame("geno"="fy3 ctrl", "count" = list$fy3ctrl$mut)
dat2 <- data.frame("geno"="WT ctrl", "count" = list$fy3ctrl$WT)
dat3 <- data.frame("geno"="fy3 sppRNA","count" = list$fy3sppRNA$mut)
dat4 <- data.frame("geno"="WT sppRNA","count" = list$fy3sppRNA$WT)
dat <- rbind(dat1,dat3,dat2,dat4)

pdf("fy3_sppRNA.pdf")
ggplot(dat, aes(x=geno, y=count, fill=geno)) + 
                        geom_boxplot(width=0.5) + 
                        geom_signif(comparisons = list(c("fy3 ctrl", "fy3 sppRNA")), map_signif_level=TRUE) +
                        geom_signif(comparisons = list(c("WT ctrl", "WT sppRNA")), map_signif_level=TRUE) +
                        theme_classic() +  
                        labs(y = "TPM at canonical PAS")
dev.off()




#CPSF30

dat1 <- data.frame("geno"="CPSF30 ctrl", "count" = list$cpsf30ctrl$mut)
dat2 <- data.frame("geno"="WT ctrl", "count" = list$cpsf30ctrl$WT)
dat3 <- data.frame("geno"="CPSF30 sppRNA","count" = list$cpsf30sppRNA$mut)
dat4 <- data.frame("geno"="WT sppRNA","count" = list$cpsf30sppRNA$WT)
dat <- rbind(dat1,dat3,dat2,dat4)

pdf("cpsf30_sppRNA.pdf")
ggplot(dat, aes(x=geno, y=count, fill=geno)) + 
                        geom_boxplot(width=0.5) + 
                        geom_signif(comparisons = list(c("CPSF30 ctrl", "CPSF30 sppRNA")), map_signif_level=TRUE) +
                        geom_signif(comparisons = list(c("WT ctrl", "WT sppRNA")), map_signif_level=TRUE) +
                        theme_classic() +  
                        labs(y = "TPM at canonical PAS")
dev.off()
