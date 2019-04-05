
wt <- gr_wt[width(gr_wt) < 6000]
hen2<- gr_hen2[width(gr_hen2) < 6000]
wtcold <- gr_wtcold[width(gr_wtcold) < 6000]
hen2cold <- gr_hen2cold[width(gr_hen2cold) < 6000]

dat1 <- data.frame(lines = factor("Col-0" ), width = width(wt))
dat2 <- data.frame(lines = factor("hen2-2"), width = width(hen2))
dat3 <- data.frame(lines = factor("COl.0 Cold"), width = width(wtcold))
dat4 <- data.frame(lines = factor("hen2-2 Cold"), width = width(hen2cold))


dat <- rbind(dat2,dat1)
datcold <- rbind(dat4,dat3)


pdf("violin_size_TUs_22C.pdf")
ggplot(dat, aes(x=lines, y=width, fill=lines)) + geom_violin() + geom_boxplot(width=0.05) + coord_flip() +theme_classic() + scale_fill_manual(values=c("#4f83e2", "#e2574d"))
dev.off()


pdf("violin_size_TUs_4C.pdf")
ggplot(datcold, aes(x=lines, y=width, fill=lines)) + geom_violin() + geom_boxplot(width=0.05) + coord_flip() +theme_classic() + scale_fill_manual(values=c("#4f83e2", "#e2574d"))
dev.off()


