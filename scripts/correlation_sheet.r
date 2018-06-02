library(cluster)
require(lattice)
library(corrplot)
args <- commandArgs(TRUE)
dat_path <- args[1]
pdf_path <- args[2]

corr <- read.table(dat_path)
rownames(corr) <- 0:(nrow(corr)-1)
colnames(corr) <- rownames(corr)
#colnames(corr) <- 0 : rows(corr) - 1;
bwr <- c("#4080FF", "#4080FF", "red")
#bwr <- c("white", "#4080FF", "red")
col3 <- colorRampPalette(bwr, space = "rgb")

pdf(pdf_path, height=25, width=25)
lattice.options(axis.padding=list(factor=0.5))

x.scale <- list(cex=2, alternating=2, col='black',rot=90)
y.scale <- list(cex=2, alternating=1, col='black')

q  = t(as.matrix(corr))
q2 = q[,] 
#corrplot(q2,col=col3(100), tl.col="black", tl.srt=45)
corrplot(q2,col=col3(100), tl.col="black", tl.srt=45, order = "FPC")

#image(q2, col=col3(100))
dev.off()
