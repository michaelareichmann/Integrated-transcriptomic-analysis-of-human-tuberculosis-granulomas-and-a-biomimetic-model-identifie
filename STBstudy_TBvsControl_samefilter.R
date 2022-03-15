library(limma)
library(ggplot2)
library(dplyr)
library(edgeR)
library(RUVSeq)
library(RColorBrewer)
library(rgl)
library(EDASeq)

###Load data and metadata

data <- read.csv('/Users/mr4e14/Documents/STBstudy/STBstudy_tximport_kallisto/STBstudy_tximport_results.unfiltered.csv',row.names = 1)
metadata <- read.csv('/Users/mr4e14/Documents/Rdocuments/R_STBstudy_sheets/STB_sample_description_complete.csv')
head(metadata)

samples <- metadata$Sample
samples <- as.character(samples)
samples
data <- data[,samples]
colnames(data)
head(data)

chip <- factor(metadata$Chip)
chip 
day <- factor(metadata$Day)
day
sex <- factor(metadata$Sex)
sex
location <- factor(metadata$Location)
location
race <- factor(metadata$Race)
race
Group <- factor(metadata$Condition)
Group

design <- model.matrix(~Group)
design
colnames(design)

ncol(design)

FILTER=8

DGEList <- DGEList(data,samples = metadata$Sample, genes=row.names(data))

table(rowSums(DGEList$counts==0)==24)
cpm <- cpm(DGEList)
colnames(cpm) <- metadata$Sample
keepers <- rowSums(cpm>FILTER)>4

DGEList.filtered <- DGEList[keepers,]
colnames(DGEList.filtered) <- metadata$Sample
counts_filtered<-DGEList.filtered$counts

dim(counts_filtered)

eds<-DGEList(counts=counts_filtered, genes=rownames(counts_filtered))
eds<-calcNormFactors(eds, method="TMM")

DGEList.filtered.norm <- calcNormFactors(eds, method = "TMM")

v <- voomWithQualityWeights(DGEList.filtered.norm, design=design, normalization="none", plot=TRUE)

fit <- lmFit(v,design)
efit <- eBayes(fit)
ncol(fit)
topTable(efit, coef=3, adjust="BH")
topTable <- topTable(efit, coef=3, adjust="BH", number=10000, sort.by="P")
topTable
write.csv(topTable, file=paste("STBstudy_3conditions_TBvsControl_limma.All_results.","filter_",FILTER,"_TMM",".csv",sep=""))

summary(decideTests(efit))

#Plot PCA
Y <- apply(cpm, 1, function(y) scale(y, center=TRUE, scale=FALSE))
s <- svd(Y)
col.cell <- c("slateblue1","turquoise", "salmon1")[metadata$Condition]
col.cell

pdf("PDF_3conditions_limma_filter_8_TMM_with_legend.pdf")
plot(s$u[,1], s$u[,2], pch=19, cex=4, xlab="PC1", ylab="PC2", col=col.cell,
     ylim=c((min(s$u[,2])*1.2),(max(s$u[,2])*1.2)),xlim=c((min(s$u[,1])*1.2),(max(s$u[,1])*1.2)),
     cex.lab=1.5)
#text(s$u[,1], s$u[,2], labels=colnames(cpm),adj=c(0,2),cex=0.5)

par(xpd = TRUE)

legend(0.12, 0.65, xpd = TRUE,
       legend = c("Control", "Sarcoid", "TB"), 
       col = unique(metadata$Condition),
       fill = c("slateblue1", "turquoise", "salmon1"),
       border = FALSE,
       bty = "n",
       cex = 1.8)
dev.off()

png("PCA_3conditions_limma_filter_8_TMM_with_legend.png",
    width = 10*300,
    height = 10*300,
    res = 300)
plot(s$u[,1], s$u[,2], pch=19, cex=4, xlab="PC1", ylab="PC2",main="PCA plot TMM Normalization",col=col.cell,
     ylim=c((min(s$u[,2])*1.2),(max(s$u[,2])*1.2)),xlim=c((min(s$u[,1])*1.2),(max(s$u[,1])*1.2)))
text(s$u[,1], s$u[,2], labels=colnames(cpm),adj=c(0,2),cex=0.5)

par(xpd = TRUE)

legend(0.2, 0.64, xpd = TRUE,
       legend = c("TB", "Sarcoid", "Control"), 
       col = unique(metadata$Condition),
       fill = c("salmon1", "turquoise", "slateblue1"),
       border = FALSE,
       bty = "n",
       cex = 1.6)
dev.off()

##
col.cell <- c("blue3", "red3", "springgreen3", "yellow", "magenta2", "dodgerblue2")[metadata$Chip]
col.cell

png("PCA_chip_with_legend.png",
    width = 10*300,
    height = 10*300,
    res = 300)

pdf("PCA_chip_with_legend.pdf")

plot(s$u[,1], s$u[,2], pch=19, cex=4, xlab="PC1", ylab="PC2", col=col.cell,
     ylim=c((min(s$u[,2])*1.2),(max(s$u[,2])*1.2)),xlim=c((min(s$u[,1])*1.2),(max(s$u[,1])*1.2)),
     cex.lab=1.3)
#text(s$u[,1], s$u[,2], labels=colnames(cpm),adj=c(0,2),cex=0.5)

par(xpd = TRUE)

legend(0.19, 0.66, xpd = TRUE,
       legend = c("Chip 1", "Chip 2", "Chip 3", "Chip 4", "Chip 5", "Chip 6"), 
       col = unique(metadata$Chip),
       fill = c("blue3", "red3", "springgreen3", "yellow", "magenta2", "dodgerblue2"),
       border = FALSE,
       bty = "n",
       cex = 1.3)
dev.off()

##
col.cell <- c("springgreen3", "yellow", "magenta2")[metadata$Day]
col.cell

png("PCA_day_with_legend.png",
    width = 10*300,
    height = 10*300,
    res = 300)

pdf("PCA_day_with_legend.pdf")
plot(s$u[,1], s$u[,2], pch=19, cex=4, xlab="PC1", ylab="PC2", col=col.cell,
     ylim=c((min(s$u[,2])*1.2),(max(s$u[,2])*1.2)),xlim=c((min(s$u[,1])*1.2),(max(s$u[,1])*1.2)),
     cex.lab=1.3)
#text(s$u[,1], s$u[,2], labels=colnames(cpm),adj=c(0,2),cex=0.5)

par(xpd = TRUE)

legend(0.19, 0.66, xpd = TRUE,
       legend = c("Day X", "Day Y", "Day Z"), 
       col = unique(metadata$Day),
       fill = c("springgreen3", "yellow", "magenta2"),
       border = FALSE,
       bty = "n",
       cex = 1.3)
dev.off()

##
col.cell <- c("red3","blue3")[metadata$Sex]
col.cell

png("PCA_gender_with_legend.png",
    width = 10*300,
    height = 10*300,
    res = 300)

pdf("PCA_gender_with_legend.pdf")
plot(s$u[,1], s$u[,2], pch=19, cex=4, xlab="PC1", ylab="PC2", col=col.cell,
     ylim=c((min(s$u[,2])*1.2),(max(s$u[,2])*1.2)),xlim=c((min(s$u[,1])*1.2),(max(s$u[,1])*1.2)),
     cex.lab=1.3)
#text(s$u[,1], s$u[,2], labels=colnames(cpm),adj=c(0,2),cex=0.5)

par(xpd = TRUE)

legend(0.18, 0.66, xpd = TRUE,
       legend = c("Male", "Female"), 
       col = unique(metadata$Sex),
       fill = c("blue3", "red3"),
       border = FALSE,
       bty = "n",
       cex = 1.3)
dev.off()

##
col.cell <- c("blue3","red3")[metadata$Location]
col.cell

png("PCA_location_with_legend.png",
    width = 10*300,
    height = 10*300,
    res = 300)

pdf("PCA_location_with_legend.pdf")
plot(s$u[,1], s$u[,2], pch=19, cex=4, xlab="PC1", ylab="PC2", col=col.cell,
     ylim=c((min(s$u[,2])*1.2),(max(s$u[,2])*1.2)),xlim=c((min(s$u[,1])*1.2),(max(s$u[,1])*1.2)),
     cex.lab=1.3)
#text(s$u[,1], s$u[,2], labels=colnames(cpm),adj=c(0,2),cex=0.5)

par(xpd = TRUE)

legend(0.11, 0.66, xpd = TRUE,
       legend = c("Mediastinum", "Neck"), 
       col = unique(metadata$Location),
       fill = c("blue3","red3"),
       border = FALSE,
       bty = "n",
       cex = 1.3)
dev.off()

##
col.cell <- c("blue3", "red3", "springgreen3", "yellow", "magenta2")[metadata$Race]
col.cell

png("PCA_race_with_legend.png",
    width = 10*300,
    height = 10*300,
    res = 300)

pdf("PCA_race_with_legend.pdf")
plot(s$u[,1], s$u[,2], pch=19, cex=4, xlab="PC1", ylab="PC2", col=col.cell,
     ylim=c((min(s$u[,2])*1.2),(max(s$u[,2])*1.2)),xlim=c((min(s$u[,1])*1.2),(max(s$u[,1])*1.2)),
     cex.lab=1.3)
#text(s$u[,1], s$u[,2], labels=colnames(cpm),adj=c(0,2),cex=0.5)

par(xpd = TRUE)

legend(0.15, 0.66, xpd = TRUE,
       legend = c("Filipino", "Georgian", "Indian", "S.African", "Unknown"), 
       col = unique(metadata$Location),
       fill = c("blue3", "red3", "springgreen3", "yellow", "magenta2"),
       border = FALSE,
       bty = "n",
       cex = 1.3)
dev.off()

#Plot IQR vs Median TMM normalised
col.cell <- c("slateblue1", "turquoise", "salmon1")

log_counts<-log(counts_filtered+1)

IQR<-apply(cpm, 2, IQR)
Median<-apply(cpm, 2, median)
diff1<-mean(Median)-min(Median)
diff2<-max(Median)-mean(Median)
diff3<-mean(IQR)-min(IQR)
diff4<-max(IQR)-mean(IQR)
Xlim=c(mean(Median)-2*diff1,mean(Median)+2*diff2)
Ylim=c(mean(IQR)-2*diff3,mean(IQR)+2*diff4)

#plot(Median, IQR, main="IQR vs. Median TMM normalisation", type="n", xlim=Xlim,ylim=Ylim, col=col.cell)
#text(Median, IQR, labels=names(IQR),col=col.cell, cex=0.5)
plot(Median, IQR, main="IQR vs. Median TMM normalisation", type="n", xlim=c(-0.2, 0.3),ylim=c(3.5, 20), col=col.cell, )
text(Median, IQR, labels=names(IQR),col=col.cell, cex=1.6)

#  Make boxes for StDev.
Median_mean<-mean(Median)
c_sd1_mean<-sd(Median)
c_sd2_mean<-2*sd(Median)
c_sd3_mean<-3*sd(Median)
IQR_mean<-mean(IQR)
c_sd1_IQR<-sd(IQR)
c_sd2_IQR<-2*sd(IQR)
c_sd3_IQR<-3*sd(IQR)

x0_c<-Median_mean-c_sd1_mean
y0_c<-IQR_mean-c_sd1_IQR
x1_c<-Median_mean+c_sd1_mean
y1_c<-IQR_mean+c_sd1_IQR

x0_c.2<-Median_mean-c_sd2_mean
y0_c.2<-IQR_mean-c_sd2_IQR
x1_c.2<-Median_mean+c_sd2_mean
y1_c.2<-IQR_mean+c_sd2_IQR

x0_c.3<-Median_mean-c_sd3_mean
y0_c.3<-IQR_mean-c_sd3_IQR
x1_c.3<-Median_mean+c_sd3_mean
y1_c.3<-IQR_mean+c_sd3_IQR

segments(x0_c,y0_c, x1=x1_c, y1=y0_c, col="blue")
segments(x0_c,y0_c, x1=x0_c, y1=y1_c, col="blue") 
segments(x1_c,y0_c, x1=x1_c, y1=y1_c, col="blue") 
segments(x0_c,y1_c, x1=x1_c, y1=y1_c, col="blue")

segments(x0_c.2,y0_c.2, x1=x1_c.2, y1=y0_c.2, col="red") 
segments(x0_c.2,y0_c.2, x1=x0_c.2, y1=y1_c.2, col="red") 
segments(x1_c.2,y0_c.2, x1=x1_c.2, y1=y1_c.2, col="red")
segments(x0_c.2,y1_c.2, x1=x1_c.2, y1=y1_c.2, col="red")

#  This portion is out of range so I removed it.
segments(x0_c.3,y0_c.3, x1=x1_c.3, y1=y0_c.3, col="green")
segments(x0_c.3,y0_c.3, x1=x0_c.3, y1=y1_c.3, col="green")
segments(x1_c.3,y0_c.3, x1=x1_c.3, y1=y1_c.3, col="green")
segments(x0_c.3,y1_c.3, x1=x1_c.3, y1=y1_c.3, col="green")











