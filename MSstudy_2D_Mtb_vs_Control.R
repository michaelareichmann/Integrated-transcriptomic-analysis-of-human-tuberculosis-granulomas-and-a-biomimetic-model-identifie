library(tximport)
library("biomaRt")
library(GSEABase)
library(GSVA)
library(Biobase)
library(edgeR)
library(limma)
library(ggplot2)
library(dplyr)
library(tidyverse) 
library(Biostrings)
library(ensembldb)
library(EnsDb.Hsapiens.v86) 
library(sleuth) 
library(RColorBrewer) 
library(reshape2) 
library(trelliscopejs) 
library(genefilter)
library(matrixStats) 
library(RUVSeq)
library(rgl)
library(EDASeq)

# Anotations table --------------------------------------------------------

listTables(EnsDb.Hsapiens.v86)
listColumns(EnsDb.Hsapiens.v86, "tx")
Tx <- transcripts(EnsDb.Hsapiens.v86, columns=c(listColumns(EnsDb.Hsapiens.v86,"tx"), "gene_name"))
Tx <- as_tibble(Tx)
Tx <- dplyr::rename(Tx, target_id = tx_id)
Tx <- dplyr::select(Tx, target_id, gene_name)
head(Tx)


# Load the data -----------------------------------------------------------

DATA="REPEAT_MSstudy_tximport_AlgCollagen" # name of the folder for the results 
c1="2D_Control"
c2="2D_Mtb"

metadata <- read.csv("C:/Users/mr4e14/Documents/Rdocuments/R_MSstudy_sheets/Microsphere_sample_description_batch_2D.csv")
col.cell <- brewer.pal(n=6,"Set2")[metadata$Treatment]

sampleLabels <- metadata$Sample
sampleLabels

donor <- factor(metadata$Donor)
donor 
day <- factor(metadata$Day)
day
Group <- factor(metadata$Treatment)
Group

design <- model.matrix(~donor+Group)

design
ncol(design)

cbind(metadata,Group=Group)
#contrast <- makeContrasts(TB=GroupMtb,levels=design)
#contrast

colnames(design)
design

path <- file.path("C:/Users/mr4e14/Documents/REPEAT_MSstudy/REPEAT_MSstudy_kallisto_rawfastq_unstranded_cdna1_resultsONLY",metadata$Sample, "abundance.h5")
all(file.exists(path)) 
head(path)

# Import Kallisto transcript counts into R using Tximport ----
# copy the abundance files to the working directory and rename so that each sample has a unique name
Txi_gene <- tximport(path, 
                     type = "kallisto", 
                     tx2gene = Tx, 
                     txOut = FALSE, #How does the result change if this =FALSE vs =TRUE?
                     countsFromAbundance = "lengthScaledTPM",
                     ignoreTxVersion=TRUE)

# these are the counts after adjusting for transcript length
# Next, these are your transcript per million (TPM) values
head(Txi_gene$abundance) 
colSums(Txi_gene$counts)
barplot(colSums(Txi_gene$counts),main="MSstudy counts unfiltered", ylab="Counts", xlab="Sample", col=col.cell)

metadata$Sample
raw.data <- Txi_gene$counts
colnames(raw.data) <- metadata$Sample
dim(raw.data)
write.csv(raw.data,file="2D_Mtb_vs_2D_Control_tximport_results.unfiltered.csv")

# DEG ---------------------------------------------------------------------
FILTER=0.8

# iterative ---------------------------------------------------------------

DGEList <- DGEList(Txi_gene$counts,samples = metadata$Sample, genes=row.names(Txi_gene$counts))

table(rowSums(DGEList$counts==0)==12) #How many genes have no counts across all 12 samples?
cpm <- cpm(DGEList)
colnames(cpm) <- metadata$Sample
keepers <- rowSums(cpm>FILTER)>4
DGEList.filtered <- DGEList[keepers,]
colnames(DGEList.filtered) <- metadata$Sample
dim(DGEList.filtered)
counts_filtered<-DGEList$counts
DGEList.filtered.norm <- calcNormFactors(DGEList.filtered, method = "TMM")

#voomWithQualityWeights

v <- voomWithQualityWeights(DGEList.filtered.norm, design=design, normalization="none", plot=TRUE)
fit <- lmFit(v,design)
fit <- eBayes(fit)
All_results<-topTable(fit, coef=7, number=nrow(fit), adjust.method="BH", sort.by="P")

counts_per_M<-cpm(DGEList.filtered.norm)
counts_per_M <- as.data.frame((counts_per_M))
All_results <- merge(All_results,counts_per_M,by="row.names")


DE_results<-All_results[All_results[,"adj.P.Val"]<=0.05,]
dim(DE_results)

colnames(DGEList.filtered.norm) <- metadata$Sample
dim(DGEList.filtered.norm)

write.csv(All_results,file=paste(c2,"_vs_",c1,"_limma.All_results.","filter_",FILTER,"_TMM",".csv",sep=""),row.names=FALSE)
write.csv(DE_results,file=paste(c2,"_vs_",c1,"_limma.DE_results.","filter_",FILTER,"_TMM",".csv",sep=""),row.names=FALSE)

#Plot PCA
Y <- apply(cpm, 1, function(y) scale(y, center=TRUE, scale=FALSE))
s <- svd(Y)

pdf("PCA.pdf")
plot(s$u[,1], s$u[,2], pch=19, cex=4, xlab="PC1", ylab="PC2",main="PCA plot TMM Normalization",col=col.cell,
     ylim=c((min(s$u[,2])*1.2),(max(s$u[,2])*1.2)),xlim=c((min(s$u[,1])*1.2),(max(s$u[,1])*1.2)))
text(s$u[,1], s$u[,2], labels=colnames(cpm),adj=c(0,2),cex=1)
dev.off()

#Plot IQR vs Median TMM normalised
log_counts<-log(counts_filtered+1)

IQR<-apply(cpm, 2, IQR)
Median<-apply(cpm, 2, median)
diff1<-mean(Median)-min(Median)
diff2<-max(Median)-mean(Median)
diff3<-mean(IQR)-min(IQR)
diff4<-max(IQR)-mean(IQR)
Xlim=c(mean(Median)-2*diff1,mean(Median)+2*diff2)
Ylim=c(mean(IQR)-2*diff3,mean(IQR)+2*diff4)

plot(Median, IQR, main="IQR vs. Median TMM normalisation", type="n", xlim=Xlim,ylim=Ylim, col=col.cell)
text(Median, IQR, labels=names(IQR),col=col.cell, cex=0.5)

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









