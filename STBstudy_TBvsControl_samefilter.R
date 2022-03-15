library(limma)
library(ggplot2)
library(dplyr)
library(edgeR)
library(RUVSeq)
library(RColorBrewer)
library(rgl)
library(EDASeq)

###Load data and metadata
## data is raw counts, derived from tximport, annotated using EnsDb.Hsapiens.v86
## data and metadata located on GEO (GSE174443)

data <- read.csv('rawcounts.csv',row.names = 1)
metadata <- read.csv('md.csv')
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
write.csv(topTable, file=paste("All_results","filter_",FILTER,"_TMM",".csv",sep=""))

summary(decideTests(efit))

