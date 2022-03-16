library(tximport)
library(edgeR)
library(limma)
library(ggplot2)
library(dplyr)
library(tidyverse) 
library(ensembldb)
library(EnsDb.Hsapiens.v86) 

### Anotations table 
listTables(EnsDb.Hsapiens.v86)
listColumns(EnsDb.Hsapiens.v86, "tx")
Tx <- transcripts(EnsDb.Hsapiens.v86, columns=c(listColumns(EnsDb.Hsapiens.v86,"tx"), "gene_name"))
Tx <- as_tibble(Tx)
Tx <- dplyr::rename(Tx, target_id = tx_id)
Tx <- dplyr::select(Tx, target_id, gene_name)
head(Tx)

###Load metadata, detailed on GEO (GSE174566)
metadata <- read.csv("md.csv")

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
contrast <- makeContrasts(TB=GroupMtb,levels=design)
contrast

colnames(design)
design

path <- file.path("rawfastq",metadata$Sample, "abundance.h5")
all(file.exists(path)) 
head(path)

# Import Kallisto transcript counts into R using Tximport
Txi_gene <- tximport(path, 
                     type = "kallisto", 
                     tx2gene = Tx, 
                     txOut = FALSE, #How does the result change if this =FALSE vs =TRUE?
                     countsFromAbundance = "lengthScaledTPM",
                     ignoreTxVersion=TRUE)

head(Txi_gene$abundance) 
colSums(Txi_gene$counts)

metadata$Sample
raw.data <- Txi_gene$counts
colnames(raw.data) <- metadata$Sample
dim(raw.data)
write.csv(raw.data,file="2D_Mtb_vs_2D_Control.csv")

# DGE
FILTER=0.8

DGEList <- DGEList(Txi_gene$counts,samples = metadata$Sample, genes=row.names(Txi_gene$counts))

table(rowSums(DGEList$counts==0)==12) 
cpm <- cpm(DGEList)
colnames(cpm) <- metadata$Sample
keepers <- rowSums(cpm>FILTER)>4
DGEList.filtered <- DGEList[keepers,]
colnames(DGEList.filtered) <- metadata$Sample
dim(DGEList.filtered)
counts_filtered<-DGEList$counts
DGEList.filtered.norm <- calcNormFactors(DGEList.filtered, method = "TMM")

#voomWithQualityWeights
v <- voomWithQualityWeights(DGEList.filtered.norm, design=design, normalize.method="none", plot=TRUE)
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

write.csv(All_results,file=paste("MSstudy_2D_Mtb_vs_Control_All_results.csv"),row.names=FALSE)
write.csv(DE_results,file=paste("MSstudy_2D_Mtb_vs_Control_DE_results.csv"),row.names=FALSE)
