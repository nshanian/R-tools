#set directory
setwd("~/<path_to_working_directory>/")
#source("https://bioconductor.org/biocLite.R")
#annotate genomic features within 5kb
#biocLite("EnsDb.Hsapiens.v86")
library(EnsDb.Hsapiens.v86)
library("AnnotationDbi")
#biocLite("org.Hs.eg.db")
library("org.Hs.eg.db")
#biocLite("ChIPpeakAnno")
library(ChIPpeakAnno)
#biocLite(GenomicRanges)
library(GenomicRanges)

annoData <- annoGR(EnsDb.Hsapiens.v86, feature="gene")
info(annoData)
annoData[1:2]

prop <- read.csv("file_name_report.csv", header=TRUE)
head(prop)
propgr <- makeGRangesFromDataFrame(prop, keep.extra.columns=TRUE)
propgr

prop.anno <- annotatePeakInBatch(propgr, AnnotationData=annoData, output="both", maxgap=1000L, select="all")

prop.anno$gene_name <- annoData$gene_name[match(prop.anno$feature, names(annoData))]


head(prop.anno$feature)
prop.anno$symbol <- mapIds(org.Hs.eg.db, 
                          keys=prop.anno$feature,
                          column="SYMBOL",
                          keytype="ENSEMBL",
                          multiVals="first")
prop.anno

prop.anno$entrezID <- mapIds(org.Hs.eg.db, 
                            keys=prop.anno$feature,
                            column="ENTREZID",
                            keytype="ENSEMBL",
                            multiVals="first")
prop.anno


#profile of ChIP peaks binding to TSS regions
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ggupset)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix <- getTagMatrix(propgr, windows=promoter)
tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red")

plotAvgProf(tagMatrix, xlim=c(-3000, 3000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")

plotAvgProf(tagMatrix, xlim=c(-3000, 3000), conf = 0.95, resample = 1000)

peakAnno <- annotatePeak(propgr, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")

plotAnnoPie(peakAnno)
plotAnnoBar(peakAnno)
upsetplot(peakAnno)

