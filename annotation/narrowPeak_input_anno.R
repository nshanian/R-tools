#set directory
setwd("~/Desktop/SCFA data/ChIP-seq/SW480/input_butyr_GS")
#source("https://bioconductor.org/biocLite.R")
#annotate genomic features within 5kb
#biocLite("EnsDb.Hsapiens.v86")
library(EnsDb.Hsapiens.v86)
library(AnnotationDbi)
#biocLite(org.Hs.eg.db)
library(org.Hs.eg.db)
#biocLite(ChIPpeakAnno)
library(ChIPpeakAnno)
#biocLite(GenomicRanges)
library(GenomicRanges)
library(clusterProfiler)
#biocLite(rtracklayer)
library(rtracklayer)

annoData <- annoGR(EnsDb.Hsapiens.v86, feature="gene")
info(annoData)
annoData[1:2]

but <- read.csv("but_input_GS.narrowPeak.csv", header=TRUE)
head(but)

butgr <- makeGRangesFromDataFrame(but, keep.extra.columns=TRUE)
head(butgr)

covplot(butgr, weightCol="signalValue")

but.anno <- annotatePeakInBatch(butgr, AnnotationData=annoData, output="both", maxgap=1000L, select="all")

but.anno$gene_name <- annoData$gene_name[match(but.anno$feature, names(annoData))]
head(but.anno$feature)

but.anno$symbol <- mapIds(org.Hs.eg.db, 
                          keys=but.anno$feature,
                          column="SYMBOL",
                          keytype="ENSEMBL",
                          multiVals="first")
but.anno

but.anno$entrezID <- mapIds(org.Hs.eg.db, 
                            keys=but.anno$feature,
                            column="ENTREZID",
                            keytype="ENSEMBL",
                            multiVals="first")
but.anno

# write annotated output
write.table(but.anno, file="but_input_GS.anno.bed", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
head(but.anno)


#profile of ChIP peaks binding to TSS regions
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ggupset)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix <- getTagMatrix(butgr, windows=promoter)
tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red")

plotAvgProf(tagMatrix, xlim=c(-3000, 3000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")

plotAvgProf(tagMatrix, xlim=c(-3000, 3000), conf = 0.95, resample = 1000)

peakAnno <- annotatePeak(butgr, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")

plotAnnoPie(peakAnno)
plotAnnoBar(peakAnno)
upsetplot(peakAnno)


