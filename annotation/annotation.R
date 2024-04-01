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

but <- read.csv("file_name_chip_report.csv", header=TRUE)
head(but)
butgr <- makeGRangesFromDataFrame(but, keep.extra.columns=TRUE)
butgr

but.anno <- annotatePeakInBatch(butgr, AnnotationData=annoData, output="both", maxgap=1000L, select="all")

but.anno$gene_name <- annoData$gene_name[match(but.anno$feature, names(annoData))]

#k9ac2.anno
#k9ac2.anno.1kb <- k9ac2.anno[which(k9ac2.anno$shortestDistance <= 1000),]
#head(prop.anno.1kb$feature)
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
diffbut <- data.frame(seqnames=seqnames(but.anno),
  starts=(start(but.anno)-1),
  ends=end(but.anno),
  ENSG=but.anno$feature,
  ENSG.strand=but.anno$feature_strand, 
  insideFeature=but.anno$insideFeature, 
  distancetoFeature=but.anno$distancetoFeature,
  shortestDistance=but.anno$shortestDistance,
  symbol=but.anno$symbol, 
  entrezID=but.anno$entrezID,
  fold=but.anno$Fold,
  conc=but.anno$Conc,
  conc_0mM=but.anno$Conc_0.mM,
  conc_1mM=but.anno$Conc_1.mM,
  FDR=but.anno$FDR)
head(diffbut)
write.table(diffbut, file="H4K12bu.anno.bed", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
dim(diffbut)


# do pathway analysis in R
BiocManager::install("ReactomePA")
library(ReactomePA)
manual(ReactomePA)
head(dffbut$entrezID)

dffup <- dffbut[which(dffbut$fold < 0), ]
dim(dffup)

dffdown <- dffbut[which(dffbut$fold > 0), ]
dim(dffdown)

de <- (dffbut$entrezID)
length(de)
de <- de[!is.na(de)]
head(de)
length(de)
x <- enrichPathway(gene=de,pvalueCutoff=0.2, readable=TRUE)
x
x$Description

class(x)
dev.off()
barplot(x, showCategory=8)
dotplot(x, showCategory=8)
?dotplot

deup <- dffup$entrezID
deup <- deup[!is.na(deup)]
length(deup)
xup <- enrichPathway(gene=deup,pvalueCutoff=0.05, readable=TRUE)
xup  # 77 pathways found @0.2, 1 pathway found @0.05
barplot(xup, showCategory=8)
?barplot
dotplot(xup, showCategory=15)
xup$Description


dedown <- dffdown$entrezID
dedown <- dedown[!is.na(dedown)]
length(dedown)
xdown <- enrichPathway(gene=dedown, pvalueCutoff=0.2, readable=TRUE)
xdown  # no enriched pathways
dotplot(xdown, showCategory=15)
xdown$Description

