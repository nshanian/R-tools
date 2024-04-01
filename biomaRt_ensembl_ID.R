#set directory
setwd("~/Desktop/SCFA data/RNA-seq/SW480/prop/")

library(biomaRt)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

ensembl_gene_id_version <-read.csv("SW480_prop_ENS_ID.csv", header=FALSE)
head(ensembl_gene_id_version)


library(stringr)

df <- subset(data.frame(str_split_fixed(ensembl_gene_id_version$V1, "\\.",2)), select = -c(X2))
colnames(df) <- c('gene_id')
head(df)

write.csv(df, file = 'SW480_prop_ENS_ID_SPLIT_NAMES.csv', row.names = FALSE)
read.csv('SW480_prop_ENS_ID_SPLIT_NAMES.csv')



library(dplyr)

nth(strsplit(ensembl_gene_id_version, split ="\\."))

getBM(attributes = ensemble_ID,
      filters = ensemble_ID,
      values = ensemble_ID,
      mart = mart)


gene_ids <- str_replace(ensemble_ID,
                        pattern = ".[0-9]+$",
                        replacement = "")

getBM(attributes = c('ensembl_gene_id'),
      filters = 'ensembl_gene_id', 
      values = gene_ids,
      mart = mart)

#####
ensembl_gene_id_version
strsplit(ensembl_gene_id_version$ensembl_gene_id_version, '\\.')
