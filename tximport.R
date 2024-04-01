library(tximport)
library(readr)

# locate the directory containing the files
dir <- system.file("extdata", package = "tximportData")
list.files(dir)

# vector pointing to the quantification files

samples <- read.table(file.path(dir, "samples.txt"), header = TRUE)
samples

files <- file.path(dir, "rsem", "genes", paste0(samples$sample, ".genes.results"))
files


names(files) <- paste0("sample", 1:6)
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE, importer=read.delim)
head(txi.rsem$counts)


write.csv(txi.rsem, 'tximport.csv')
