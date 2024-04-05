## R data analysis tools

This repository R tools for statisitcal analysis and visualization of next-generation sequencing (NGS) data such as ChIP-seq, ATAC-seq and RNA-seq.

The R scripts include the following analyses: 

 - ChIP-seq data visualization (`ChIPseeker.R`)
 - Genomic peak annotation (`annotation.R`, `narrowPeak_anno.R`) 
 - Accessing Ensembl annotation (`biomaRt`)
 - Import and summarize transcript-level abundance (`tximport.R`)
 - Data transformation and visualization (`tidyverse`, `ggplot`)

### Documentation and References:

EnsDb.Hsapiens.v86: https://bioconductor.org/packages/release/data/annotation/html/EnsDb.Hsapiens.v86.html

AnnotationDbi: https://www.bioconductor.org/packages//2.10/bioc/html/AnnotationDbi.html

org.Hs.eg.db: https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html

GenomicRanges: https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html

ChIPpeakAnno: https://bioconductor.org/packages/release/bioc/html/ChIPpeakAnno.html

ChIPseeker: https://bioconductor.org/packages/release/bioc/html/ChIPseeker.html

biomaRt: https://bioconductor.org/packages/release/bioc/html/biomaRt.html

tximport: https://bioconductor.org/packages/release/bioc/html/tximport.html

ggplot2: https://cran.r-project.org/web/packages/ggplot2/index.html

ggplot2 tidyverse: https://ggplot2.tidyverse.org

tidyverse: https://www.tidyverse.org/

The scripts `run_spp.R` and `run_spp_nudups.R` are by:

[1] Anshul Kundaje, Computer Science Dept., MIT, ENCODE Consortium, http://code.google.com/p/phantompeakqualtools, Feb 2013

[2] Kharchenko PK, Tolstorukov MY, Park PJ, Design and analysis of ChIP-seq experiments for DNA-binding proteins Nat Biotechnol. 2008

Please see README.txt for more information.

