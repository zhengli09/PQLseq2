# PQLseq2: A Faster Re-implementation of PQLseq

PQLseq implements a Poisson and Binomial mixed model for heritability estimation and differential analysis of RNA sequencing (RNAseq) data and Bisulfite sequencing (BSseq) data in the presence of individual relatedness and population structure. PQLseq2 provides a faster re-implementation of PQLseq. By fully implementing its core algorithm in C++, PQLseq2 circumvents repetitive type casting and object duplications in PQLseq, achieving a 15-fold increase in computational speed while providing consistent estimates for all the parameters as compared to PQLseq.

## Installation

``` r
if(!require(devtools))
  install.packages(devtools)
devtools::install_github("zhengli09/PQLseq2")
library(PQLseq2)
?pqlseq2
```

## Refer to `PQLseq` for detailed algorithm

Shiquan Sun*, Jiaqiang Zhu*, Sahar Mozaffari, Carole Ober, Mengjie Chen, and Xiang Zhou#. *Heritability estimation and differential analysis of count data with generalized linear mixed models in genomic sequencing studies.* Bioinformatics 35, no. 3 (2019): 487-496.
