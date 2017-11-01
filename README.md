# GusMap Version 0.1.0

Genotyping Uncertainty with Sequencing data and linkage MAPping (GUSMap).

An R package for constructing genetic linkage maps using low (or high) coverage sequencing data without requiring filtering with respect to read depth. This release accompanies the publication by Bilton et al. (2017).

### Installation:

The easiest way to install GUSMap in R is using the devtools package.

```
install.packages(devtools)
library(devtools)
install_github("tpbilton/GUSMap")
```

Note: Some of the functions are coded in C and therefore an appropriate C complier is needed for the package to work. For windows OS, Rtools (https://cran.r-project.org/bin/windows/Rtools/) provides a complier. 

### Data Format:

The key data required to use this package are:

- Matrix of allele counts for the reference allele. The rows are the individuals and columns are the SNPs. Entries must be a positive integer value.
- Matrix of allele counts for the alternate allele. The rows are the individuals and columns are the SNPs. Entries must be a positive integer value.
- Config Vector: This vector gives the segregation type for each SNP and the length of the vector must equal the number of SNPs. Enteries are 1 for both-informative SNP (ABxAB), 2 for paternal-informative SNP (ABxAA), 2 for paternal-informative SNP (ABxBB), 4 for a maternal-informative SNP (AAxAB), and 5 for a maternal-informative SNP (BBxAB). 

### References:

Bilton, T.P., Schofield, M.R., Black, M.A., Chagne D., Wilcox P., Dodds K.G. (2017). Accounting for errors in low coverage sequencing data when constructing genetic maps using biparental outcrossed populations. Unpublished manuscript.

