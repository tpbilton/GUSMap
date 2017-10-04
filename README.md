# GusMap Version 0.1.0

Genotyping Uncertainty with Sequencing data and linkage MAPping (GUSMap).

An R package for performing linkage mapping using low (or high) coverage sequencing data with requiring filtering with respect to read depth. This release accompanies the publication by Bilton et al. (2017).

### Installation:

The easiest way to install GUSMap in R is using the devtools package.

```
install.packages(devtools)
library(devtools)
install_github("tpbilton/GUSMap")
```

Note: Some of the functions are coded in C and therefore an appropriate C complier is needed for the package to work. For windows OS, Rtools (https://cran.r-project.org/bin/windows/Rtools/) provides a complier. 

### Data Format:

The format of the data required are:

- Genon Matrix: Rows are the individuals and columns are the SNPs. Entries of the matrix must be an integer value between 0 and 2 corresponding to the number of reference alleles called and NA for missing genotypes.
- Depth Matrix: Rows are the individuals and columns are the SNPs. Entries must be a positive integer value. 
- Config Vector: This vector gives the segregation type for each SNP and the length of the vector must each the number of SNPs. Enteries are 1 for both-informative SNP, 2 for paternal-informative SNP, 3 for maternal-informative SNP and 4 for a uniformative SNP. 

### References:

Bilton, T.P., Schofield, M.R., Black, M.A., Chagne D., Wilcox P., Dodds K.G. (2017). Comstructing accurate genetic maps using low coverage sequencing data in multiparental outcrossed populations. Unpublished manuscript.

