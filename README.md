# GusMap Version 0.1

Genotyping Uncertainty with Sequencing data and linkage MAPping (GusMap).

An R package for performing linkage mapping using low (or high) coverage sequencing data with requiring filtering with respect to read depth. This release accompanies the publication by Bilton et al. (2017).  

### Installation:

The easiest way to install GusMap in R is using the devtools package.

```
install.packages(devtools)
library(devtools)
install_github("tpbilton/GusMap")
```

Note: Some of the functions are coded in C and therefore an appropriate C complier is needed for the package to work. For windows OS, Rtools (https://cran.r-project.org/bin/windows/Rtools/) provides a complier. 

### References:

Bilton, T.P., Schofield, M.R., Black, M.A., Chagne D., Wilcox P., Dodds K.G. (2017). Multilocus genetic linkage mapping in multiparental outcrossed populations using low coverage sequencing data. Unpublished manuscript.

