# GUSMap Version 1.0.0

[![gplv3+](https://img.shields.io/badge/license-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl.html)

Genotyping Uncertainty with Sequencing data and linkage MAPping (GUSMap).

An R package for performing linkage mapping using low and/or high coverage sequencing data without requiring filtering with respect to read depth.


### Installation:

The easiest way to install GUSMap in R is using the devtools package.

```
install.packages("devtools")
library(devtools)
install_github("tpbilton/GUSbase")
install_github("tpbilton/GUSMap", ref="Development")
```

Note: Some of the functions are coded in C and therefore an appropriate C compiler is needed for the package to work. For windows OS, Rtools (https://cran.r-project.org/bin/windows/Rtools/) provides a compiler. 

### Tutorials:

The best way to understand how to sue GUSMap to construct linkage maps is to read through the [Introduction](inst/doc/Introduction.html) Tutorial contained in the package.

### Development

There are plans to contiune extending this package. If you have any suggests for improvement or additional features you like to see, I'd suggest posting a question under the Issues. I'm also happy for people to suggest changes via a pull request provided it fits within the favour of the package.


### Citation:

To cite this R package:

Bilton, T.P., Schofield, M.R., Black, M.A., Chagn&#233;, D., Wilcox, P.L., & Dodds, K.G. (2018). Accounting for errors in low coverage high-throughput sequencing data when constructing genetic maps using biparental outcrossed populations. *Genetics*, *209*(1), 65--76. doi:[10.1534/genetics.117.300627](http://www.genetics.org/content/209/1/65) 