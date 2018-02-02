# GUSMap Version 0.1.0

Genotyping Uncertainty with Sequencing data and linkage MAPping (GUSMap).

An R package for constructing genetic linkage maps using low and/or high coverage sequencing data without requiring filtering with respect to read depth. It accounts for errors associated with low sequencing depth and miscalled bases. This release accompanies the paper by Bilton et al. (2017).

[![gplv3+](https://img.shields.io/badge/license-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl.html)

### Installation:

The easiest way to install GUSMap in R is using the devtools package.

```
install.packages(devtools)
library(devtools)
install_github("tpbilton/GUSMap")
```

Note: Some of the functions are coded in C and therefore an appropriate C compiler is needed for the package to work. For windows OS, Rtools (https://cran.r-project.org/bin/windows/Rtools/) provides a compiler. 

### Reading Data into GUSMap:

The format of the data required by GUSMap is a matrix of allele counts for the reference and alternate allele and segregation type of each SNP. This data can
produced from a VCF in GUSMap using a couple of functions. There two steps to the process.

1. The VCF file needs to be converted into RA (reference/alternate) format using the function VCFtoRA(). The following example (using the chromosome 11 SNPs of the Manuka dataset used
by Bilton et al. (2017)) shows how this is done in GUSMap.
```
vcffile <- Manuka11()
RAfile <- VCFtoRA(vcffile)
```
The function Manuka() just returns the path to the chromosome 11 SNPs in the Manuka data set that comes with the package.
The function VCFtoRA() creates two files, an RA file and a partially completed pedigree file (we shall touch on the later in the next step). More information about the VCFtoRA() function and 
RA format data is available using the following command.
```
?VCFtoRA
```

2. The pedigree file created from VCFtoRA() is in the form of

|SampleID |IndividualID | Mother | Father | Family |
| ------- | ----------- | ------ | ------ | ------ |
|P1_S1    |             |        |        |        |
|P1_S2    |             |        |        |        |
|P2_S1    |             |        |        |        |
|P2_S2    |             |        |        |        | 
|O1       |             |        |        |        |
|O2       |             |        |        |        |
| ...     |             |        |        |        |
|O177     |             |        |        |        |

The SampleID column is the sample ID used in the VCF to represent each sample. For the Manuka data set, there are 177 offspring (which are denoted O1 to O177), 
while the samples P1_S1 and P1_S2 are two samples of one parent and P2_S1 and P2_S2 are two samples of the other parent. The IndividualID column is used to 
denote the individual to which each sample belongs, where each individual must have a unique identifier. This is useful as it allows us to have multiple samples 
on the same individual. The Mother column is used to specify which individual is the mother of that individual (and likewise the father column). Note that if the 
Mother or Father is unknown, then the entry in these columns should be left blank. The Family column is just used to allow the user to specify a name for each specific family.

For our Manuka dataset, a completed pedigree file should look like the following:

|SampleID | IndividualID | Mother | Father | Family |
| ------- | ------------ | ------ | ------ | ------ |
|P1_S1    |  1           |        |        |        |
|P1_S2    |  1           |        |        |        |
|P2_S1    |  2           |        |        |        | 
|P2_S2    |  2           |        |        |        |
|O1       |  3           | 2      | 1      | Manuka | 
|O2       |  4           | 2      | 1      | Manuka | 
| ...     | ...          | ...    | ...    | ...    |
|O177     | 179          | 2      | 1      | Manuka | 

3. The final step is to read in and process the RA data. For the Manuka example, the relavent code is the following:
```
MKdata <- readRA(RAfile, gform = "reference", pedfile=paste0(strsplit(RAfile,split="\\.")[[1]][1],"_ped.csv"))
```
The function readRA() processes the RA data file and extract list with the required data.
```
str(MKdata)
```
From the output, one can see that the object MKdata is a list with various components, which can be accessed as follows.
```
MKdata$depth_Ref
```

### Analysis in GUSMap:

To show how GUSMap is used to estimate the adjacent recombination fractions, we use the low depth (average read depth below 6) chromosome 11 from the Manuka data set in Bilton et al. (2017). First, it was noted that there were a number of SNPs that were erroreously placed on the chromsome. The following produces a vector with the  position of these erroreous SNPs.
```
badSnps <- rep(FALSE, ncol(MKdata$genon[[1]])) 
badSnps[c(3:14,16,18:23,25:27,29,33:35,39:42,44,46:49, 
   53,58:60,63,65:67,69:74,76,77,85,104,114:129, 
   152,162,163,173:180,187:190,206:208,210,212,213, 
   216,228,272,285,286,298,304,314,318,320,342,344,345, 
   354,358,370,373,380,389,390,397,408,411,421,427, 
   437,446,469,493,498:500,506,510,516,519,522,531, 
   543,553,556,569,582,585,601,617,620,623,628,635, 
   636:643,647:650,669:680)] <- TRUE 
```
We now perform the analysis using the low depth depth as in Bilton et al. (2017). To extract these SNPs, we use the code,
```
ind_rd6 <- which(colMeans(MKdata$depth_Ref[[1]] + MKdata$depth_Alt[[1]],na.rm=T ) < 6 & !badSnps) 
```
Given the order of these SNPs, the parental phase (or OPGP) can be infered using the function infer_OPGP_FS.
```
OPGP_rd6 <- infer_OPGP_FS(MKdata$depth_Ref[[1]][,ind_rd6], MKdata$depth_Alt[[1]][,ind_rd6], 
                          MKdata$config[[1]][ind_rd6], epsilon=0.001, reltol=1e-3) 
```
Once the OPGPs have been correctly inferred, the adjacent recombination fractions can be estimated using the function rf_est_FS().
```
rf_est <- rf_est_FS(depth_Ref=list(MKdata$depth_Ref[[1]][,ind_rd6]), 
                    depth_Alt=list(MKdata$depth_Alt[[1]][,ind_rd6]), 
                    OPGP=list(OPGP_rd6)) 
```
The recombination fraction estimates can be extract using,
```
rf_est$rf
```
and the estimate of the sequencing error rate can be extract using,
```
rf_est$epsilon
```
while the value of the log-likelihood can be extract using,
```
rf_est$loglik
```
Once the recombination fractions have been estimated, they can be written to a file as follows.
```
write(rf_est$rf, file="Manuka_chr11_rd6_rf.txt", ncolumns=1)
```

### References:

Bilton, T.P., Schofield, M.R., Black, M.A., Chagn&#233, D., Wilcox, P.L., Dodds, K.G. (2017). Accounting for errors in low coverage high-throughput sequencing data when constructing genetic maps using biparental outcrossed populations. Unpublished manuscript.

### License

GUSMap is Copyright 2017 Timothy P. Bilton, released under the GNU General Public License version 3.

