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
vcffile <- Manuka()
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
|P1_S1    |             |        |        |        |
|P1_S2    |             |        |        |        |
|P2_S1    |             |        |        |        |
|P2_S2    |             |        |        |        | 
|O1       |             |        |        |        |
|O2       |             |        |        |        |  
| ...     |             |        |        |        |
|O177     |             |        |        |        |

The SampleID column is the sample ID used in the VCF to represent each sample. For the Manuka data set, there are 177 offspring (which are denoted O1 to O177) 
while the samples P1_S1 and P1_S2 are two samples of one parent and P2_S1 and P2_S2 are two samples of the other parent. The IndividualID column is used to 
denote the individual to which each sample belongs, where each individual must have a unique identifier. This is useful as it allows us to have multiple samples 
on the same individual. The Mother column is used to specify which individual is the mother of that individual (and likewise the father column). Note that if the 
Mother or Father is unknown, then the entry in these columns should be left blank. The Family column is just used to allow the user to specify a name for each specific family.

For our Manuka dataset, a completed pedigree file should look like the following:

|SampleID | IndividualID | Mother | Father | Family |
| ------- | ------------ | ------ | ------ | ------ |
|P1_S1	  |  1           |        |        |        |
|P1_S2    |  1           |        |        |        |
|P2_S1    |  2           |        |        |        | 
|P2_S2    |  2           |        |        |        |  
|O1       |  3           | 2      | 1      | Manuka | 
|O2       |  4           | 2      | 1      | Manuka | 
| ...     | ...          | ...    | ...    | ...    |
|O177     | 179          | 2      | 1      | Manuka | 

3. The final step is to read in and process the RA data. For the Manuka example, the relavent code is the following:
```
MKdata <- readRA(RAfile, gform = "reference", pedfile=paste0(strsplit(RAfile,split=".")[[1]][1],"_ped.vcf"))
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

### References:

Bilton, T.P., Schofield, M.R., Black, M.A., Chagn&#233;, D., Wilcox, P.L., Dodds, K.G. (2017). Accounting for errors in low coverage high-throughput sequencing data when constructing genetic maps using biparental outcrossed populations. Unpublished manuscript.

### License

GUSMap is Copyright 2017 Timothy P. Bilton, released under the GNU General Public License version 3.

