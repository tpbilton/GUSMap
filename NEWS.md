# GUSMap 2.2.2

## Bugs fixed
* Bug in the `makeFS` function showing the number of SNPs filtered for PVALUE was incorrect. Now fixed
* Fixed some bugs in the BC object. Would have been copy and paste errors

# GUSMap 2.2.1

## Updates
* Major changes to the `makeFS` function. Added functionality for determining how many SNPs are discarded based on the different filtering criteria
* Modified how the genotypes in the parents are being called. Now allow for 5% error in calling a homozygous genotype (i.e., if the major/minor allele is between 0% and 5%, then the genotype call is homozygous for the prominent allele)

## Bugs fixed:
* Issue with `computeMap` when `inferSNPs=TRUE`. Now fixed.
* BIN filter was not being implemented correctly. Suspect this would have resulted in fewer SNPs being returned at the end of the process. Now fixed.

# GUSMap 2.2.0

## New Features 
* Added new class IC for analysis based on intercross (F2) populations.
* New vignette for the IC class
* Functionality for infering the halpotypes of the progeny from the final map. These haplotypes can be included in the file
generated by the `$writeLM` function.

## Updates
* Completed update of documentation.

# GUSMap 2.1.0

## New Features 
* Completed updated of new `BC` class for linkage mapping using a pseudo-testcross approach. Similar to the FS class but with some differences to handle back-cross types of analyses.
* Additional argument `inferSNPs` to the `makeFS` function. If `TRUE`, the segregation type of SNPs with unknown segregation is inferred from the progeny genotype information.
* New function `addSNPs` to the FS object for adding additional SNPs to the linkage groups (including the SI SNPs inferred when `inferSNPs = TRUE`)

## Updates
* Updated functions in FS object to handle the new versions of the cometPlot in the GUSbase package.
* Version bugs have also been fixed.
* Added additional filter to the `makeFS` function
* Documentation has been updated and a new vignette has been created for the new `BC` class.

# GUSMap 2.0.1

A few bugs were fixed in this version relating to inferring the segregation of SNPs when the read depth in the parents was insufficient.

# GUSMap 2.0.0

* This verison accompanies the Thesis by Bilton (2019).

## New Features

* Functions relating to a backcross population have been created.
   
   * `BC` class for performing linkage mapping in backcross populations analogous to the `FS` class
   * `makeBC` which creates a backcross population
   * Added functionality to infer segregation type in cases where the depth on the parental genotypes is too low.

## Improvements

* Optimization of the HMM using `optim` now using first derivatives to improve computational performance and can now handle sex-specific cases and/ro SNP specific error parameters
* The EM algorithm has also been generalized to handle sex-specific cases and/ro SNP specific error parameters

## Bugs

A number of bugs have been fixed from the previous release.

# GUSMap 1.0.0

This version of GUSMap has been completed restructured from the previous and now includes function for undertakingthe whole linkage mmapping process from filtering the data to creating linkage groups, ordering linakge groups and computing linkage maps.

## New Features

* GUSMap now depends on the new package [GUSbase](https://github.com/tpbilton/GUSbase) for reading in data from VCF files.
* GUSMap now structured around creating an FS (and R6) object from the `makeFS` function which stores the data in a safe way and contains various functions needed for performing the linkage mapping process. 
* Functions available for an FS object:

   * Creating linkage groups
   * Editing linkage groups (e.g., remove SNPs or linkage groups, merging linkage groups)
   * Ordering linkage groups using the MDS algorithm by Preedy and Hackett (2016).
   * Various plot function for visualizing the results (including some interactive heatmaps)
   * Writing results from a file

* Created a vignette as an introduction to GUSMap demonstrating how to use it.
  
## Improvements

* The computational time for the optimization of the likelihoods to compute the recombination fraction estimates and inferring the parental phase has been considerably reduced by:

   1. Computing the derivatives directly (instead of using finite differencing) in the optim routine.
   2. Implementing OpenMP parallalization in the computation of the likelihood.
   3. Optimizing some of the code.
   
## Changes:

* The `simFS` function for simulating full-sib families has been changed. It now returns an FS object with the simulated data.
* The `rf_est` and `infer_OPGP_FS` function are no long available directly to the user but are called from the `$computeMap` function in a FS object.
* The transformation functions have been moved to the GUSbase package

## Bug fixes:

* Inference of parental phase now allows for backcross approaches where there are no SNPs segregating in one of the parental meioses.

# GUSMap 0.1.1

* This verison accompanies the publication by Bilton et al. (2018).

## Added Features

* Function for converting VCF files into RA format. The VCF file is required to have some information regarding read depths.
* Reading in an RA file with a pedigree file.
* Inclusion of the EM algorithm (in addition to using optim) for optimizing the likelihoods.

# GUSMap 0.1.0

## README 

* This is the first release of GUSMap.
* As this verison is only an early verison of the package, there are only a few features at present (See Features). More features are intended to be added on as the package is developed. In particular, the ability to group SNPs and order the SNPs within groups is to be included at a later stage. Suggestions for particular improvements are welcome.

## Features

* Phase inference on a set of ordered SNPs for a full-sib family.
* Recombination fraction estimation for an order set of SNPs using sequencing data without requiring the need to filter with respect to read depth for outcrossed full-sib families.
* Simulation of sequencing data for a full-sib family and functions to convert data into various other software formats.

# References

* Bilton, T.P., Schofield, M.R., Black, M.A., Chagne, D., Wilcox, P.L., Dodds, K.G. (2018). "Accounting for Errors in Low coverage high-throughput sequencing data when constructing genetic maps using biparental outcrossed populations." *Genetics*, *209*(1), 65--76. doi: [10.1534/genetics.117.300627](http://www.genetics.org/content/209/1/65)
* Bilton, T.P. (2019) Developing statistical methods for genetic analysis of genotypes from genotyping-by-sequencing data. (Unpublished doctoral dissertation). University of Otago, Dunedin, New Zealand.
* Preedy, K.F., Hackett, C.A. (2016). A rapid marker ordering approach for high-density genetic linkage maps in experimental autotetraploid populations using multidimensional scaling. *Theoretical and Applied Genetics*, *129*(11), 2117-2132. doi: [10.1007/s00122-016-2761-8](https://link.springer.com/article/10.1007%2Fs00122-016-2761-8) 
