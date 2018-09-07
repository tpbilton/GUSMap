
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
* Preedy, K.F., Hackett, C.A. (2016). A rapid marker ordering approach for high-density genetic linkage maps in experimental autotetraploid populations using multidimensional scaling. *Theoretical and Applied Genetics*, *129*(11), 2117-2132. doi: [10.1007/s00122-016-2761-8](https://link.springer.com/article/10.1007%2Fs00122-016-2761-8) 
