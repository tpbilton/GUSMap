# Optimise call to rfEst for smaller case

Run code and `save` object to file just before the call, so we can just focus
on improving that section of code.

## Adding OpenMP

Loops in the following C functions were parallelised with OpenMP:

* `EM_HMM`
* `EM_HMM_UP`
* `ll_fs_scaled_err_c`
* `score_fs_scaled_err_c`

After making these changes the code was benchmarked:

![benchmark results](gusmap-smaller-openmp-compare.png)
