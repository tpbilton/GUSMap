# Adding OpenMP

## General

Full changes can be viewed
[here](https://github.com/chrisdjscott/GUSMap/compare/cd3b05ed9321235eb42d31a677538e9acc61c06f...6053b15a4c89b98ff867cdd8bbbdd58797bf1ba5).

*Note*: ignore the changes about `doSNOW` and `doParallel`, this was just to
get around a problem with Pan and will be changed back. Also, there are extra
print statements that will be removed.

### Setting the number of threads in the C code

The number of threads can be set by passing the `nThreads` argument. Currently,
if this is not set it will default to using all available OpenMP cores. This
will either be the number of cores detected by the operating system, or the
value of the `OPENMP_NUM_THREADS` environment variable, if it is set.

To handle this, in the code we need to call an OpenMP function to get the
number of available cores and therefore need to include the OpenMP header
file. In case the code is compiled without OpenMP support, we add a default,
serial version of this function too:

```c
#ifdef _OPENMP
    #include <omp.h>
#else
    inline int omp_get_max_threads() { return 1; }
#endif

SEXP my_function(...) {
    ...

    // if nThreads is set to zero then use everything
    if (nThreads_c <= 0) {
        nThreads_c = omp_get_max_threads();
    }

    ...
}
```

Note the `#ifdef _OPENMP` - the code will only get in here if the compiler
supports OpenMP and OpenMP is enabled (e.g. in the `Makevars` file).

### Enabling OpenMP at compile time

In the `src/Makevars` file we add the following lines to enable OpenMP:

```
PKG_LIBS = $(SHLIB_OPENMP_CFLAGS)
PKG_CFLAGS = $(SHLIB_OPENMP_CFLAGS)
```

### Calling `rf_est` from R

Added `nThreads` argument to `rf_est` call.

* default (if `nThreads` is not specified) is to use either all cores on the system or, if the `OMP_NUM_THREADS` environment variable is set, it will use that number:
  
  ```r
  MK_fs$rf_est(method="EM")
  MK_fs$rf_est(method="EM", nThreads=0)
  ```

* otherwise, set `nThreads` to the number of threads you want to use:

  ```r
  # use 4 threads
  MK_fs$rf_est(method="EM", nThreads=4)
  ```

## em.c

[Link to file](https://github.com/chrisdjscott/GUSMap/blob/opt_em/src/em.c)

### computeProb

We had to change the way the inputs were passed into this function. Originally
pointers were passed in and incremented during the loop to access the next
data element. This only works if the loop is run in serial. In order to run in
parallel, we pass in the array itself, which requires us to add the dimensions
as an argument too.

Instead of:

```c
double computeProb(double *ppAA, double *ppBB, ...) {
    ...
    for (ind = 0; ind < nInd; ind++) {
        int snp;
        for(snp = 0; snp < nSnps; snp++) {
            ...
            // the loops must run sequentially in order to access the correct element
            *ppAA = ...;
            ppAA++;
            ...
        }
    }
    ...
}
```

we end up with:

```c
double computeProb(int nInd, int nSnps, double pAA[nInd][nSnps], double pBB[nInd][nSnps], ...) {
    ...
    #pragma omp parallel for num_threads(nThreads_c)
    for (ind = 0; ind < nInd; ind++) {
        int snp;
        for(snp = 0; snp < nSnps; snp++) {
            ...
            // now the loops can be run in any order (or in parallel)
            pAA[ind][snp] = ...;
            ...
        }
    }
    ...
}
```

It is now safe to run the loop in parallel.


### EM_HMM

Added timings around the loops to find out which ones were worthwhile to convert.

```c
// Compute the forward and backward probabilities for each individual
for(fam = 0; fam < noFam_c; fam++) {
    // reduction for llval and mark variables as private
    #pragma omp parallel for reduction(+:llval) num_threads(nThreads_c) \
                         private(sum, s1, s2, alphaDot, snp, w_new, betaDot)
    for(ind = 0; ind < nInd_c[fam]; ind++) {
        // this variable is implicitly private
        int indx = ind + indSum[fam];
        ...
        // each variable does its own sum, with reduction performed at the end
        llval = llval + log(sum);
        ...
    }
}
```

Also converted some of the later loops, over the number of SNPs, for example:

```c
// Error parameter:
if(seqError_c){
    sumA = 0;
    sumB = 0;
    // reduction for multiple variables
    #pragma omp parallel for reduction(+:sumA,sumB) \
                             private(fam, ind, a, b, s1) \
                             num_threads(nThreads_c)
    for(snp = 0; snp < nSnps_c; snp++){
        for(fam = 0; fam < noFam_c; fam++){
            for(ind = 0; ind < nInd_c[fam]; ind++){
                int indx = ind + indSum[fam];
                ...
                sumA = sumA + ...;
                sumB = sumB + ...;
            }
        }
    }
}
```


### EM_HMM_UP

Very similar to `EM_HMM`.


## likelihoods.c

[Link to file](https://github.com/chrisdjscott/GUSMap/blob/opt_em/src/likelihoods.c)

### ll_fs_scaled_err_c

There is one main loop in this function and we parallelise that one:

```c
// Now compute the likelihood
#pragma omp parallel for reduction(+:llval) num_threads(nThreads_c) \
                         private(sum, s1, alphaDot, alphaTilde, snp, s2, w_new)
for(ind = 0; ind < nInd_c; ind++) {
    ...
    llval = llval + log(sum);
    ...
}
```


## score.c

[Link to file](https://github.com/chrisdjscott/GUSMap/blob/opt_em/src/score.c)

### score_fs_scaled_err_c

Similar to likelihood function, we have the one main loop but also compute
the gradient vector.

More complicated because different iterations of the main loop may write to
the same element of the gradient vector, which requires special handling when
converting to run in parallel.

```c
// Now compute the likelihood and score function
#pragma omp parallel for reduction(+:llval) num_threads(nThreads_c) \
                         private(snp, snp_der)
for(ind = 0; ind < nInd_c; ind++){
    // most variables defined within the scope of the parallel loop, so
    // don't need a long private declaration
    int s1, s2;
    double phi[4][nSnps_c], phi_prev[4][nSnps_c];
    double alphaTilde[4], alphaDot[4], sum, sum_der, w_new, w_prev;
    double delta;

    ...

    // add contributions to the score vector
    for(snp_der = 0; snp_der < snp-1; snp_der++){
        sum_der = 0;
        for(s2 = 0; s2 < 4; s2++)
            sum_der = sum_der + phi[s2][snp_der]/w_prev;
        // make this an "atomic" operation to avoid multiple processes
        // writing to the same element of "score_c" at the same time
        #pragma omp atomic
        score_c[snp_der] += sum_der;
    }

    ...
```

In some cases `atomic` operations may not be most efficient, since they can
cause threads to have to wait, instead of running calculations. We tried an
alternative implementation that allocates thread local storage but it turned
out that the `atomic` version was slightly faster in this case. See
[here](https://github.com/chrisdjscott/GUSMap/tree/openmp-thread-local/benchmarks/compare_atomic)
for details.
