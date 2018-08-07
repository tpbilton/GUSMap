# Optimise call to rfEst for smaller case

Run code and `save` object to file just before the call, so we can just focus
on improving that section of code.

## Adding OpenMP

Loops in the following C functions were parallelised with OpenMP:

* `EM_HMM`
* `EM_HMM_UP`
* `ll_fs_scaled_err_c`
* `score_fs_scaled_err_c`

After making these changes the code was benchmarked (note: timings are just
for the call to `MK_fs$rf_est`):

![benchmark results](gusmap-smaller-openmp-compare.png)

## OpenMP notes

OpenMP is an API for shared memory multiprocessing (i.e. within a node)

### Setting the number of threads to use

* using an environment variable: `export OMP_NUM_THREADS=8`
* in Slurm:
  ```
  #!/bin/bash
  #SBATCH --cpus-per-task=8
  # other SBATCH directives

  export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
  srun <my_program>
  ```

### Simple parallel loop

For a very simple loop:

```c
for (int i = 0; i < size; i++) {
    my_array[i] = i * i;
}
```

it can be as easy as adding one line:

```c
#pragma omp parallel for
for (int i = 0; i < size; i++) {
    my_array[i] = i;
}
```

The `#pragma omp parallel for` tells the compiler that the following loop
should be run in parallel. Note, if the compiler does not support OpenMP it
ignores this line and the loop will run in serial (there won't be an error).

Some loops cannot be parallelised, for example if each iteration uses a value
computed at a previous iteration, e.g.:

```c
my_array[0] = 0;
for (int i = 1; i < size; i++) {
    my_array[i] = my_array[i - 1] + i;
}
```

If you try to parallelise this loop then you will probably get a completely
wrong result (it may not fail though).

### Private vs shared variables

Variables that are inside a parallel loop are either `shared` or `private` and
the default is for them to be `shared`. Consider the following loop:

```c
int my_array[size];
int my_value;
for (int i = 0; i < size; i++) {
    my_value = i * i;
    my_array[i] = my_value;
}
```

By default, both `my_array` and `my_value` are shared between parallel
processes. This is fine for `my_array`; we want it to be shared so that
multiple processes can populate its values at the same time (note, however,
that at each iteration a different element will be populated, no two processes
will ever try to populate the same element of the array).

However, `my_value` is also shared. This is bad, because multiple processes
will be writing to that value at the same time and the result will probably be
wrong. There are (at least) two ways around this. First, one could declare the
variable as `private`:

```c
int my_array[size];
int my_value;
#pragma omp parallel for private(my_value)
for (int i = 0; i < size; i++) {
    my_value = i * i;
    my_array[i] = my_value;
}
```

Secondly, one could define the variables that are private to a loop only
within the scope of that loop, making them implicitly private.

```c
int my_array[size];
#pragma omp parallel for
for (int i = 0; i < size; i++) {
    int my_value;
    my_value = i * i;
    my_array[i] = my_value;
}
```

If you want to be really safe, or you have lots of variables and it is
difficult to work out which should be private, you can make it so that all
variables must be explicitly declared as private or shared (excluding those
variables that are defined only within the scope of the loop):

```c
int my_array[size];
int my_value;
#pragma omp parallel for default(none) private(my_value) shared(my_array)
for (int i = 0; i < size; i++) {
    my_value = i * i;
    my_array[i] = my_value;
}
```

### Reductions

OpenMP support doing parallel reductions; a common operation is summing some
values in a loop. This is achieved with the `reduction` clause, which supports
a number of operators.

```c
int my_sum = 0;
#pragma omp parallel for reduction(+:my_sum)
for (int i = 0; i < size; i++) {
    my_sum += my_array[i];
}
```

Here, each thread will have its own copy of `my_sum`, initialised to zero. At
the end of the loop the values from each thread will be reduced to one value
using the operator specified in the `reduction` clause (`+` here).

Note, if you did not use the `reduction` clause, then multiple processes could
write to `my_sum` at the same time and cause the result to be wrong.

### Atomic operations

When an operation is marked as `atomic` it is guaranteed that only one process
will perform that operation at a given time. Therefore it is best to try to
avoid `atomic` operations if possible, since they can result in a process
spending time waiting for other processes to complete the operation before it
can proceed. However, sometimes `atomic` operations are required.

```c
#pragma omp parallel for
for (int i = 0; i < size1; i++) {
    for (int j = 0; j < size2; j++) {
        int value = my_array[i] * j;
        #pragma omp atomic
        score[j] += value;
    }
}
```

Here, multiple processes, which are running the `i` loop in parallel, could
try to access the same element of the `score` array at the same time. The
`atomic` clause ensures this does not happen.
