# load dev version
# assume GUSMAP_DIR is set
gusmap_src_dir <- Sys.getenv('GUSMAP_DIR')
cat("Loading GUSMap from:", gusmap_src_dir, "\n")
library(devtools)
load_all(gusmap_src_dir)
library(tictoc)

# set number of thread from environment variable OMP_NUM_THREADS

num_threads <- max(strtoi(Sys.getenv('OMP_NUM_THREADS')), 1)
cat("Running with num_threads =", num_threads, "\n")

# Code for testing GUSMap

tic("RTIME total")

library(GUSMap)
library(onemap)

load(file = "MK_fs.rdata")
save(MK_fs, file = "MK_fs.rdata")

## Compute the rf estimates
## using optim
tic("RTIME rf_est (optim)")
MK_fs$rf_est(method="optim")
tmp <- sum(haldane(MK_fs$para$rf_p[[1]])) ## should be 76.09481
cat(tmp, 76.09481, "\n")
tmp <- MK_fs$para$ep[[1]]                 ## should be 0.003167973
cat(tmp, 0.003167973, "\n")
toc()

toc()
