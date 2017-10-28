##### File for implementing the SPIN algorithm for ordering a set of loci in the context of LD
##### The following adjustments are made to the orriginal implementation of SPIN in the R package seriation
##### a) Approximation step of solving the LASP is changed to be solved using the hugerian algorithm
##### b)
####
### Author: Timothy Bilton
### Date: 31/01/17


############# Package required ###################
library(seriation) # Contains an implementation of the SPIN algorithm.
library(clue)      # Contains function solve_LSAP to solve the LASP step


## Define new method of seriation in the seriation package.
## Specifically designed for LD analysis
## Note: Based on the seriate:dist method
seriate.distLD <- function (x, method = "SPIN_NH_LD", control = NULL, ...) 
{
  if (!all(x >= 0)) 
    stop("Negative distances not supported!")
  control <- c(control, list(...))
  if (any(is.na(x))) 
    stop("NAs not allowed in x!")
  if (any(x < 0)) 
    stop("No negative values allowed in x!")
  if (!is.character(method) || (length(method) != 1L)) 
    stop("Argument 'method' must be a character string.")
  method <- get_seriation_method("distLD", method)
  if (!is.null(control$verbose) && control$verbose) 
    cat(method$name, ": ", method$description, "\n", sep = "")
  ## Changes from the seriate:dist method
  orderL <- method$fun(x, control = control) # is a list of two objects
  # return the element in the list orderL
  return(new_order=list(ser_permutation(ser_permutation_vector(orderL[[1]], method = method$name)),
                        Best_energy=orderL[[2]]))
}


#### Define function for the new seriation method which has the new implementation of the SPIN neighbourhood algorithm for LD analysis.
SPIN_NH_LD <- function(x,control){
  param <- seriation:::.get_parameters(control, list(sigma = seq(20, 1, length.out = 10), step = 100, W_function = NULL, verbose = FALSE))
  ## If weight matrix is not pre-defined, create a function to compute the weight matrix proposed by Tsafrir et al (2005) 
  W_function <- if (is.null(param$W_function)) 
    create_W
  else param$W_function
  ## extract parameter values
  sigma <- param$sigma
  step <- param$step
  verbose <- param$verbose
  D <- as.matrix(x)
  n <- nrow(D)
  ## Generate weight matrix
  W <- W_orig <- W_function(n, sigma[1], verbose)

  ## Run the neighbourhood algorithm
  # Step 1:
  energy_prev <- 100
  energy_new <- 1000
  i = 0
  # iterate over the remaining steps
  while((abs(energy_prev - energy_new) > 1e-40) & (i<step)) {
    i = i + 1
    energy_prev <- energy_new
    if (verbose) 
      cat("Iteration", i, "... ")
    # Step 2:
    M <- D %*% W
    # Step 3 (Different from the seriate method SPIN_NH) 
    P <- permutation_vector2matrix(order(as.vector(solve_LSAP(M))))
    # Step 4: (Different from the seriate method SPIN_NH) 
    energy_new <- sum(diag(P %*% M))
    if (verbose) 
      cat("new energy: ", energy_new, "\n")
    W <- crossprod(P, W_orig)
  }
  ## Once algorithm has converged, extract the new order and return the new energy
  o <- permutation_matrix2vector(P)
  names(o) <- names(x)[o]
  return(list(o,energy_new))
}
## Set the new seriation method in seriation package
set_seriation_method(kind="distLD", name="SPIN_NH_LD", definition=SPIN_NH_LD,
                     description="SPIN (Neignbourhood algorithm for LD)")

## Function for running the Neighbourhood algorithm using a different starting 
orderLoci <- function(LDmat,nIter=10,sigma=10){
  nSnps <- nrow(LDmat)
  distMat <- as.dist(abs(1-LDmat))
  
  ## try initial order
  out <- seriate.distLD(distMat,method="SPIN_NH_LD",control=list(sigma=sigma,verbose=T))
  best_energy <- out[[2]]
  best_order <- get_order(out[[1]])
  
  ## Try randomizing starting order
  for(i in 1:nIter){
    initOrd <- sample(1:nSnps,size=nSnps)
    out <- seriate.distLD(distMat,method="SPIN_NH_LD",control=list(sigma=sigma,verbose=T))
    if(best_energy > out[[2]]){
      best_energy <- out[2]
      best_order <- get_order(out[[1]])
    }
  }
  return(list(best_order,best_energy))
}


create_W <- function (n, sigma, verbose = FALSE)
{
  w <- function(i, j, n, sigma) exp(-1 * (i - j)^2/n/sigma)
  W <- outer(1:n, 1:n, FUN = w, n = n, sigma = sigma)
  for (i in 1:1000) {
    W <- sweep(W, MARGIN = 1, STATS = rowSums(W), "/")
    W <- sweep(W, MARGIN = 2, STATS = colSums(W), "/")
    if (round(rowSums(W), 20) == 1 && round(colSums(W), 20) ==
        1)
      break
  }
  if (verbose)
    cat("It took", i, "iterations to make W doubly stochastic!\n")
  if (i > 999)
    warning("Weight matrix did not converge to doubly stochastic in 1000 itermation!")
  W
}

