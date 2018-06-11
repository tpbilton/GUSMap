


# 
# FS <- R6Class("FS",
#   public = list(
#     indx = NULL,
#     initialize = function(indx = NA) {
#       self$indx <- indx
#     },
#     ## Now include a function to compute the 2-point estimates
#     rf_2pt <- function(nClust, inferSNPs = TRUE){
#       
#       if(length(c(unlist(obj$group), unlist(obj$group_infer))) == 0)
#         stop("There are no SNPs in the data set.")
#       
#       ## Cacluate the number of SNPs
#       if(length(obj$group_infer$BI) > 0 & inferSNPs){
#         indx_BI <- c(private$group$BI, private$group_infer$BI)
#         nSnps_BI <- length(indx_BI)
#       }
#       else{
#         indx_BI <- c(obj$group$BI)
#         nSnps_BI <- length(indx_BI)
#       }
#       if(length(obj$group_infer$SI) > 0 & inferSNPs){
#         indx_MI <- c(obj$group$MI, obj$group_infer$SI)
#         nSnps_MI <- length(indx_MI)
#       }
#       else{
#         indx_MI <- c(obj$group$MI)
#         nSnps_MI <- length(indx_MI)
#       }
#       indx_PI <- c(obj$group$PI)
#       nSnps_PI = length(indx_PI)
#       
#       ## set up the config vector
#       config <-  obj$config
#       config[which(!is.na(obj$config_infer))] <- obj$config_infer[which(!is.na(obj$config_infer))]
#       if(any(is.na(config)))
#         stop("There are some missing segregation types in the data.")
#       
#       ## Set up the Clusters
#       cl <- makeCluster(nClust)
#       registerDoSNOW(cl)
#       
#       cat("\nComputing 2-point recombination fraction estimates ...\n")
#       cat("Paternal SNPs\n")
#       ## Paternal informative SNPs
#       rf.PI <- foreach(snp1 = iter(seq(length.out=nSnps_PI)), .combine=comb) %dopar% {
#         rf <- replicate(2,numeric(nSnps_PI),simplify=F)
#         for(snp2 in seq_len(snp1-1)){
#           ind = indx_PI[c(snp1,snp2)]
#           rf.est1 <- rf_est_FS(ref=list(obj$ref[,ind]),alt=list(obj$alt[,ind]),
#                                OPGP=list(c(5,5) + 2*(config[ind]==3)), epsilon=NULL)
#           rf.est2 <- rf_est_FS(ref=list(obj$ref[,ind]),alt=list(obj$alt[,ind]),
#                                OPGP=list(c(5,6) + 2*(config[ind]==3)), epsilon=NULL)
#           rf.ind <- switch(which.min(c(rf.est1$loglik,rf.est2$loglik)), rf.est1, rf.est2)
#           rf[[1]][snp2] <- rf.ind$rf
#           rf[[2]][snp2] <- rf.ind$LOD
#         }
#         return(rf)
#       }
#       for(i in 1:2){
#         rf.PI[[i]][upper.tri(rf.PI[[i]])] <- t(rf.PI[[i]])[upper.tri(rf.PI[[i]])]
#       }
#       
#       cat("Maternal SNPs\n")
#       ## Maternal informative SNPs
#       rf.MI <- foreach(snp1 = iter(seq(length.out=nSnps_MI)), .combine=comb) %dopar% {
#         rf <- replicate(2,numeric(nSnps_MI),simplify=F)
#         for(snp2 in seq_len(snp1-1)){
#           ind = indx_MI[c(snp1,snp2)]
#           rf.est1 <- rf_est_FS(ref=list(obj$ref[,ind]),alt=list(obj$alt[,ind]),
#                                OPGP=list(c(9,9) + 2*(config[ind]==5)), epsilon=NULL)
#           rf.est2 <- rf_est_FS(ref=list(obj$ref[,ind]),alt=list(obj$alt[,ind]),
#                                OPGP=list(c(9,10) + 2*(config[ind]==5)), epsilon=NULL)
#           rf.ind <- switch(which.min(c(rf.est1$loglik,rf.est2$loglik)), rf.est1, rf.est2)
#           rf[[1]][snp2] <- rf.ind$rf
#           rf[[2]][snp2] <- rf.ind$LOD
#         }
#         return(rf)
#       }  
#       for(i in 1:2){
#         rf.MI[[i]][upper.tri(rf.MI[[i]])] <- t(rf.MI[[i]])[upper.tri(rf.MI[[i]])]
#       }
#       
#       cat("Both informative SNPs\n")
#       ### Both Informative
#       rf.BI <- foreach(snp1 = iter(seq(length.out=nSnps_BI)), .combine=comb) %dopar% {
#         rf <- replicate(2,numeric(nSnps_BI),simplify=F)
#         for(snp2 in seq_len(snp1-1)){
#           ind = indx_BI[c(snp1,snp2)]
#           temp1 <- rf_est_FS(0.1,ref=list(obj$ref[,ind]),alt=list(obj$alt[,ind]), OPGP=list(c(1,1)), epsilon=NULL)
#           temp2 <- rf_est_FS(0.4,ref=list(obj$ref[,ind]),alt=list(obj$alt[,ind]), OPGP=list(c(1,1)), epsilon=NULL)
#           rf.est1 <- switch(which.min(c(temp1$loglik,temp2$loglik)),temp1,temp2) 
#           temp1 <- rf_est_FS(0.1,ref=list(obj$ref[,ind]),alt=list(obj$alt[,ind]), OPGP=list(c(1,2)), epsilon=NULL)
#           temp2 <- rf_est_FS(0.4,ref=list(obj$ref[,ind]),alt=list(obj$alt[,ind]), OPGP=list(c(1,2)), epsilon=NULL)
#           rf.est2 <- switch(which.min(c(temp1$loglik,temp2$loglik)),temp1,temp2) 
#           temp1 <- rf_est_FS(0.1,ref=list(obj$ref[,ind]),alt=list(obj$alt[,ind]), OPGP=list(c(1,4)), epsilon=NULL)
#           temp2 <- rf_est_FS(0.4,ref=list(obj$ref[,ind]),alt=list(obj$alt[,ind]), OPGP=list(c(1,4)), epsilon=NULL)
#           rf.est4 <- switch(which.min(c(temp1$loglik,temp2$loglik)),temp1,temp2) 
#           rf.ind <- switch(which.min(c(rf.est1$loglik, rf.est2$loglik, rf.est4$loglik) ),
#                            rf.est1,rf.est2,rf.est4)
#           rf[[1]][snp2] <- rf.ind$rf
#           rf[[2]][snp2] <- rf.ind$LOD
#         }
#         return(rf)
#       }
#       for(i in 1:2){
#         rf.BI[[i]][upper.tri(rf.BI[[i]])] <- t(rf.BI[[i]])[upper.tri(rf.BI[[i]])]
#       }
#       
#       cat("Paternal information vs Both informative\n")
#       ## Paternal and Informative SNPs
#       rf.PI.BI <- foreach(snp.ps = iter(seq(length.out=nSnps_PI)), .combine=comb) %dopar% {
#         rf <- replicate(2,numeric(nSnps_BI),simplify=F)
#         for(snp.bi in 1:nSnps_BI){
#           ind <- c(indx_PI[snp.ps],indx_BI[snp.bi])
#           rf.est1 <- rf_est_FS(ref=list(obj$ref[,ind]),alt=list(obj$alt[,ind]),
#                                OPGP=list(c(5,1) + 2*c(config[ind[1]]==3,0)), epsilon=NULL )
#           rf.est2 <- rf_est_FS(ref=list(obj$ref[,ind]),alt=list(obj$alt[,ind]),
#                                OPGP=list(c(5,2) + 2*c(config[ind[1]]==3,0)), epsilon=NULL )
#           rf.ind <- switch(which.min(c(rf.est1$loglik,rf.est2$loglik)), rf.est1, rf.est2)
#           rf[[1]][snp.bi] <- rf.ind$rf
#           rf[[2]][snp.bi] <- rf.ind$LOD
#         }
#         return(rf)
#       }
#       
#       cat("Maternal information vs Both informative\n")
#       ## Maternal and Informative SNPs
#       rf.MI.BI <- foreach(snp.ms = iter(seq(length.out=nSnps_MI)), .combine=comb) %dopar% {
#         rf <- replicate(2,numeric(nSnps_BI),simplify=F)
#         for(snp.bi in 1:nSnps_BI){
#           ind <- c(indx_MI[snp.ms],indx_BI[snp.bi])
#           rf.est1 <- rf_est_FS(ref=list(obj$ref[,ind]),alt=list(obj$alt[,ind]),
#                                OPGP=list(c(9,1) + 2*c(config[ind[1]]==5,0)), epsilon=NULL)
#           rf.est2 <- rf_est_FS(ref=list(obj$ref[,ind]),alt=list(obj$alt[,ind]),
#                                OPGP=list(c(9,3) + 2*c(config[ind[1]]==5,0)), epsilon=NULL)
#           rf.ind <- switch(which.min(c(rf.est1$loglik,rf.est2$loglik)), rf.est1, rf.est2)
#           rf[[1]][snp.bi] <- rf.ind$rf
#           rf[[2]][snp.bi] <- rf.ind$LOD
#         }
#         return(rf)
#       }
#       
#       ## For the non-informative computations
#       ## Really done so that we can check that there is no miss identification of the group
#       cat("Maternal information vs Paternal informative\n")
#       rf.MI.PI <- foreach(snp.ms = iter(seq(length.out=nSnps_MI)), .combine=comb) %dopar% {
#         rf <- replicate(2,numeric(nSnps_PI),simplify=F)
#         for(snp.pi in 1:nSnps_PI){
#           ind <- c(indx_MI[snp.ms],indx_PI[snp.pi])
#           rf.est1 <- rf_est_FS(ref=list(obj$ref[,ind]),alt=list(obj$alt[,ind]),
#                                OPGP=list(c(9,9) + 2*(config[ind] %in% c(3,5))), epsilon=NULL)
#           rf.est2 <- rf_est_FS(ref=list(obj$ref[,ind]),alt=list(obj$alt[,ind]),
#                                OPGP=list(c(9,10) + 2*(config[ind] %in% c(3,5))), epsilon=NULL)
#           rf.ind <- switch(which.min(c(rf.est1$loglik,rf.est2$loglik)), rf.est1, rf.est2)
#           rf[[1]][snp.pi] <- rf.ind$rf
#           rf[[2]][snp.pi] <- rf.ind$LOD
#         }
#         return(rf)
#       }
#       stopCluster(cl) 
#       
#       ## Build the rf and LOD matrices
#       origOrder <- order(c(indx_BI,indx_PI,indx_MI))
#       rf.mat <- rbind(cbind(rf.BI[[1]],t(rf.PI.BI[[1]]),t(rf.MI.BI[[1]])),
#                       cbind(rf.PI.BI[[1]], rf.PI[[1]], t(rf.MI.PI[[1]])),
#                       cbind(rf.MI.BI[[1]], rf.MI.PI[[1]], rf.MI[[1]]))[origOrder,origOrder]
#       LOD.mat <- rbind(cbind(rf.BI[[2]],t(rf.PI.BI[[2]]),t(rf.MI.BI[[2]])),
#                        cbind(rf.PI.BI[[2]], rf.PI[[2]], t(rf.MI.PI[[2]])),
#                        cbind(rf.MI.BI[[2]], rf.MI.PI[[2]], rf.MI[[2]]))[origOrder,origOrder]
#       obj$rf <- rf.mat
#       obj$LOD <- LOD.mat
#       return(obj)
#     }
#   ), 
#   ### Private stuff
#   private = list(
#     genon = NULL,
#     ref = NULL,
#     alt = NULL,
#     chrom = NULL,
#     pos = NULL,
#     indID = NULL,
#     config = NULL,
#     config_infer = NULL,
#     group = NULL,
#     group_infer = NULL,
#     nInd = NULL,
#     nSnps = NULL,
#     noFam = NULL
#   )
# )
# 
# multi <- R6Class("multi",
#   inherit = FS,
#   public = list(
#     LGs_mat = NULL,
#     LGs_pat = NULL,
#     initialize = function(R6obj){
#       LGs_mat = list()
#       LGs_pat = list()
#       private$init_private(R6obj)
#     },
#   ),
#   private = list(
#     init_private = function(R6obj){
#       private$genon        <- R6obj$.__enclos_env__$private$genon
#       private$ref    <- R6obj$.__enclos_env__$private$ref
#       private$alt    <- R6obj$.__enclos_env__$private$alt
#       private$chrom        <- R6obj$.__enclos_env__$private$chrom
#       private$pos          <- R6obj$.__enclos_env__$private$pos
#       private$indID        <- R6obj$.__enclos_env__$private$indID
#       private$config       <- R6obj$.__enclos_env__$private$config
#       private$config_infer <- R6obj$.__enclos_env__$private$config_infer
#       private$group        <- R6obj$.__enclos_env__$private$group
#       private$group_infer  <- R6obj$.__enclos_env__$private$group_infer
#       private$nInd         <- R6obj$.__enclos_env__$private$nInd
#       private$nSnps        <- R6obj$.__enclos_env__$private$nSnps
#       private$noFam        <- R6obj$.__enclos_env__$private$noFam
#     }
#   )
# )
# 
# 
# 
# 
