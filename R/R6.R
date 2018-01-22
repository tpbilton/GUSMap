

### R6 class for data aligned to reference assembly
RA <- R6Class("RA",
  public = list(
    initialize = function(List){
      private$genon     <- List$genon
      private$depth_Ref <- List$depth_Ref
      private$depth_Alt <- List$depth_Alt
      private$chrom     <- List$chrom
      private$pos       <- List$pos
      private$SNP_Names <- List$SNP_Names
      private$indID     <- List$indID
      private$nSnps     <- List$nSnps
      private$nInd      <- List$nInd
      private$gform     <- List$gform
    },
    #### Diagonostic functions ####
    ## Ratio of alleles for heterozygous genotype calls (observed vs expected)
    HeteroPlot = function(model=c("random", "beta-binom"), alpha=NULL, filename="HeteroPlot", cex=1){
      HeteroPlot(private$depth_Ref, private$depth_Alt, model=model, alpha=alpha, filename=filename, cex=cex)
    }
    ###############################
  ),
  private = list(
    genon = NULL,
    depth_Ref = NULL,
    depth_Alt = NULL,
    chrom = NULL,
    pos = NULL,
    SNP_Names = NULL,
    indID = NULL,
    nSnps = NULL,
    nInd = NULL,
    gform = NULL
  )
)

### R6 class for creating a data format for a full-sib families
FS <- R6Class("FS",
  inherit = RA,
  public = list(
    ## variables
    LG_mat    = NULL,
    LG_pat    = NULL,
    seg_group = NULL,
    filter    = NULL,
    seg_num   = NULL,
    OPGP      = NULL,
    ## initialize function
    initialize = function(R6obj){
      self$LG_mat       <- list()
      self$LG_pat       <- list()
      self$seg_group    <- list()
      self$filter       <- list()
      self$seg_num      <- list()
      self$OPGP         <- list()
      private$genon     <- R6obj$.__enclos_env__$private$genon
      private$depth_Ref <- R6obj$.__enclos_env__$private$depth_Ref
      private$depth_Alt <- R6obj$.__enclos_env__$private$depth_Alt
      private$chrom     <- R6obj$.__enclos_env__$private$chrom
      private$pos       <- R6obj$.__enclos_env__$private$pos
      private$SNP_Names <- R6obj$.__enclos_env__$private$SNP_Names
      private$indID     <- R6obj$.__enclos_env__$private$indID
      private$nSnps     <- R6obj$.__enclos_env__$private$nSnps
      private$nInd      <- R6obj$.__enclos_env__$private$nInd
      private$gform     <- R6obj$.__enclos_env__$private$gform
    },
    #############################################################
    ## Function for computing the 2-point rf estimates
    rf_2pt = function(nClust=4){
      ## do some checks
      if(!is.numeric(nClust) || length(nClust)!=1 || nClust < 0 || is.infinite(nClust) )
        stop("Number clusters specifies for the parallelization needs to be a positive finite numeric number")
      
      ## If there is only one family
      if(private$noFam == 1)
        mat <- rf_2pt_single(private$depth_Ref[[1]], private$depth_Alt[[1]],
                             private$config[[1]], private$config_infer[[1]],
                             private$group, private$group_infer,
                             nClust)
      else
        mat <- rf_2pt_multi(private$depth_Ref, private$depth_Alt,
                            private$config,private$group, nClust, private$noFam)
      ## Save the results to the object
      private$rf <- mat$rf
      private$LOD <- mat$LOD
      return(invisible())
    },
    ## Function for creating linkage groups
    createLG = function(parent, LODthres=10, nComp=10){
      ## Do some checks
      if(!is.character(parent) || length(parent) != 1 || !(parent %in% c("maternal","paternal")))
        stop("parent argument is not a string of of length one or is incorrect:
             Please select one of the following:
             maternal: Only MI SNPs
             paternal: Only PI SNPs")
      if(!is.numeric(LODthres) || LODthres < 0 || !is.finite(LODthres) || length(LODthres) != 1)
        stop("The LOD threshold (argument 2) needs to be a finite numeric number.")
      if(!is.integer(nComp) || nComp < 0 || !is.finite(nComp) || length(nComp) != 1)
        stop("The number of comparsion points (argument 3) needs to be a finite numeric number.")
      
      ## Create the groups
      if(parent == "maternal")
        self$LG_mat <- createLG(private$group, private$LOD, "maternal", LODthres, nComp)
      else if(parent == "paternal")
        self$LG_pat <- createLG(private$group, private$LOD, "paternal", LODthres, nComp)
    },
    ## Function for plotting linkage groups
    plotLG = function(mat=c("rf","LOD"), LG, filename=NULL, names=NULL, chrS=2, lmai=2, chrom=T){
      
      if(mat == "rf")
        plotLG(mat=private$rf, LG=LG, filename=filename, names=names, chrS=chrS, lmai=lmai, chrom=chrom)
      else if(mat == "LOD")
        plotLG(mat=private$LOD, LG=LG, filename=filename, names=names, chrS=chrS, lmai=lmai, chrom=chrom)
      else
        stop("Matrix to be plotted not found.")
    },
    ## Ratio of alleles for heterozygous genotype calls (observed vs expected)
    # Slightly different than the one in RA class
    HeteroPlot = function(model=c("random"), alpha=NULL, filename="HeteroPlot", cex=1){
      depth_Ref <- depth_Alt <- NULL
      for(fam in 1:private$noFam){
        depth_Ref <- rbind(depth_Ref,private$depth_Ref[[fam]])
        depth_Alt <- rbind(depth_Alt,private$depth_Alt[[fam]])
      }
      HeteroPlot(depth_Ref, depth_Alt, model=model, alpha=alpha, filename=filename, cex=cex)
    }
    
    ##############################################################
  ),
  private = list(
    config       = NULL,
    config_infer = NULL,
    group        = NULL,
    group_infer  = NULL,
    noFam        = NULL,
    rf           = NULL,
    LOD          = NULL,
    ## Function for updating the private variables
    updatePrivate = function(List){
      if(!is.null(List$genon))
        private$genon        = List$genon
      if(!is.null(List$depth_Ref))
        private$depth_Ref    = List$depth_Ref
      if(!is.null(List$depth_Alt))
        private$depth_Alt    = List$depth_Alt
      if(!is.null(List$chrom))
        private$chrom        = List$chrom
      if(!is.null(List$pos))
        private$pos          = List$pos
      if(!is.null(List$SNP_Names))
        private$SNP_Names    = List$SNP_Names
      if(!is.null(List$indID))
        private$indID        = List$indID
      if(!is.null(List$nSnps))
        private$nSnps        = List$nSnps
      if(!is.null(List$nInd))
        private$nInd         = List$nInd
      if(!is.null(List$config))
        private$config       = List$config
      if(!is.null(List$config_infer))
        private$config_infer = List$config_infer
      if(!is.null(List$group))
        private$group        = List$group
      if(!is.null(List$group_infer))
        private$group_infer  = List$group_infer
      if(!is.null(List$noFam))
        private$noFam        = List$noFam
      if(!is.null(List$rf))
        private$rf           = List$rf
      if(!is.null(List$LOD))
        private$LOD          = List$LOD
    }
  )
)


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
#           rf.est1 <- rf_est_FS(depth_Ref=list(obj$depth_Ref[,ind]),depth_Alt=list(obj$depth_Alt[,ind]),
#                                OPGP=list(c(5,5) + 2*(config[ind]==3)), epsilon=NULL)
#           rf.est2 <- rf_est_FS(depth_Ref=list(obj$depth_Ref[,ind]),depth_Alt=list(obj$depth_Alt[,ind]),
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
#           rf.est1 <- rf_est_FS(depth_Ref=list(obj$depth_Ref[,ind]),depth_Alt=list(obj$depth_Alt[,ind]),
#                                OPGP=list(c(9,9) + 2*(config[ind]==5)), epsilon=NULL)
#           rf.est2 <- rf_est_FS(depth_Ref=list(obj$depth_Ref[,ind]),depth_Alt=list(obj$depth_Alt[,ind]),
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
#           temp1 <- rf_est_FS(0.1,depth_Ref=list(obj$depth_Ref[,ind]),depth_Alt=list(obj$depth_Alt[,ind]), OPGP=list(c(1,1)), epsilon=NULL)
#           temp2 <- rf_est_FS(0.4,depth_Ref=list(obj$depth_Ref[,ind]),depth_Alt=list(obj$depth_Alt[,ind]), OPGP=list(c(1,1)), epsilon=NULL)
#           rf.est1 <- switch(which.min(c(temp1$loglik,temp2$loglik)),temp1,temp2) 
#           temp1 <- rf_est_FS(0.1,depth_Ref=list(obj$depth_Ref[,ind]),depth_Alt=list(obj$depth_Alt[,ind]), OPGP=list(c(1,2)), epsilon=NULL)
#           temp2 <- rf_est_FS(0.4,depth_Ref=list(obj$depth_Ref[,ind]),depth_Alt=list(obj$depth_Alt[,ind]), OPGP=list(c(1,2)), epsilon=NULL)
#           rf.est2 <- switch(which.min(c(temp1$loglik,temp2$loglik)),temp1,temp2) 
#           temp1 <- rf_est_FS(0.1,depth_Ref=list(obj$depth_Ref[,ind]),depth_Alt=list(obj$depth_Alt[,ind]), OPGP=list(c(1,4)), epsilon=NULL)
#           temp2 <- rf_est_FS(0.4,depth_Ref=list(obj$depth_Ref[,ind]),depth_Alt=list(obj$depth_Alt[,ind]), OPGP=list(c(1,4)), epsilon=NULL)
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
#           rf.est1 <- rf_est_FS(depth_Ref=list(obj$depth_Ref[,ind]),depth_Alt=list(obj$depth_Alt[,ind]),
#                                OPGP=list(c(5,1) + 2*c(config[ind[1]]==3,0)), epsilon=NULL )
#           rf.est2 <- rf_est_FS(depth_Ref=list(obj$depth_Ref[,ind]),depth_Alt=list(obj$depth_Alt[,ind]),
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
#           rf.est1 <- rf_est_FS(depth_Ref=list(obj$depth_Ref[,ind]),depth_Alt=list(obj$depth_Alt[,ind]),
#                                OPGP=list(c(9,1) + 2*c(config[ind[1]]==5,0)), epsilon=NULL)
#           rf.est2 <- rf_est_FS(depth_Ref=list(obj$depth_Ref[,ind]),depth_Alt=list(obj$depth_Alt[,ind]),
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
#           rf.est1 <- rf_est_FS(depth_Ref=list(obj$depth_Ref[,ind]),depth_Alt=list(obj$depth_Alt[,ind]),
#                                OPGP=list(c(9,9) + 2*(config[ind] %in% c(3,5))), epsilon=NULL)
#           rf.est2 <- rf_est_FS(depth_Ref=list(obj$depth_Ref[,ind]),depth_Alt=list(obj$depth_Alt[,ind]),
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
#     depth_Ref = NULL,
#     depth_Alt = NULL,
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
#       private$depth_Ref    <- R6obj$.__enclos_env__$private$depth_Ref
#       private$depth_Alt    <- R6obj$.__enclos_env__$private$depth_Alt
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
