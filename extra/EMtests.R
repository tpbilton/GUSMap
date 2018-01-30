FS1 <- simFS(rVec_f = 0.01,rVec_m = 0.01,epsilon = 0.001,config = c(2,1,1,4,2,4,1,1,4,1,2,4), nInd = 2,meanDepth = 6)
FS2 <- simFS(rVec_f = 0.01,rVec_m = 0.01,epsilon = 0.001,config = c(2,1,1,4,2,4,1,1,4,1,2,4), nInd = 3,meanDepth = 6)

FS1$depth_Ref <- matrix(as.numeric(FS1$depth_Ref), nrow=2)
FS2$depth_Ref <- matrix(as.numeric(FS2$depth_Ref), nrow=3)

OPGP <- list(FS1$OPGP, FS2$OPGP)
depth_Ref <- list(FS1$depth_Ref, FS2$depth_Ref)
depth_Alt <- list(FS1$depth_Alt, FS2$depth_Alt)
init_r <- matrix(0.01,nrow=2,ncol=11)
nInd = as.integer(c(2,3))
nSnps = as.integer(12)
noFam = as.integer(2)
sexSpec = FALSE
epsilon = 0.001

dyn.load("em.so")


OPGPmat = do.call(what = "rbind",OPGP)
depth_Ref_mat = do.call(what = "rbind",depth_Ref)
depth_Alt_mat = do.call(what = "rbind",depth_Alt)

EMout <- .Call("EM_HMM", init_r, epsilon, depth_Ref_mat, depth_Alt_mat, OPGPmat,
               noFam, nInd, nSnps, sexSpec)