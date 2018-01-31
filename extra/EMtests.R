FS1 <- simFS(rVec_f = 0.1,rVec_m = 0.1,epsilon = 0.000,config = c(2,1,1,4,2,4,1,1,4,1,2,4), nInd = 50,meanDepth = 10, seed2=96843)
FS2 <- simFS(rVec_f = 0.1,rVec_m = 0.1,epsilon = 0.000,config = c(2,1,1,4,2,4,1,1,4,1,2,4), nInd = 30,meanDepth = 10,seed1 = 94732, seed2 = 1903423)

FS1$depth_Ref <- matrix(as.numeric(FS1$depth_Ref), nrow=50)
FS2$depth_Ref <- matrix(as.numeric(FS2$depth_Ref), nrow=30)

OPGP <- list(FS1$OPGP, FS2$OPGP)
depth_Ref <- list(FS1$depth_Ref, FS2$depth_Ref)
depth_Alt <- list(FS1$depth_Alt, FS2$depth_Alt)
init_r <- matrix(0.1,nrow=2,ncol=11)
nInd = as.integer(c(50,30))
nSnps = as.integer(12)
noFam = as.integer(2)
sexSpec = FALSE
epsilon = 0.001

para = c(1000,1e-10)

OPGPmat = do.call(what = "rbind",OPGP)
depth_Ref_mat = do.call(what = "rbind",depth_Ref)
depth_Alt_mat = do.call(what = "rbind",depth_Alt)

dyn.load("em.so")
EMout <- .Call("EM_HMM", init_r, epsilon, depth_Ref_mat, depth_Alt_mat, OPGPmat, noFam, nInd, nSnps, F, F, para)
OPTIMout <- rf_est_FS(init_r=0.05,epsilon = NULL,depth_Ref=depth_Ref,depth_Alt=depth_Alt,OPGP=OPGP, noFam=2,sexSpec=F)
round(EMout[[1]],6)
round(OPTIMout[[1]],6)

EMout <- .Call("EM_HMM", init_r, epsilon, matrix(as.numeric(FS1$depth_Ref),nrow=1), FS1$depth_Alt, FS1$OPGP, as.integer(1), as.integer(50), nSnps, F, F, para)

OPTIMout <-rf_est_FS(init_r= 0.1,epsilon=NULL,depth_Ref=list(FS1$depth_Ref),depth_Alt=list(FS1$depth_Alt), OPGP=list(FS1$OPGP))
