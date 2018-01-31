FS1 <- simFS(rVec_f = 0.1,rVec_m = 0.1,epsilon = 0.000,config = c(2,1,1,4,2), nInd = 5,meanDepth = 10, seed2=74632)
FS2 <- simFS(rVec_f = 0.1,rVec_m = 0.1,epsilon = 0.000,config = c(2,1,1,4,2), nInd = 5,meanDepth = 10,seed1 = 94732, seed2 = 1109584)

FS1$depth_Ref <- matrix(as.numeric(FS1$depth_Ref), nrow=FS1$nInd)
FS2$depth_Ref <- matrix(as.numeric(FS2$depth_Ref), nrow=FS2$nInd)

OPGP <- list(FS1$OPGP, FS2$OPGP)
depth_Ref <- list(FS1$depth_Ref, FS2$depth_Ref)
depth_Alt <- list(FS1$depth_Alt, FS2$depth_Alt)
init_r <- matrix(0.05,nrow=2,ncol=FS2$nSnps)
nInd = as.integer(c(FS1$nInd,FS2$nInd))
nSnps = as.integer(FS2$nSnps)
noFam = as.integer(2)
sexSpec = FALSE
epsilon = 0

para = c(200,1e-40)

OPGPmat = do.call(what = "rbind",OPGP)
depth_Ref_mat = do.call(what = "rbind",depth_Ref)
depth_Alt_mat = do.call(what = "rbind",depth_Alt)

dyn.unload("em.dll")
dyn.load("em.dll")
EMout <- .Call("EM_HMM", init_r, epsilon, depth_Ref_mat, depth_Alt_mat, OPGPmat, noFam, nInd, nSnps, F, F, para)
OPTIMout <- rf_est_FS(init_r=0.05,epsilon = NULL,depth_Ref=depth_Ref,depth_Alt=depth_Alt,OPGP=OPGP, noFam=2,sexSpec=F)
round(EMout[[1]],6)
round(OPTIMout[[1]],6)

EMout <- .Call("EM_HMM", init_r, epsilon, matrix(as.numeric(FS1$depth_Ref),nrow=FS1$nInd), FS1$depth_Alt, FS1$OPGP, as.integer(1), as.integer(FS1$nInd), nSnps, F, F, para)

OPTIMout <-rf_est_FS(init_r= 0.1,epsilon=NULL,depth_Ref=list(FS1$depth_Ref),depth_Alt=list(FS1$depth_Alt), OPGP=list(FS1$OPGP))
