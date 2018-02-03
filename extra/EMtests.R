FS1 <- simFS(rVec_f = 0.1,rVec_m = 0.1,epsilon = 0.01,config = c(2,1,1,4,2,4,4,1,1), nInd = 100,meanDepth = 10, seed2=74632)
FS2 <- simFS(rVec_f = 0.1,rVec_m = 0.1,epsilon = 0.01,config = c(2,1,1,4,2,4,4,1,1), nInd = 3,meanDepth = 10,seed1 = 94732, seed2 = 1109584)

FS1$depth_Ref <- matrix(as.numeric(FS1$depth_Ref), nrow=FS1$nInd)
FS2$depth_Ref <- matrix(as.numeric(FS2$depth_Ref), nrow=FS2$nInd)

OPGP <- list(FS1$OPGP, FS2$OPGP)
depth_Ref <- list(FS1$depth_Ref, FS2$depth_Ref)
depth_Alt <- list(FS1$depth_Alt, FS2$depth_Alt)
init_r <- matrix(0.1,nrow=2,ncol=FS2$nSnps)
nInd = as.integer(c(FS1$nInd,FS2$nInd))
nSnps = as.integer(FS2$nSnps)
noFam = as.integer(2)
sexSpec = FALSE
epsilon = 0.001
config = list(FS1$config, FS2$config)

ps <- sort(unique(unlist(lapply(OPGP,function(x) which(x %in% 1:8)))))[-1] - 1
ms <- sort(unique(unlist(lapply(OPGP,function(x) which(x %in% c(1:4,9:12))))))[-1] - 1
npar <- c(length(ps),length(ms))
ss_rf <- logical(2*(nSnps-1))
ss_rf[ps] <- TRUE
ss_rf[ms + nSnps-1] <- TRUE


para = c(2000,1e-30)

configmat = do.call(what = "rbind",config)
OPGPmat = do.call(what = "rbind",OPGP)
depth_Ref_mat = do.call(what = "rbind",depth_Ref)
depth_Alt_mat = do.call(what = "rbind",depth_Alt)

dyn.unload("em.so")
dyn.load("em.so")
EMout <- .Call("EM_HMM", init_r, epsilon, depth_Ref_mat, depth_Alt_mat, OPGPmat, noFam, nInd, nSnps, F, T, para, as.integer(ss_rf))

OPGP = list(FS1$OPGP)
ps <- sort(unique(unlist(lapply(OPGP,function(x) which(x %in% 1:8)))))[-1] - 1
ms <- sort(unique(unlist(lapply(OPGP,function(x) which(x %in% c(1:4,9:12))))))[-1] - 1
npar <- c(length(ps),length(ms))
ss_rf <- logical(2*(nSnps-1))
ss_rf[ps] <- TRUE
ss_rf[ms + nSnps-1] <- TRUE
EMout <- .Call("EM_HMM_UP", init_r, epsilon, depth_Ref[[1]], depth_Alt[[1]], config[[1]], as.integer(1), nInd[[1]], nSnps, T, para, as.integer(ss_rf))

infer_OPGP_FS(FS1$depth_Ref,FS1$depth_Alt,FS1$config)

OPTIMout <- GUSMap:::rf_est_FS_UP(epsilon = 0.001,depth_Ref=depth_Ref[[1]],depth_Alt=depth_Alt[[1]],config=config[[1]], reltol=1e-10)


OPTIMout <- rf_est_FS(init_r=0.2,epsilon = 0.001,depth_Ref=depth_Ref,depth_Alt=depth_Alt,OPGP=OPGP, noFam=2,sexSpec=F)
round(EMout[[1]],6)
round(OPTIMout[[1]],6)

EMout <- .Call("EM_HMM", init_r, epsilon, matrix(as.numeric(FS1$depth_Ref),nrow=FS1$nInd), FS1$depth_Alt, FS1$OPGP, as.integer(1), as.integer(FS1$nInd), nSnps, F, F, para)

OPTIMout <-rf_est_FS(init_r= 0.1,epsilon=NULL,depth_Ref=list(FS1$depth_Ref),depth_Alt=list(FS1$depth_Alt), OPGP=list(FS1$OPGP))
