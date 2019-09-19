set.seed(6745)
config <- list(list(sample(c(1,2,4), size=30, replace=TRUE)))
F1data <- GUSMap:::simFS(0.01, config=config, meanDepth=5, nInd=100, epsilon=0.001)

ref <- F1data$extractVar(c("ref"))[[1]]
alt <- F1data$extractVar(c("alt"))[[1]]
config <- F1data$extractVar(c("config"))[[1]]
nSnps <- F1data$extractVar(c("nSnps"))[[1]]
nInd <- F1data$extractVar(c("nInd"))[[1]]
noFam <- F1data$extractVar(c("noFam"))[[1]]
bcoef_mat <- Kab <- vector(mode="list", length=noFam)
for(fam in 1:noFam){
  bcoef_mat[[fam]] <- choose(ref[[fam]]+alt[[fam]],ref[[fam]])
  Kab[[fam]] <- bcoef_mat[[fam]]*(1/2)^(ref[[fam]]+alt[[fam]])
}
seqErr=T; multiErr=T; nThreads=as.integer(3); extra = 0


OPGP <- as.integer(GUSMap:::infer_OPGP_FS(ref[[1]], alt[[1]], config[[1]], method="EM"))
MLE_EM1 <- GUSMap:::rf_est_FS(init_r = 0.01,ep = 0.001,ref=ref, alt=alt, OPGP = list(OPGP),method="EM",sexSpec = F, seqErr = T, multiErr = F)
MLE_EM2 <- GUSMap:::rf_est_FS(init_r = 0.01,ep = 0.001,ref=ref, alt=alt, OPGP = list(OPGP),method="EM",sexSpec = F, seqErr = T, multiErr = T)
MLE_EM3 <- GUSMap:::rf_est_FS(init_r = 0.01,ep = 0.001,ref=ref, alt=alt, OPGP = list(OPGP),method="EM",sexSpec = T, seqErr = T, multiErr = F)
MLE_EM4 <- GUSMap:::rf_est_FS(init_r = 0.01,ep = 0.001,ref=ref, alt=alt, OPGP = list(OPGP),method="EM",sexSpec = T, seqErr = T, multiErr = T)

#MLE_optim <- GUSMap:::rf_est_FS(init_r = 0.01,ep = 0.001,ref=ref, alt=alt, OPGP = list(OPGP),method="optim")

para <- c(GUSbase::logit2(c(MLE_EM2$rf)),GUSbase::logit(MLE_EM2$ep))
#para <- c(GUSbase::inv.logit2(rep(0.01,2*nSnps-2)),GUSbase::inv.logit(rep(0.001,1)))
#para <- c(GUSbase::logit2(rep(0.01,nSnps-1)),GUSbase::logit(rep(0.005,nSnps)))
tt <- GUSMap:::score_fs_mp_scaled_err(para,ref=ref,alt=alt,bcoef_mat=bcoef_mat,Kab=Kab,
                       nInd=nInd,nSnps=nSnps,OPGP=list(OPGP),noFam=noFam,
                       seqErr=seqErr,extra=ep,nThreads=nThreads,multiErr=T)

para <- c(GUSbase::logit2(c(MLE_EM4$rf_p,MLE_EM4$rf_m)),GUSbase::logit(MLE_EM4$ep))
para[para < -23] <- -23
tt <- GUSMap:::score_fs_mp_ss_scaled_err(para,ref=ref,alt=alt,bcoef_mat=bcoef_mat,Kab=Kab,
                                         nInd=nInd,nSnps=nSnps,OPGP=list(OPGP),noFam=noFam,
                                         seqErr=seqErr,extra=ep,nThreads=1,multiErr=T)
