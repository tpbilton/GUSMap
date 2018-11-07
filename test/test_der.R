set.seed(6745)
config <- list(list(sample(c(1,2,4), size=30, replace=TRUE)))
F1data <- GUSMap:::simFS(0.01, config=config, meanDepth=10, nInd=5, epsilon=0.005)

ref <- F1data$extractVar(list("ref"))[[1]]
alt <- F1data$extractVar(list("alt"))[[1]]
config <- F1data$extractVar(list("config"))[[1]]
nSnps <- F1data$extractVar(list("nSnps"))[[1]]
nInd <- F1data$extractVar(list("nInd"))[[1]]
noFam <- F1data$extractVar(list("noFam"))[[1]]
bcoef_mat <- Kab <- vector(mode="list", length=noFam)
for(fam in 1:noFam){
  bcoef_mat[[fam]] <- choose(ref[[fam]]+alt[[fam]],ref[[fam]])
  Kab[[fam]] <- bcoef_mat[[fam]]*(1/2)^(ref[[fam]]+alt[[fam]])
}
seqErr=T; multiErr=T; nThreads=1; extra = 0


OPGP <- as.integer(GUSMap:::infer_OPGP_FS(ref[[1]], alt[[1]], config[[1]]))

para <- c(GUSbase::inv.logit2(rep(0.01,nSnps-1)),GUSbase::inv.logit2(rep(0.005,nSnps)))
tt <- GUSMap:::score_fs_mp_scaled_err(para,ref=ref,alt=alt,bcoef_mat=bcoef_mat,Kab=Kab,
                       nInd=nInd,nSnps=nSnps,OPGP=OPGP,noFam=noFam,
                       seqErr=seqErr,extra=ep,nThreads=nThreads,multiErr=multiErr)