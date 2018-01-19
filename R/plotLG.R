#### Function for plotting linkage groups (or a single linkage group)
plotLG <- function(mat, LG, filename=NULL, names=NULL, chrS=2, lmai=2, chrom=T){

  b <- ncol(mat) + 1
  if(chrom)
    chrom.ind <- unlist(lapply(LG, function(x) c(x,b)))[-length(unlist(LG))+length(LG)]
  else
    chrom.ind <- 1:ncol(mat)

  ## Subset the matrix
  mat <- cbind(mat,rep(0,b-1))
  mat <- rbind(mat,rep(0,b))
  mat <- mat[chrom.ind,chrom.ind]
  ## work out where the breaks are
  breaks <- which(chrom.ind==b)
  npixels <- length(chrom.ind)
  if(chrom){
    if(!is.null(filename))
      png(filename,width=npixels+72*lmai,height=npixels,res=72)
    par(xaxt='n',yaxt='n',mai=c(0,lmai,0,0),bty='n',ann=F)
    image(1:npixels,1:npixels,mat,zlim=c(0,0.5),col=heat.colors(100))
    if(is.null(names))
      mtext(paste("LG",1:length(LG),"  "),
                  at=floor(apply(cbind(c(0,breaks),c(breaks,npixels)),1,median)),side=2, line=0,cex=chrS,las=1)
    else
      mtext(names, at=floor(apply(cbind(c(0,breaks),c(breaks,npixels)),1,median)),side=2, line=0,cex=chrS,las=1)
    abline(h=breaks)
    abline(v=breaks)
    if(!is.null(filename))
      dev.off()
  }
  else{
    npixels <- length(chrom.ind)
    if(!is.null(filename))
      png(filename,width=npixels,height=npixels)
    par(xaxt='n',yaxt='n',mar=c(0,0,0,0),bty='n',ann=F)
    image(1:npixels,1:npixels,mat,zlim=c(0,0.5),col=heat.colors(100))
    abline(h=breaks)
    abline(v=breaks)
    if(!is.null(filename))
      dev.off()
  }
}
