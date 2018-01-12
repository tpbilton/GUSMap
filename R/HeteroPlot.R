

HeteroPlot <- function(depth_Ref, depth_Alt, model="random", alpha=NULL, filename="HeteroPlot", cex=1){

  ## Run some tests on the input variables
  if(model == "beta-binomal" && (alpha < 0 | !is.numeric(alpha)) )
     stop("Parameter for the beta-binomal model is not a positive numeric value.")
  
  ## Compute the observed depth count
  heteCall <- which(depth_Ref != 0 & depth_Alt != 0)
  hete <- cbind(depth_Ref[heteCall], depth_Alt[heteCall])
  hete <- cbind(apply(hete,1,min),apply(hete,1,max))
  ## compute the expected depth count
  heteDepth <- rowSums(hete)
  expHete <- matrix(nrow=length(heteDepth), ncol=2)
  ### simulated the expected ratios
  for(i in 1:length(heteDepth)){
    if(model == "beta-binom"){
      newCall <- 0
      while((newCall %% heteDepth[i]) == 0){
        p <- rbeta(1, alpha, alpha)
        newCall <- rbinom(1, heteDepth[i], p)
      }
      expHete[i,] <- sort(c(newCall,heteDepth[i]-newCall))
    }
    else {
      newCall <- c(0,0)
      while(any(newCall == 0))
        newCall <- table(  factor(rbinom(heteDepth[i], size=1, prob=0.5),levels=0:1,labels=c("0","1")))
      expHete[i,] <- sort(newCall)
    }
  }
  maxCount <- max(hete,expHete,100)
  obsHeteTab <- table(factor(hete[,1],levels=1:maxCount,labels=1:maxCount),
                      factor(hete[,2],levels=1:maxCount,labels=1:maxCount))
  obsHeteTab <- log(obsHeteTab)
  expHeteTab <- table(factor(expHete[,1],levels=1:maxCount,labels=1:maxCount),
                      factor(expHete[,2],levels=1:maxCount,labels=1:maxCount))
  expHeteTab <- log(expHeteTab)
  
  HeteMat <- matrix(nrow=maxCount+1,ncol=maxCount+1)
  HeteMat[lower.tri(HeteMat)] <- t(expHeteTab)[-which(upper.tri(t(expHeteTab)))]
  HeteMat[upper.tri(HeteMat)] <- obsHeteTab[-which(lower.tri(t(expHeteTab)))]
  HeteMat[which(is.infinite(HeteMat))] <- NA
  
  ### Produce the plot
  png(paste0(filename,".png"), width=maxCount*3+100,height=maxCount*3+100)
  par(mar = c(5.1,5.1,4.1,3.1))
  newCol <- colorRampPalette(c("red","orange","yellow","green","cyan","blue"))(200)
  image(1:nrow(HeteMat),1:ncol(HeteMat),HeteMat, xlab="No Reads",ylab="No Reads",col=newCol,zlim=c(0,max(HeteMat,na.rm=T)),
        cex.lab=cex,cex.axis=cex)
  abline(0,1)
  abline(0,4, lty=2)
  mtext("Observed",side=3,cex=cex, font=2, line=1)
  mtext(paste0("Expected (",ifelse(model=="beta-binom","Beta-Binomial","Random")," model)"), side=4,cex=cex, font=2,line=1)
  legend_image <- as.raster(matrix(rev(newCol), ncol = 1))
  rasterImage(legend_image,xleft = maxCount-20, ybottom = 10, ytop = 60, xright = maxCount-10 )
  text(x = maxCount-28, y = 10 + 10 * 0:5, labels = format(exp(seq(0,max(HeteMat,na.rm=T),length.out=6)),digits=1))
  dev.off()
  return(invisible())
}