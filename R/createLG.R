
createLG <- function(group, LOD, parent, LODthres, nComp, masked){

  ## Initialize the LG
  LG <- list()
  
  ## check that the parent argument is correct
  if(!is.character(parent) || length(parent) != 1 || 
     !(parent %in% c("maternal","paternal")))
    stop("parent argument is not a string of of length one or is incorrect:
         Please select one of the following:
         maternal: Only MI SNPs
         paternal: Only PI SNPs")
  
  if(parent == "maternal")
    unmapped <- sort(group$MI[which(group$MI %in% which(!masked))])
  else if(parent == "paternal")
    unmapped <- sort(group$PI[which(group$PI %in% which(!masked))])
  
  if(length(unmapped) < 2)
    stop("There are no SNPs available to create linkage groups with.")
  
  ## Run algorithm for generating the linkage groups
  finish = FALSE
  while(!finish){
    newLG <- unmapped[sort(which(LOD[unmapped,unmapped]==max(LOD[unmapped,unmapped],na.rm=T),arr.ind=T)[1,])]
    unmapped <- unmapped[-which(unmapped%in%newLG)]
    compLOD <- LOD[newLG,unmapped]
    meanLOD <- apply(compLOD,2,function(x) mean(sort(x,decreasing = T)[1:nComp], na.rm=T))
    maxMeanLOD <- max(meanLOD)
    while(maxMeanLOD > LODthres){
      newSNP <- which(meanLOD == maxMeanLOD)
      newLG <- c(newLG,unmapped[newSNP])
      unmapped <- unmapped[-newSNP]
      if(length(unmapped) > 0){
        compLOD <- matrix(LOD[newLG,unmapped],ncol=length(unmapped))
        meanLOD <- apply(compLOD,2,function(x) mean(sort(x,decreasing = T)[1:nComp], na.rm=T))
        maxMeanLOD <- max(meanLOD)
      }
      else {
        finish = TRUE
        maxMeanLOD = -999
        next
      }
    }
    meanLOD <- apply(LOD[unmapped,unmapped],2,function(x) mean(sort(x,decreasing = T)[1:nComp], na.rm=T))
    finish <- !any(meanLOD > LODthres)
    LG <- c(LG,list(newLG))
  }
  return(LG)
}
