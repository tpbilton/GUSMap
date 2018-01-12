
createLGs <- function(group, LOD, parent, LODthres, nComp){
  
  ## Initialize the LGs
  LGs <- list()
  
  ## check that the parent argument is correct
  if(!is.character(parent) || length(parent) != 1 || 
     !(parent %in% c("maternal","paternal")))
    stop("parent argument is not a string of of length one or is incorrect:
         Please select one of the following:
         maternal: Only MI SNPs
         paternal: Only PI SNPs")
  
  if(parent == "maternal")
    unmapped <- sort(c(group$MI))
  else if(parent == "paternal")
    unmapped <- sort(c(group$PI))
  
  ## Run algorithm for generating the linkage groups
  finish = FALSE
  while(!finish){
    newLG <- unmapped[sort(which(LOD[unmapped,unmapped]==max(LOD[unmapped,unmapped],na.rm=T),arr.ind=T)[1,])]
    unmapped <- unmapped[-which(unmapped%in%newLG)]
    compLOD <- LOD[newLG,unmapped]
    meanLOD <- apply(compLOD,2,function(x) mean(sort(x,decreasing = T)[1:nComp], na.rm=T))
    maxMeanLOD <- max(meanLOD)
    while(maxMeanLOD > LODthres ){
      newSNP <- which(meanLOD == maxMeanLOD)
      newLG <- c(newLG,unmapped[newSNP])
      unmapped <- unmapped[-newSNP]
      compLOD <- LOD[newLG,unmapped]
      meanLOD <- apply(compLOD,2,function(x) mean(sort(x,decreasing = T)[1:nComp], na.rm=T))
      maxMeanLOD <- max(meanLOD)
    }
    meanLOD <- apply(LOD[unmapped,unmapped],2,function(x) mean(sort(x,decreasing = T)[1:nComp], na.rm=T))
    finish <- !any(meanLOD > LODthres)
    LGs <- c(LGs,list(newLG))
  }
  return(LGs)
}
