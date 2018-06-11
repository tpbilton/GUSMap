#### R script contains functions for performing Monte Carlo simulations
### Author: Timothy Bilton
### Date: 20/02/17
### Edited: 25/08/17


## Function for simulating sequencing data for a single full-sib family
simFS <- function(rVec_f, rVec_m=rVec_f, epsilon=0, config, nInd, meanDepth, thres=NULL, noFam=1, noChr, NoDS=1,
                  formats=list(gusmap=F,onemap=F,lepmap=F,joinmap=F,crimap=F), rd_dist="Neg_Binom",
                  filename="sim", direct="./", seed1=1, seed2=1){
  
  ## perform some checks for data input
  if(!is.numeric(noFam) || length(noFam) != 1 || noFam < 1 || is.finite(noFam))
    stop("Number of families parameter must be a single finite positive numeric number")
  if(!is.numeric(noChr) || length(noChr) != 1 || noChr < 1 || is.finite(noChr))
    stop("Number of chromosomes parameter must be a single finite positive numeric number")
  if( !is.numeric(rVec_f) || !is.numeric(rVec_m) || any(rVec_f < 0) || any(rVec_m < 0) ||
      any(rVec_f > 0.5) || any(rVec_m > 0.5) )
    stop("Recombination factions are required to be a numeric number between 0 and 0.5")
  if( !is.numeric(epsilon) || length(epsilon) != 1 || epsilon < 0 || epsilon > 1 )
    stop("Sequencing error parameter is not single numeric number between 0 and 1")
  if(length(nInd) != noFam || any(!is.numeric(nInd)) || any(nInd < 1) || any(nInd != round(nInd)) || any(!is.finite(nInd)) )
    stop("Number of individuals or number of SNPs are not a positive integer")
  if(length(config) != noChr || !is.list(config))
    stop("Segregation information needs to be a list equal to the nunber of chromosome.")
  else{
    for(chr in 1:onChr){
      if( length(config[[chr]]) != noFam || !is.list(config[[chr]]))
        stop(paste0("Segregation information for chromosome ",chr," is not equal to the number of families")
      if(!is.numeric(config[[chr]]) || !is.vector(config[[chr]]) || any(!(config[[chr]] == round(config[[chr]]))) || any(config[[chr]] < 1) || any(config[[chr]] > 9) )
        stop(paste0("Segregation information for chromosome ",chr," in familiy ",fam," needs to be an integer vector equal to the number of SNPs with entires from 1 to 9")
    }
  }
  if( !is.numeric(meanDepth) || meanDepth <= 0 || !is.finite(meanDepth) )
    stop("The mean of the read depth distribution is not a finitie positive number.")
  if( !is.numeric(NoDS) || NoDS < 1 || NoDS != round(NoDS) || !is.finite(NoDS))
    stop("Imput for the number of data sets needs to be a finite positive number")
  if( !is.null(thres) & !is.numeric(thres) )
    stop("The read depth threshold value is not a finite numeric number")
  if( !is.numeric(seed1) || !is.numeric(seed2) )
    stop("Seed values for the randomziation need to be numeric values")
  
  ## Compute the number of SNPs
  nSnps <- lapply(config, function(x) as.integer(length(x)))
  
  ## booleen variable to indicate whether to write out files
  writeFiles <- any(c(isTRUE(formats$gusmap),isTRUE(formats$onemap),isTRUE(formats$lepmap),isTRUE(formats$crimap),isTRUE(formats$joinmap)))
  
  ## If to write files: check that file directory is in correct format
  if(writeFiles){
    if(!is.character(filename) || length(filename)!=1)
      stop("Name to write data sets to need to be a string of length 1")
    if(!is.character(direct) || length(direct)!=1)
      stop("Name to write data sets to need to be a string of length 1")
  }
  
  OPGP <- replicate(n = noChr, vector(mode = "list", length = noFam), simplify = F)
  genon <- ref <- alt <- vector(mode = "list", length=noFam)
  for(chr in 1:noChr){
    for(fam in 1:noFam){
  
      ## Create list of the recombination fraction and 1 minus the recombination fraction
      ## for each SNP
      rVec_f <- sapply(rep(rVec_f,length.out=nSnps[[chr]]-1), function(r) c(r,1-r),simplify=F)
      rVec_m <- sapply(rep(rVec_m,length.out=nSnps[[chr]]-1), function(r) c(r,1-r),simplify=F)
      
      ## simulate the parental haplotypes
      set.seed(seed1)
      parHap <- matrix(rep(paste0(rep("A",nSnps[[chr]])),4),nrow=4)
      parHap[cbind(sample(1:2,size=sum(config[[chr]][[fam]] %in% c(2,3)),replace=T),which(config[[chr]][[fam]]%in%c(2,3)))] <- "B"
      parHap[cbind(sample(3:4,size=sum(config[[chr]][[fam]] %in% c(4,5)),replace=T),which(config[[chr]][[fam]] %in% c(4,5)))] <- "B"
      parHap[cbind(c(rep(3,sum(config[[chr]][[fam]] %in% c(3,7,9))),rep(4,sum(config[[chr]][[fam]] %in% c(3,7,9)))),which(config[[chr]][[fam]] %in% c(3,7,9)))] <- "B"
      parHap[cbind(c(rep(1,sum(config[[chr]][[fam]] %in% c(5,8,9))),rep(2,sum(config[[chr]][[fam]] %in% c(5,8,9)))),which(config[[chr]][[fam]] %in% c(5,8,9)))] <- "B"
      parHap[cbind(c(sample(1:2,size=sum(config[[chr]][[fam]]==1),replace=T),sample(3:4,size=sum(config[[chr]][[fam]]==1),replace=T)),rep(which(config[[chr]][[fam]]==1),2))] <- "B"
      if(any(parHap[1:2,] == "B") && parHap[1,which(apply(parHap[1:2,],2,function(x) !(all(x=='A'))))[1]] == 'B')
        parHap[1:2,] <- parHap[2:1,]
      if(any(parHap[3:4,] == "B") && parHap[3,which(apply(parHap[3:4,],2,function(x) !(all(x=='A'))))[1]] == 'B')
        parHap[3:4,] <- parHap[4:3,]
      OPGP[[chr]][[fam]] <- parHapToOPGP(parHap)
      
      #### Simulate the data sets
      set.seed(seed2)
      for(sim in 1:NoDS){
        
        #### Simulate the true Meiosis for each individual at each SNP.
        mIndx <- matrix(c(sample(c(0,1),size=2*nInd[fam],replace=T)),ncol=1)
        for(i in 1:(nSnps[[chr]]-1)){
          newmIndx_f <- numeric(nInd[fam])
          newmIndx_m <- numeric(nInd[fam])
          for(j in 1:(nInd[fam])){
            newmIndx_f[j] <- sample(c(0,1),size=1,prob=c(rVec_f[[i]][(mIndx[j,i]==0)+1],rVec_f[[i]][(mIndx[j,i]==1)+1]))
            newmIndx_m[j] <- sample(c(0,1),size=1,prob=c(rVec_m[[i]][(mIndx[j+nInd[fam],i]==0)+1],rVec_m[[i]][(mIndx[j+nInd[fam],i]==1)+1]))
          }
          mIndx <- cbind(mIndx,c(newmIndx_f,newmIndx_m))
        }
        # Determine the true genotype calls
        geno <- rbind(sapply(1:nSnps[[chr]],function(x) parHap[mIndx[1:nInd[fam],x]+1,x]),sapply(1:nSnps[[chr]],function(x) parHap[mIndx[1:nInd[fam]+nInd[fam],x]+3,x]))
        geno <- sapply(1:nSnps[[chr]],function(y) {
          tempGeno <- geno[,y]
          sapply(1:nInd[fam], function(x) paste(sort(c(tempGeno[x],tempGeno[x+nInd[fam]])),collapse=""))
        })
        geno <- (geno=="AA")*2 + (geno=="AB")*1
        
        ### Now generate the sequencing data
        # 1: Simulate Depths
        depth <- matrix(0,nrow=nInd[fam], ncol=nSnps[[chr]])
        if(rd_dist=="NegBinom")
          depth[which(!is.na(geno))] <- rnbinom(sum(!is.na(geno)),mu=meanDepth,size=2) 
        else   
          depth[which(!is.na(geno))] <- rpois(sum(!is.na(geno)),meanDepth)
        # 2: simulate sequencing genotypes (with sequencing error rate of epsilon)
        aCounts <- matrix(rbinom(nInd[fam]*nSnps[[chr]],depth,geno/2),ncol=nSnps[[chr]])
        bCounts <- depth - aCounts
        aCountsFinal <- matrix(rbinom(nInd[fam]*nSnps[[chr]],aCounts,prob=1-epsilon),ncol=nSnps[[chr]]) + matrix(rbinom(nInd[fam]*nSnps[[chr]],bCounts,prob=epsilon),ncol=nSnps[[chr]])
        SEQgeno <- aCountsFinal/depth
        SEQgeno[which(SEQgeno^2-SEQgeno<0)] <- 0.5
        SEQgeno <- 2* SEQgeno  ## GBS genotype call
          
        ## Write data to file
        if(writeFiles)
          genoToOtherFormats(SEQgeno,aCountsFinal,depth-aCountsFinal,config[[chr]][[fam]],formats=formats,filename=paste0(filename,sim,"_fam",fam,"_chr",chr),direct=direct,thres=thres,sim=sim)
        
      }
      ## Write simulation parameters to a file
      if(writeFiles)
        dput(list(nInd=nInd,nSnps=nSnps,NoDS=NoDS,rVec_f=unlist(lapply(rVec_f,function(x) x[1])),
                  rVec_m=unlist(lapply(rVec_m,function(x) x[1])),
                  config=config,OPGP=OPGP,meanDepth=meanDepth,rd_dist=rd_dist),
             paste0(trim_fn(paste0(direct,"/",filename)),"_fam",fam,"_chr",chr,"info.txt"))
      ## return simulated data and parameter values ued to generate the data
      else{
        genon[[fam]] <- cbind(genon[[fam]], SEQgeno)
        ref[[fam]]   <- cbind(ref[[fam]], aCountsFinal)
        alt[[fam]]   <- cbind(alt[[fam]], matrix(as.integer(depth-aCountsFinal), nrow=nInd, ncol=nSnps[[chr]]))
      }
    }
  }
  if(!writeFiles){
    chrom <- sapply(1:noChr, function(x) rep(x,nSnps[[x]]), simplify = F)
    pos  <-  sapply(1:noChr, function(x) 1:(nSnps[[x]]), simplify = F)
    nSnps <- as.integer(unlist(nSnps))

    obj <- RA$new(
      list(genon = genon, ref = ref, alt = alt, chrom = chrom, pos = pos,
           SNP_Names = NULL, indID = 1:nSnps, nSnps = nSnps, nInd = lapply(as.list(nInd), as.integer), gform = "simFS", AFrq = NULL)
    )
    newObj <- FS$new(obj)
    return(newObj)
  }
}

### Function for writing simulated sequencing data to various software formats
genoToOtherFormats <- function(genon,depth_Ref,depth_Alt,config,formats,filename,direct,thres=NULL,ratioThres=TRUE,sim=sim){
  
  ## specify which formats to use
  gusmap <- isTRUE(formats$gusmap)
  onemap <- isTRUE(formats$onemap)
  lepmap <- isTRUE(formats$lepmap)
  joinmap <- isTRUE(formats$joinmap)
  crimap <- isTRUE(formats$crimap)
  
  depth <- depth_Ref + depth_Alt
  
  ## Write data to gusmap format
  newfile <- paste0(trim_fn(direct),"/",trim_fn(filename))
  if(gusmap){
    write.table(genon,paste0(newfile,"_genon_SEQ.txt"),row.names=F,col.names=F)
    write.table(depth_Ref,paste0(newfile,"_depth_Ref_SEQ.txt"),row.names=F,col.names=F)
    write.table(depth_Alt,paste0(newfile,"_depth_Alt_SEQ.txt"),row.names=F,col.names=F)
  }
  
  ## write data to other formats is required
  if(any(onemap,lepmap,joinmap,crimap)){
    
    nSnps <- ncol(genon); nInd <- nrow(genon)
    
    # Set genotypes below certain depths to zero
    if(!is.null(thres) & is.numeric(thres) & length(thres) == 1)
      genon[which(depth < max(1,thres))] <- NA
    # Set genotypes with distortion of allele ratios
    if(ratioThres){
      genon[which( (depth_Ref > 9 & depth_Alt == 1) | (depth_Alt > 9 & depth_Ref == 1))] <- NA
    }
      
    # specify the genotypes
    AAgeno <- which(genon==2,arr.ind=T)
    BBgeno <- which(genon==0,arr.ind=T)
    badBB <- matrix(BBgeno[which(BBgeno[,2] %in% (which(config %in% c(2,4)))),],ncol=2,byrow=F)
    badAA <- matrix(AAgeno[which(AAgeno[,2] %in% (which(config %in% c(3,5)))),],ncol=2,byrow=F)
    # Check where there are any partially informative markers with BB or AA
    genon[badBB] <- 1
    genon[badAA] <- 1
    AAgeno <- which(genon==2)
    BBgeno <- which(genon==0)
    ABgeno <- which(genon==1)
    NAgeno <- which(is.na(genon))
    
    #### OneMap
    if(onemap){
      onemapData <- genon
      
      ## change the format of the genotypes for onemap
      onemapData[NAgeno] <- "-"
      onemapData[AAgeno] <- 'a'
      onemapData[ABgeno] <- 'ab'
      onemapData[BBgeno] <- 'b'
      
      ## can't have any BB's in onemap's data set
      for(snp in 1:ncol(onemapData)){
        if(config[snp] %in% c(3,5))
          onemapData[which(onemapData[,snp]=="b"),snp] <- "a"
      }
            
      # Vector giving the three different marker types
      mType <- c('B3.7','D1.10','D1.10','D2.15','D2.15')
      
      ## Write out the onemap file
      cat(nInd,' ',nSnps,'\n',
          sapply(1:nSnps,function(x) paste('*M',x,' ',mType[config[x]],'\t',paste(onemapData[,x],collapse=","),'\n',sep="")),
          sep="",file=paste0(newfile,'_OneMap.txt'))
      
    }
    
    #### LepMap and Crimap
    if(lepmap|crimap){
      
      # copy dataset to a new data set
      newGenon <- genon
      
      ## Compute the offsprings genons
      newGenon[NAgeno] <- "0 0"
      newGenon[AAgeno] <- "1 1"
      newGenon[ABgeno] <- "1 2"
      newGenon[BBgeno] <- "2 2"
      
      parGenon <- sapply(1:length(config),function(x) {
        if(config[x] == 1){
          p1 <- "1 2"
          p2 <- "1 2"
        }
        else if(config[x] == 2){
          p1 <- "1 2"
          p2 <- "1 1"
        }
        else if(config[x] == 3){
          p1 <- "1 2"
          p2 <- "2 2"
        }
        else if(config[x] == 4){
          p1 <- "1 1"
          p2 <- "1 2"
        }
        else if(config[x] == 5){
          p1 <- "2 2"
          p2 <- "1 2"
        }
        else if(config[x] == 6){
          p1 <- "1 1"
          p2 <- "1 1"
        }
        else if(config[x] == 7){
          p1 <- "1 1"
          p2 <- "2 2"
        }
        else if(config[x] == 8){
          p1 <- "2 2"
          p2 <- "1 1"
        }
        else if(config[x] == 9){
          p1 <- "2 2"
          p2 <- "2 2"
        }
        return(c(p1,p2))
      })
      
      if(lepmap) {#### Output file for lepmap
        # write to output
        LPout <- cbind(rep("FS",(nInd+2)),1:(nInd+2),c(0,0,rep(1,nInd)),c(0,0,rep(2,nInd)),c(1,2,rep(0,nInd)),rep(0,nInd+2),
                       rbind(parGenon,newGenon))
        write.table(x=LPout,file=paste0(newfile,"_LepMap.txt"),
                    sep="\t",col.names=F,row.names=F,quote=F)
      }
      
      if(crimap){ #### Output file for crimap
        newfile2 <- paste0(trim_fn(paste0(direct,"/chr1_",filename)),".gen")
        cat("1\n",nSnps,"\n",paste0("M",1:nSnps,collapse=" "),"\n1\n",nInd+2,"\n",sep="",file=newfile2)
        # Add Father
        cat(paste(nInd+1,0,0,0),"\n",
            paste(parGenon[2,],collapse=" "),"\n",sep="",file=newfile2,append=T)
        # Add Mother
        cat(paste(nInd+2,0,0,1),"\n",
            paste(parGenon[1,],collapse=" "),"\n",sep="",file=newfile2,append=T)
        # add offspring
        cat(paste0(paste(1:nInd,nInd+1,nInd+2,3,"\n"),sapply(1:nInd,function(x) paste(paste(newGenon[x,],collapse=" "),"\n")),sep=""),
            sep="", append=T,file=newfile2)
      }
    }
    
    #### JoinMap
    if(joinmap){
      
      # copy dataset to a new data set
      joinmapData <- genon
      usnps <- which(apply(joinmapData,2,function(x) all(is.na(x))))
      
      ## convert the genotype calls to joinmap format
      for(snp in 1:nSnps){
        newCol <- character(nInd)
        if(config[snp] == 1){
          newCol[which(genon[,snp]==2)] <- 'hh'
          newCol[which(genon[,snp]==1)] <- 'hk'
          newCol[which(genon[,snp]==0)] <- 'kk'
          newCol[which(is.na(genon[,snp]))] <- '--'
          joinmapData[,snp] <- newCol
        }
        else if(config[snp] == 2){
          newCol[which(genon[,snp]==2)] <- 'nn'
          newCol[which(genon[,snp]==1)] <- 'np'
          newCol[which(genon[,snp]==0)] <- 'np' 
          newCol[which(is.na(genon[,snp]))] <- '--'
          joinmapData[,snp] <- newCol
        }
        else if(config[snp] == 3){
          newCol[which(genon[,snp]==2)] <- 'np'
          newCol[which(genon[,snp]==1)] <- 'np'
          newCol[which(genon[,snp]==0)] <- 'nn' 
          newCol[which(is.na(genon[,snp]))] <- '--'
          joinmapData[,snp] <- newCol
        }
        else if(config[snp] == 4){
          newCol[which(genon[,snp]==2)] <- 'll'
          newCol[which(genon[,snp]==1)] <- 'lm'
          newCol[which(genon[,snp]==0)] <- 'lm'
          newCol[which(is.na(genon[,snp]))] <- '--'
          joinmapData[,snp] <- newCol
        }
        else if(config[snp] == 5){
          newCol[which(genon[,snp]==2)] <- 'lm'
          newCol[which(genon[,snp]==1)] <- 'lm'
          newCol[which(genon[,snp]==0)] <- 'll'
          newCol[which(is.na(genon[,snp]))] <- '--'
          joinmapData[,snp] <- newCol
        }
      }

      ##form the first 4 lines of the file
      cat('name = in.loc\n','popt = CP\n','nloc = ',nSnps-length(usnps),'\n','nind = ',nInd,'\n\n',
          file=paste0(newfile,'_JoinMap.loc'),sep="")
      
      cat(sapply((1:nSnps)[which(!(1:nSnps %in% usnps))], function(x) {
        paste0('M',x,'  ',switch(config[x],'<hkxhk>','<nnxnp>','<nnxnp>','<lmxll>','<lmxll>'),'\n',paste(joinmapData[,x],collapse=" "),'\n')
      }), file=paste0(newfile,'_JoinMap.loc'),sep="",append=T)
    }
  }
  return(invisible())
}

## Function for trimming away white spaces and forward slashes for input path and file names
trim_fn <- function(x) return( gsub("^\\/|\\/$", "", trimws(x)) )


