
library(Rcpp)
Rcpp::sourceCpp("viterbi.cpp")

OPGPtoParHap <- function(OPGP){
  parHap <- matrix("A", nrow=4, ncol=length(OPGP))
  parHap[1,which(OPGP %in% c(2,4,6,8,11,12,15,16))] <- "B"
  parHap[2,which(OPGP %in% c(1,3,5,7,11,12,15,16))] <- "B"
  parHap[3,which(OPGP %in% c(3,4,7,8,10,12,14,16))] <- "B"
  parHap[4,which(OPGP %in% c(1,2,7,8,9,11,14,16))] <- "B"
  return(parHap)
}

writeLM = function(BCobj, file, direct = "./", LG = NULL, what = NULL, inferGeno = TRUE){
  ## do some checks
  if(BCobj$.__enclos_env__$private$noFam == 1){
    if(!is.vector(file) || !is.character(file) || length(file) != 1)
      stop("Filename input is invalid")
    if(!is.character(direct) || !is.vector(direct) || length(direct) != 1)
      stop("Invalid input for the path to the directory where file is to be written.")
    filename <- paste0(tail(strsplit(file,split=.Platform$file.sep)[[1]],1),"_GUSMap.txt")
    outpath <- dts(normalizePath(direct, winslash=.Platform$file.sep, mustWork=T))
    outfile = file.path(outpath,filename)
    if(!file.create(outfile,showWarnings = F))
      stop("Unable to create output file.")
    ## check for NULL what input
    if(is.null(what)){
      if(!is.null(BCobj$.__enclos_env__$private$LG_map)) what = "map"
      else if(!is.null(BCobj$.__enclos_env__$private$LG_mat_bi) || !is.null(BCobj$.__enclos_env__$private$LG_pat_bi)) what = "LG-bi"
    }
    ## Find out which LGs to write to file
    if(what == "map"){
      if(is.null(BCobj$.__enclos_env__$private$LG_map) && is.null(BCobj$.__enclos_env__$private$para))
        stop("No linkage maps available to write to file.")
      if(is.null(LG)){
        LGlist <- BCobj$.__enclos_env__$private$LG_map
        LG <- 1:length(BCobj$.__enclos_env__$private$LG_map)
      } else if(isValue(LG, type="pos_integer", minv=1, maxv=length(BCobj$.__enclos_env__$private$LG_map)))
        stop("LGs to write to file is invalid.")
      else
        LGlist <- BCobj$.__enclos_env__$private$LG_map[LG]
      if(!is.vector(inferGeno) || length(inferGeno) != 1 || !is.logical(inferGeno) || is.na(inferGeno))
        stop("Argument `inferGeno` is invalid. Should be a logical value")
      if(length(unlist(LGlist)) == 0)
        stop("No maps have been computed for the specified linkage groups.")
      ## Infer genotypes if required
      if(inferGeno){
        nInd = unlist(BCobj$.__enclos_env__$private$nInd)
        OPGPs = BCobj$.__enclos_env__$private$para$OPGP[LG]
        matgroups <- which(unlist(lapply(OPGPs[LG], function(x) all(x %in% c(1:4,9:12)))))
        patgroups <- which(unlist(lapply(OPGPs[LG], function(x) all(x %in% c(1:8)))))
        geno = list()
        for(lg in 1:length(LG)){
          snps = LGlist[[lg]]
          if(length(snps) == 0)
            next
          else{
            nSnps = length(snps)
            ref = BCobj$.__enclos_env__$private$ref[[1]][,snps]
            alt = BCobj$.__enclos_env__$private$alt[[1]][,snps]
            # output from viterbi_fs_err: 0 = 00, 1 = 10, 2 = 01, 3 = 11 (mat/pat)
            ms = viterbi_fs_err(BCobj$.__enclos_env__$private$para$rf[[lg]], 
                                rep(BCobj$.__enclos_env__$private$para$ep[[lg]], length.out=nSnps),
                                nInd, nSnps, ref, alt,
                                OPGPs[[lg]])
            ## 1 = on paternal chrosome, 0 = maternal chromsome
            infer_pat = (ms > 2) + 1
            infer_mat = apply(ms, 2, function(x) x %in% c(1,3)) + 1
            parHap = OPGPtoParHap(OPGPs[[lg]])
            geno_temp = sapply(1:nInd, function(x) rbind(parHap[cbind(infer_mat[x,]+2,1:nSnps)], parHap[cbind(infer_pat[x,],1:nSnps)]), simplify = FALSE)
            geno[[lg]] = t(do.call("rbind", geno_temp))
            #geno_pat = sapply(1:nInd, function(x) parHap[cbind(infer_pat[x,]+1,1:nSnps)])
            #geno_mat = sapply(1:nInd, function(x) parHap[cbind(infer_mat[x,]+3,1:nSnps)])
          }
        }
      }
      ## compute the information to go in file
      nLG = length(LGlist)
      LGno <- rep(LG, unlist(lapply(LGlist, length)))
      LGindx <- sequence(lapply(LGlist, unlist(length)))
      chrom <- BCobj$.__enclos_env__$private$chrom[unlist(LGlist)]
      pos <- BCobj$.__enclos_env__$private$pos[unlist(LGlist)]
      depth <- colMeans(BCobj$.__enclos_env__$private$ref[[1]][,unlist(LGlist)] + BCobj$.__enclos_env__$private$alt[[1]][,unlist(LGlist)])
      config_temp <- BCobj$.__enclos_env__$private$config[[1]]
      config_temp[!is.na(BCobj$.__enclos_env__$private$config_infer[[1]])] <- BCobj$.__enclos_env__$private$config_infer[[1]][!is.na(BCobj$.__enclos_env__$private$config_infer[[1]])] 
      segType <- (c("BI","PI","PI","MI","MI","U"))[config_temp[unlist(LGlist)]]
      temp <- BCobj$.__enclos_env__$private$ref[[1]][,unlist(LGlist)] > 0 | BCobj$.__enclos_env__$private$alt[[1]][,unlist(LGlist)] > 0
      callrate <- colMeans(temp)
      matgroups <- which(unlist(lapply(BCobj$.__enclos_env__$private$para$OPGP[LG], function(x) all(x %in% c(1:4,9:12)))))
      patgroups <- which(unlist(lapply(BCobj$.__enclos_env__$private$para$OPGP[LG], function(x) all(x %in% c(1:8)))))
      rf_m <- c(unlist(lapply(BCobj$.__enclos_env__$private$para$rf[matgroups], function(y) {if(!is.na(y[1])) return(format(round(cumsum(c(0,y)),6),digits=6,scientific=F))} )),
                rep(0,length(unlist(unlist(LGlist[patgroups])))))
      rf_p <- c(rep(0,length(unlist(unlist(LGlist[matgroups])))),
                unlist(lapply(BCobj$.__enclos_env__$private$para$rf[patgroups], function(y) {if(!is.na(y[1])) return(format(round(cumsum(c(0,y)),6),digits=6,scientific=F))} )))
      ep <- format(rep(round(unlist(BCobj$.__enclos_env__$private$para$ep[LG]),8), unlist(lapply(LGlist, length))),digits=8,scientific=F)
      ## add geno information if required
      if(inferGeno){
        temp = do.call("rbind", geno)
        temp2 = split(temp, rep(1:ncol(temp), each = nrow(temp)))
        names(temp2) = paste0(c("MAT_","PAT_"), rep(BCobj$.__enclos_env__$private$indID[[1]], rep(2,nInd)))
        out <- c(list(LG=LGno, LG_POS=LGindx, CHROM=chrom, POS=pos, TYPE=segType, RF_PAT=rf_p, RF_MAT=rf_m,
                    ERR=ep, MEAN_DEPTH=depth, CALLRATE=callrate), temp2)
      } else 
        out <- list(LG=LGno, LG_POS=LGindx, CHROM=chrom, POS=pos, TYPE=segType, RF_PAT=rf_p, RF_MAT=rf_m,
                    ERR=ep, MEAN_DEPTH=depth, CALLRATE=callrate)
      ## write out file
      data.table::fwrite(out, file = outfile, sep="\t", nThread = 1)  
      cat(paste0("Linkage analysis results written to file:\nFilename:\t",filename,"\n")) 
      return(invisible(NULL))
    } else if(what == "LG-bi"){
      stop("Yet to be implemented")
      if(is.null(BCobj$.__enclos_env__$private$LG_mat_bi) && is.null(BCobj$.__enclos_env__$private$LG_pat_bi))
        stop("No ordered linkage groups available to write to file.")
      LGlist <- c(BCobj$.__enclos_env__$private$LG_mat_bi,BCobj$.__enclos_env__$private$LG_pat_bi)
      if(is.null(LG)){
        LG <- 1:length(LGlist)
      } else if(isValue(LG, type="pos_integer", minv=1, maxv=length(LGlist)))
        stop("LGs to write to file is invalid.")
      else
        LGlist <- LGlist[LG]
      ## compute the information to go in file
      nLG = length(LGlist)
      LGno <- rep(LG, unlist(lapply(LGlist, length)))
      LGindx <- sequence(lapply(LGlist, unlist(length)))
      chrom <- BCobj$.__enclos_env__$private$chrom[unlist(LGlist)]
      pos <- BCobj$.__enclos_env__$private$pos[unlist(LGlist)]
      depth <- colMeans(BCobj$.__enclos_env__$private$ref[[1]][,unlist(LGlist)] + BCobj$.__enclos_env__$private$alt[[1]][,unlist(LGlist)])
      segType <- (c("BI","PI","PI","MI","MI","U"))[BCobj$.__enclos_env__$private$config[[1]][unlist(LGlist)]]
      temp <- BCobj$.__enclos_env__$private$ref[[1]][,unlist(LGlist)] > 0 | BCobj$.__enclos_env__$private$alt[[1]][,unlist(LGlist)] > 0
      callrate <- colMeans(temp)
      out <- list(LG=LGno, LGPOS=LGindx, CHROM=chrom, POS=pos, TYPE=segType, MEAN_DEPTH=depth, CALLRATE=callrate)
      ## write out file
      data.table::fwrite(out, file = outfile, sep="\t", nThread = 1)  
      cat(paste0("Linkage analysis results written to file:\nFilename:\t",filename,"\n")) 
      return(invisible(NULL))
    } else
      stop("Invalid what input.")
  }
}

writeLM(mySpecies, file="test_bc")







