

mfun <- function(rf, fun="morgan", centiM=FALSE){
  if(fun == "morgan")
    dis = rf
  else if (fun == "haldane")
    dis = -0.5 * log(1 - 2 * rf)
  else if (fun == "kosambi")
    dis = 0.25 * log((1 + 2 * rf)/(1 - 2 * rf))
  else
    stop("Map function not defined")
  if(centiM)
    return(100 * dis)
  else
    return(dis)
}
  
  