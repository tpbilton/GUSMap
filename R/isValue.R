


isValue <- function(x, type="integer", min=0, max=NULL){
  if(type=="integer" & !is.null(max))
    return(!is.vector(x) || !is.numeric(x) || any(is.na(x)) || !all(x == round(x)) || any(x < min) || any(x > max))
  else if(type=="integer" & is.null(max))
    return(!is.vector(x) || !is.numeric(x) || any(is.na(x)) || !all(x == round(x)) || any(x < min) || is.infinite(x))
}