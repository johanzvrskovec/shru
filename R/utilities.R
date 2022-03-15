#small numerical utilities etc

clipValues <- function(x,min,max) {
  if(!is.null(min)) x[which(x<min)]<-min
  if(!is.null(max)) x[which(x>max)]<-max
  return(x)
}


padList <- function(l,padding,targetLength){
  pl<-targetLength-length(l)
  if(pl>0) {c(l,rep(padding,pl))} else {l}
}

