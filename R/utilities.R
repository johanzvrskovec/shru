#small numerical utilities etc

clipValues <- function(x,min,max) {
  if(!is.null(min)) x[which(x<min)]<-min
  if(!is.null(max)) x[which(x>max)]<-max
  return(x)
}


