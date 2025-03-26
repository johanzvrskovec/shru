#small numerical and string utilities etc

clipValues <- function(x,min,max) {
  if(!is.null(min)) x[which(x<min)]<-min
  if(!is.null(max)) x[which(x>max)]<-max
  return(x)
}

padList <- function(l,padding,targetLength){
  pl<-targetLength-length(l)
  if(pl>0) {c(l,rep(padding,pl))} else {l}
}

padListRight <- function(l,padding,targetLength){
  pl<-targetLength-length(l)
  if(pl>0) {c(l,rep(padding,pl))} else {l}
}

padStringLeft <- function(s,padding,targetLength){
  pl<-targetLength-nchar(s)
  if(pl>0) {paste0(c(rep(padding,pl),s),collapse = "")} else {s}
}

replaceNa <- function(x,repl=0,includeNaN=TRUE){
  x[is.na(x)|(includeNaN&is.nan(x))]<-repl
  return(x)
}

