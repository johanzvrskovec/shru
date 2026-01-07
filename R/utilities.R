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

#based on https://stackoverflow.com/questions/19655579/a-function-that-returns-true-on-na-null-nan-in-r
#note - not vectorised; behaves as is.null()
is.invalid <- function(x, false.triggers=FALSE){
  if(is.function(x)) return(FALSE) # Some of the tests below trigger
  # warnings when used on functions
  return(
    is.null(x) ||                # Actually this line is unnecessary since
      length(x) == 0 ||            # length(NULL) = 0, but I like to be clear
      all(is.na(x)) ||
      all(x=="") ||
      (false.triggers && all(!x))
  )
}

