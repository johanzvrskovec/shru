printSEMResult <- function(resultDf = NULL){
  parsedSEMResult<-parseSEMResult(resultDf)
  dot<-generateDOT(nodeDf=parsedSEMResult$nodeDf, edgeDf=parsedSEMResult$edgeDf)
  grViz(dot)
  return(dot)
}
