printSEMResult <- function(resultDf = NULL){
  parsedSEMResult<-parseSEMResult(resultDf)
  dot<-generateDOT(nodeDf=parsedSEMResults$nodeDf, edgeDf=parsedSEMResults$edgeDf)
  grViz(dot)
}