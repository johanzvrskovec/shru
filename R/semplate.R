semplate <-c()
semplate$powers.of.two32bit <- 2^(0:(32 - 1))
semplate$bit.one<-intToBits(1)[1]
semplate$bit.zero<-intToBits(0)[1]

semplate$generateIndicatorLoadingPatternsFromFactorLoadings<-function(factorLoadings, increment, forceOneIndicatorLoading=T){
  nrow<-nrow(factorLoadings)
  ncol<-ncol(factorLoadings)
  
  indicatorLoadingPatterns<-c()
  factorLoadings<-base::abs(factorLoadings)
  
  if(forceOneIndicatorLoading){
    forced<-as.vector(apply(factorLoadings,MARGIN = 1,FUN = which.max))
  }
  
  cValue<-0
  while(cValue < max(factorLoadings)){
    toadd<-(factorLoadings>cValue)
    #if(!all(toadd) | !any(toadd)){
    if(forceOneIndicatorLoading) {
        for(iIndicator in 1 : nrow){
          toadd[iIndicator,forced[iIndicator]]<-TRUE
      }
    }
    
    #indicatorLoadingPatterns<-rbind(indicatorLoadingPatterns,list(unlist(as.data.frame(toadd)),use.names = F))
    indicatorLoadingPatterns<-rbind(indicatorLoadingPatterns,as.vector(x = unlist(as.data.frame(toadd)), mode = 'logical'))
    cValue<-(cValue+increment)
  }
  
  return(unique(indicatorLoadingPatterns))
}

semplate$generateIndicatorLoadingPatterns<-function(searchBitValues=c(0), indicatorLocks.include=c(), indicatorLocks.exclude=c()){
  
  indicatorLocks<-(indicatorLocks.include | indicatorLocks.exclude)
  bitLength<-ncol(indicatorLocks)*nrow(indicatorLocks)
  searchBitLength<-(bitLength-sum(indicatorLocks))

  num<-length(searchBitValues)
  totalBitValues<-c()
  totalBitIndicatorLoadings<-data.frame(matrix(NA, nrow = num, ncol = bitLength))
  colnames(totalBitIndicatorLoadings)<-NA
  
  #test
  #nSearchBitValue=1
  
  for(nSearchBitValue in 1:num){
    searchBitValue<-searchBitValues[nSearchBitValue]
    searchBits<-as.logical(intToBits(x = searchBitValue))
    #add locks
    nSearch=1
    #nTotal<-1 #test
    for (nTotal in 1:bitLength) {
      if(indicatorLocks[nTotal]){
        totalBitIndicatorLoadings[nSearchBitValue,nTotal]<-ifelse(indicatorLocks.include[nTotal],T,F)
      } else {
        totalBitIndicatorLoadings[nSearchBitValue,nTotal]<-searchBits[nSearch]
        nSearch<-nSearch+1
      }
    }
    
    totalBitsValue = packBits(x = as.logical(c(totalBitIndicatorLoadings[nSearchBitValue,],rep(FALSE, times=256-bitLength))),type = c("integer"))
    totalBitValues[nSearchBitValue]<-list(totalBitsValue)
  }
  
  return (list(searchBitValues=searchBitValues,
               totalBitValues=totalBitValues,
               indicatorLoadings=apply(X = totalBitIndicatorLoadings, MARGIN = 1, FUN = function(iv){
                 list(matrix(data=iv,ncol = ncol(indicatorLocks)))
                 }
                                       )
               ))
}

semplate$generateLavaanCFAModel<-function(allow_loading.table.indicator_factor, fix_loading.table.indicator_factor=NULL, fix_loading.table.residual_variance=NULL,orthogonal=FALSE, indicatorArgs=NULL, universalResidualLimitMin=0.001 ){
  
  #lavaan definition string
  lds<-""
  lds.residuals="
  "
  lds.factorSelf=""
  lds.factorOther=""
  
  cond = "
  "
  nFactors<-ncol(allow_loading.table.indicator_factor)
  nIndicators<-nrow(allow_loading.table.indicator_factor)
  
  if(is.null(indicatorArgs)){
    indicatorArgs=data.frame(code=row.names(allow_loading.table.indicator_factor), residualSizeLimitMax=NA_real_)
  }
  
  if(is.null(indicatorArgs$code)){indicatorArgs$code=row.names(allow_loading.table.indicator_factor)}
  if(is.null(indicatorArgs$residualSizeLimitMax)){indicatorArgs$residualSizeLimitMax=NA_real_}
  if(is.null(indicatorArgs$residualLimitMin)){indicatorArgs$residualLimitMin=universalResidualLimitMin}
  
  indicatorArgs[which(is.na(indicatorArgs$residualLimitMin)),c("residualLimitMin")]<-universalResidualLimitMin
  
  indicatorResidualPatternCoefficientLabels=c()
  indicatorResidualPatternCoefficientLabelDefinitions=c()
  
  #factor loadings, factor variances
  for(iFactor in 1:nFactors){
    lds.factor=paste0("
                      F",iFactor," =~ ")
    
    hasFactor=FALSE
    for(iIndicator in 1:nIndicators){
      if(allow_loading.table.indicator_factor[iIndicator,iFactor]==TRUE){
        
        if(hasFactor==TRUE)
          lds.factor=paste0(lds.factor,"+")
        
        if(is.null(fix_loading.table.indicator_factor)==FALSE) {
          #use fixed loadings
          if(!is.na(fix_loading.table.indicator_factor[iIndicator,iFactor])) {
            #vFixLoading="1*" #show this for both cases for clarity - old setting for booleanarguments
            vFixLoading=paste0("(",as.character(fix_loading.table.indicator_factor[iIndicator,iFactor]),")*")
            } else {
            if(hasFactor==TRUE)
              vFixLoading="" else
                vFixLoading="NA*" #Use NA loading for the first indicator
          }
        } else {
          #do not use fixed loadings
          if(hasFactor==TRUE)
            vFixLoading="" else
              vFixLoading="NA*"
        }
        
        lds.factor=paste0(lds.factor,vFixLoading,indicatorArgs$code[iIndicator])
        hasFactor=TRUE
      }
    }
    if(hasFactor==TRUE){
      lds=paste0(lds,lds.factor)
      lds.factorSelf=paste0(lds.factorSelf,"
                            F",iFactor,"~~1*F",iFactor)
    }
  }
  
  for(iIndicator in 1:nIndicators){
    if(is.null(fix_loading.table.residual_variance)){
    indicatorResidualPatternCoefficientLabels[iIndicator]=paste0('r',iIndicator)
    #indicatorResidualPatternCoefficientLabelDefinitions[iIndicator]=paste0('label(',indicatorResidualPatternCoefficientLabels[iIndicator],')')
    
    lds.residuals=paste0(lds.residuals,indicatorArgs$code[iIndicator],"~~",indicatorResidualPatternCoefficientLabels[iIndicator],"*",indicatorArgs$code[iIndicator],"
                         ")
    } else {
      lds.residuals=paste0(lds.residuals,indicatorArgs$code[iIndicator],"~~(",fix_loading.table.residual_variance[iIndicator],")*",indicatorArgs$code[iIndicator],"
                         ")
    }
  }
  
  indexFactor<-data.frame(index=c(1:nFactors))
  indexFactor<-indexFactor %>%
    inner_join(indexFactor, by = character())
  
  
  lds.factorOther=""
  
  for(iComb in 1:nrow(indexFactor)) {
    row<-indexFactor[iComb,]
    if(row$index.x==row$index.y || row$index.y < row$index.x)
      next()
    
    if(orthogonal)
      lds.factorOther=paste0(lds.factorOther,"
                           F",row$index.x,"~~0*F",row$index.y) else
                             lds.factorOther=paste0(lds.factorOther,"
                           F",row$index.x,"~~F",row$index.y)
  }
  
  
  if(is.null(fix_loading.table.residual_variance)){
    cond.residualSizeLimit=""
    for(iIndicator in 1:nIndicators){
      if(!is.na(indicatorArgs$residualSizeLimitMax[iIndicator])){
      #cond.residualSizeLimit=paste0(cond.residualSizeLimit,'abs(',indicatorResidualPatternCoefficientLabels[iIndicator],')','<',indicatorArgs$residualSizeLimit[iIndicator],'
                #                    ')
        cond.residualSizeLimit=paste0(cond.residualSizeLimit,'abs(',indicatorResidualPatternCoefficientLabels[iIndicator],')','<',indicatorArgs$residualSizeLimitMax[iIndicator],'
                            ')
      }
      
      if(!is.na(indicatorArgs$residualLimitMin[iIndicator])){
        cond.residualSizeLimit=paste0(cond.residualSizeLimit,indicatorResidualPatternCoefficientLabels[iIndicator],">",indicatorArgs$residualLimitMin[iIndicator],"
                                      ")
      }
    }
    cond=paste0(cond,cond.residualSizeLimit)
  }
  
  return(paste0(lds,lds.residuals,lds.factorSelf,lds.factorOther,cond))
  
}

semplate$parseGenomicSEMResult <- function(resultDf = NULL) {
  paths <- resultDf %>%
    select(lhs, op, rhs, Unstand_Est, Unstand_SE, STD_Genotype, STD_Genotype_SE, p_value)
  
  #fix variable types
  paths$Unstand_SE<-as.numeric(paths$Unstand_SE)
  paths$STD_Genotype_SE<-as.numeric(paths$STD_Genotype_SE)
  paths$p_value<-as.numeric(paths$p_value)
  
  # Latent variables: left-hand side of "=~" paths
  var.latent <- paths %>%
    filter(op == "=~") %>%
    select(nodes = lhs) %>%
    distinct %>%
    mutate(type="latent")
  
  # Manifest variables: not latent variables
  `%not_in%` <- Negate(`%in%`)
  var.manifest <- paths %>%
    filter(op != "~1", lhs %not_in% var.latent$nodes) %>%
    select(nodes = lhs) %>%
    distinct %>%
    mutate(type="manifest")
  
  # Residuals: "~~" paths with 
  var.residual <- paths %>%
    filter(op == "~~", lhs==rhs, lhs %in% var.manifest$nodes) %>%
    select(nodes = lhs) %>%
    distinct %>%
    mutate(type="residual")
  
  #all_nodes <- DiagrammeR::combine_edfs(var.latent, var.manifest, var.residual)
  variable <- rbind(var.latent, var.manifest, var.residual) %>%
    mutate(id=row_number())
  
  # Edges, labeled by the factor loading estimates
  edges <- paths %>%
    filter(op != "~1") %>%
    left_join(variable[,c("nodes","id")],by = c("lhs" = "nodes")) %>%
    mutate(from=id) %>%
    select(-id) %>%
    left_join(variable[,c("nodes","id")],by = c("rhs" = "nodes")) %>%
    mutate(to=id) %>%
    select(-id) %>%
    filter(from %not_in% variable$id[which(variable$type=="residual")]) %>%
    filter(to %not_in% variable$id[which(variable$type=="residual")]) %>%
    left_join(variable[which(variable$type=="residual"),c("nodes","id")],by = c("lhs" = "nodes")) %>%
    mutate(tofrom.residual=id) %>%
    select(-id) #%>%
    #select(-Unstand_Est)
  
  #may need more filters here
  # factor loadings
  pattern <- edges %>%
    filter(op == "=~") %>%
    distinct %>%
    mutate(type="pattern")
  
  # Regressions: "~" lines
  regression <- edges %>%
    filter(op == "~") %>%
    distinct %>%
    mutate(type="regression")
  
  # Covariances: ~~ for non manifest variable residuals
  covariance <- edges %>%
    filter(op == "~~") %>%
    filter(is.na(tofrom.residual)) %>%
    distinct %>%
    mutate(type="covariance")
  
  # Residual loadings: ~~ for non manifest variable residuals
  residual <- edges %>%
    filter(op == "~~") %>%
    filter(!is.na(tofrom.residual)) %>%
    distinct %>%
    mutate(type="residual")
  
  # Residual self loadings: ~~ for non manifest variable residuals
  # residual_self_loading <- edges %>%
  #   filter(op == "~~") %>%
  #   filter(!is.na(tofrom.residual))
  
  loading<-rbind(pattern, regression, covariance, residual) %>%
    mutate(id=row_number())
  
  toReturn<-list(variable=variable, loading=loading)
  return(toReturn)
  
}

semplate$parseGenomicSEMResultAsMatrices <- function(resultDf){
  gsemResult<-semplate$parseGenomicSEMResult(resultDf)
  factors<-gsemResult$variable[which(gsemResult$variable$type=="latent"),]
  manifestVariables<-gsemResult$variable[which(gsemResult$variable$type=="manifest"),]
  residualVariables<-gsemResult$variable[which(gsemResult$variable$type=="residual"),]
  
  patternCoefficients<-gsemResult$loading[which(gsemResult$loading$type=="pattern"),]
  residualVaraiances<-gsemResult$loading[which(gsemResult$loading$type=="residual"),]
  
  nFactors<-nrow(factors)
  nManifestVariables<-nrow(manifestVariables)
  patternCoefficients.matrix<-matrix(data = NA_real_, nrow = nManifestVariables, ncol = nFactors)
  
  for(nFactor in 1:nFactors){
    #nFactor<-1
    factor.from<-factors$id[nFactor]
    for(nManifestVariable in 1:nManifestVariables) {
      #nManifestVariable<-1
      manifestVariable.to<-manifestVariables$id[nManifestVariable]
      toadd<-patternCoefficients[which(patternCoefficients$from==factor.from & patternCoefficients$to==manifestVariable.to),c("STD_Genotype")]
      if(!is.null(toadd)&length(toadd)>0) patternCoefficients.matrix[nManifestVariable,nFactor]<-toadd
    }
  }
  rownames(patternCoefficients.matrix)<-manifestVariables$nodes
  colnames(patternCoefficients.matrix)<-factors$nodes
  
  manifestOrder<-row_number(manifestVariables$id[which(manifestVariables$id==residualVaraiances$from)])
  residualVaraiances<-residualVaraiances[manifestOrder,]
  residualVaraiances[,c("residualVariable")]<-residualVariables$nodes[which(residualVariables$id==residualVaraiances$tofrom.residual)]
  residualVaraiances.matrix<-as.matrix(residualVaraiances[,c("STD_Genotype")])
  rownames(residualVaraiances.matrix)<-residualVaraiances$residualVariable
  
  return(list(patternCoefficients=patternCoefficients.matrix, residualVariances=residualVaraiances.matrix))
}


semplate$parseGenomicSEMResultAsDOTDataframes <- function(resultDf = NULL) {
  
  gsemResult<-semplate$parseGenomicSEMResult(resultDf)
  
  # Latent variables: left-hand side of "=~" paths
  latent <- gsemResult$variable[which(gsemResult$variable$type=="latent"),] %>%
    mutate(
           shape = "oval",
           label = nodes,
           fontname = 'Garamond',
           fontsize = '10',
           fixedsize = 'true',
           width = '1.7',
           fillcolor = 'moccasin',
           rank = 'min'
    )
  
  # Manifest variables: not latent variables
  manifest <- gsemResult$variable[which(gsemResult$variable$type=="manifest"),] %>%
    mutate(
           shape = "circle", #since already not original measurements
           label = nodes,
           fontname = 'Garamond',
           fontsize = '10',
           fixedsize = 'true',
           width = '0.8',
           fillcolor = 'GhostWhite',
           rank = '1'
    )
  
  # Residuals: "~~" paths with 
  residual <- gsemResult$variable[which(gsemResult$variable$type=="residual"),]  %>%
    mutate(
      shape = "circle",
     label = paste0("U(",nodes,")"),
     fontname = 'Garamond',
     fontsize = '8',
     fixedsize = 'true',
     width = '0.6',
     fillcolor = 'aliceblue',
     rank = 'max'
    )
  
  # Nodes are prepared
  node_set <- combine_edfs(latent, manifest, residual)
  
  
  # factor loadings
  loadings <- gsemResult$loading[which(gsemResult$variable$type=="pattern"),] %>%
    mutate(type="loading", dir="forward", rel = "leading_to", minlen="5", headport="n", weight="2")
  
  # Regressions: "~" lines
  regressions <- gsemResult$loading[which(gsemResult$variable$type=="regression"),] %>%
    mutate(type="regression", dir="forward", rel = "leading_to", headport="n", weight="2")
  
  # Covariances: ~~ for non manifest variable residuals
  covariances <- gsemResult$loading[which(gsemResult$variable$type=="covariance"),] %>%
    mutate(type="covariance", rel = "related", dir="both")
  
  # Set longer length between latent factors
  covariances <- covariances %>%
    mutate(minlen = if_else( lhs %in% node_set$nodes[which(node_set$type=="latent")],"5","0.5"))
  
  # Set variance port
  covariances[which(covariances$lhs==covariances$rhs),c("headport")]<-"n"
  covariances[which(covariances$lhs==covariances$rhs),c("tailport")]<-"n"
  
  
  # Residual loadings: ~~ for non manifest variable residuals
  residual_loadings <- gsemResult$loading[which(gsemResult$variable$type=="residual"),] %>%
    mutate(type="residual_loading", rel = "related", dir="forward", from=tofrom.residual, label="1", headport="s", tailport="n")
  
  # Residual self loadings: ~~ for non manifest variable residuals
  residual_self_loadings <- gsemResult$loading[which(gsemResult$variable$type=="residual"),] %>%
    mutate(type="residual_loading", rel = "related", dir="both", from=tofrom.residual, to=tofrom.residual, headport="s", tailport="s")
  
  edge_set <- combine_edfs(loadings, regressions, covariances, residual_loadings, residual_self_loadings) %>% mutate(
    values = round(STD_Genotype, 2),
    values.se = round(STD_Genotype_SE, 2),
    label = paste0(values," (",values.se,")")) %>%
    mutate(
      tofrom.residual=id, 
      style = if_else(p_value<0.05, "solid","dashed")
      )
  
  toReturn<-list(nodeDf=node_set, edgeDf=edge_set)
  return(toReturn)
  
}

semplate$generateDOT <- function(nodeDf = NULL, edgeDf = NULL) {
  out<-c("digraph {")
  out<-c(out,"graph [layout = 'dot',
       rankdir = 'TB',
       outputorder = 'nodesfirst',
       bgcolor = 'white',
       splines = 'line',
       ranksep = '0.5',
       nodesep = '0.3']

node [fontname = 'Garamond',
      fontsize = '10',
      shape = 'circle',
      fixedsize = 'true',
      width = '0.5',
      style = 'filled',
      fillcolor = 'aliceblue',
      color = 'gray70',
      fontcolor = 'gray50']

edge [fontname = 'Helvetica',
     fontsize = '8',
     color = 'gray80',
     arrowsize = '0.4']") 
  
  
  minRankNodes<-c()
  midRankNodes<-c()
  maxRankNodes<-c()
  
  # add/process nodes
  if(!is.null(nodeDf )) {
    nNode<-nrow(nodeDf)
    for(iNode in 1 : nNode) {
      cNode<-nodeDf[iNode,]
      nodeString<-paste0("\'",cNode[c("id")],"\'"," [")
      
      if("label" %in% colnames(cNode) && !is.na(cNode[c("label")]))
        nodeString<-paste0(nodeString,"label=","\'",cNode[c("label")],"\'",",")
      
      if("style" %in% colnames(cNode) && !is.na(cNode[c("style")]))
        nodeString<-paste0(nodeString,"style=","\'",cNode[c("style")],"\'",",")
      
      if("shape" %in% colnames(cNode) && !is.na(cNode[c("shape")]))
        nodeString<-paste0(nodeString,"shape=","\'",cNode[c("shape")],"\'",",")
      
      if("fontname" %in% colnames(cNode) && !is.na(cNode[c("fontname")]))
        nodeString<-paste0(nodeString,"fontname=","\'",cNode[c("fontname")],"\'",",")
      
      if("fontsize" %in% colnames(cNode) && !is.na(cNode[c("fontsize")]))
        nodeString<-paste0(nodeString,"fontsize=","\'",cNode[c("fontsize")],"\'",",")
      
      if("fixedsize" %in% colnames(cNode) && !is.na(cNode[c("fixedsize")]))
        nodeString<-paste0(nodeString,"fixedsize=","\'",cNode[c("fixedsize")],"\'",",")
      
      if("width" %in% colnames(cNode) && !is.na(cNode[c("width")]))
        nodeString<-paste0(nodeString,"width=","\'",cNode[c("width")],"\'",",")
      
      if("color" %in% colnames(cNode) && !is.na(cNode[c("color")]))
        nodeString<-paste0(nodeString,"color=","\'",cNode[c("color")],"\'",",")
      
      if("fontcolor" %in% colnames(cNode) && !is.na(cNode[c("fontcolor")]))
        nodeString<-paste0(nodeString,"fontcolor=","\'",cNode[c("fontcolor")],"\'",",")
      
      if("fillcolor" %in% colnames(cNode) && !is.na(cNode[c("fillcolor")]))
        nodeString<-paste0(nodeString,"fillcolor=","\'",cNode[c("fillcolor")],"\'",",")
      
      
      nodeString<-paste0(nodeString,"]")
      
      out<-c(out,nodeString)
      
      
      # store rank settings
      
      if("type" %in% colnames(cNode) && cNode[c("type")]=="latent")
        minRankNodes<-c(minRankNodes,cNode[c("id")])
      
      if("type" %in% colnames(cNode) && cNode[c("type")]=="manifest")
        midRankNodes<-c(midRankNodes,cNode[c("id")])
      
      if("type" %in% colnames(cNode) && cNode[c("type")]=="residual")
        maxRankNodes<-c(maxRankNodes,cNode[c("id")])
      
      
    }
  }
  
  # add edges
  if(!is.null(edgeDf )) {
    nEdge<-nrow(edgeDf)
    for(iEdge in 1 : nEdge) {
      cEdge<-edgeDf[iEdge,]
      edgeString<-paste0("\'",cEdge[c("from")],"\'","->","\'",cEdge[c("to")],"\'"," [")
      
      if("id" %in% colnames(cEdge) && !is.na(cEdge[c("id")]))
        edgeString<-paste0(edgeString,"id=","\'",cEdge[c("id")],"\'",",")
      
      if("dir" %in% colnames(cEdge) && !is.na(cEdge[c("dir")]))
        edgeString<-paste0(edgeString,"dir=","\'",cEdge[c("dir")],"\'",",")
      
      if("weight" %in% colnames(cEdge) && !is.na(cEdge[c("weight")]))
        edgeString<-paste0(edgeString,"weight=","\'",cEdge[c("weight")],"\'",",")
      
      if("label" %in% colnames(cEdge) && !is.na(cEdge[c("label")]))
        edgeString<-paste0(edgeString,"label=","\'",cEdge[c("label")],"\'",",")
      
      if("style" %in% colnames(cEdge) && !is.na(cEdge[c("style")]))
        edgeString<-paste0(edgeString,"style=","\'",cEdge[c("style")],"\'",",")
      
      if("color" %in% colnames(cEdge) && !is.na(cEdge[c("color")]))
        edgeString<-paste0(edgeString,"color=","\'",cEdge[c("color")],"\'",",")
      
      if("len" %in% colnames(cEdge) && !is.na(cEdge[c("len")]))
        edgeString<-paste0(edgeString,"len=","\'",cEdge[c("len")],"\'",",")
      
      if("minlen" %in% colnames(cEdge) && !is.na(cEdge[c("minlen")]))
        edgeString<-paste0(edgeString,"minlen=","\'",cEdge[c("minlen")],"\'",",")
      
      # if("splines" %in% colnames(cEdge) && !is.na(cEdge[c("splines")]))
      #   edgeString<-paste0(edgeString,"splines=","\'",cEdge[c("splines")],"\'",",")
      
      if("fontname" %in% colnames(cEdge) && !is.na(cEdge[c("fontname")]))
        edgeString<-paste0(edgeString,"fontname=","\'",cEdge[c("fontname")],"\'",",")
      
      if("fontsize" %in% colnames(cEdge) && !is.na(cEdge[c("fontsize")]))
        edgeString<-paste0(edgeString,"fontsize=","\'",cEdge[c("fontsize")],"\'",",")
      
      if("penwidth" %in% colnames(cEdge) && !is.na(cEdge[c("penwidth")]))
        edgeString<-paste0(edgeString,"penwidth=","\'",cEdge[c("penwidth")],"\'",",")
      
      if("labelfloat" %in% colnames(cEdge) && !is.na(cEdge[c("labelfloat")]))
        edgeString<-paste0(edgeString,"labelfloat=","\'",cEdge[c("labelfloat")],"\'",",")
      
      if("headport" %in% colnames(cEdge) && !is.na(cEdge[c("headport")]))
        edgeString<-paste0(edgeString,"headport=","\'",cEdge[c("headport")],"\'",",")
      
      if("tailport" %in% colnames(cEdge) && !is.na(cEdge[c("tailport")]))
        edgeString<-paste0(edgeString,"tailport=","\'",cEdge[c("tailport")],"\'",",")
      
      
      edgeString<-paste0(edgeString,"]")
      
      out<-c(out,edgeString)
    }
  }
  
  # add rank settings
  if(length(minRankNodes)>0)
    out<-c(out,paste0("subgraph{rank=min;",paste0(minRankNodes,collapse = ";"),"}"))
  
  # add rank settings
  if(length(midRankNodes)>0)
    out<-c(out,paste0("subgraph{rank=same;",paste0(midRankNodes,collapse = ";"),"}"))
  
  # add rank settings
  if(length(maxRankNodes)>0)
    out<-c(out,paste0("subgraph{rank=max;",paste0(maxRankNodes,collapse = ";"),"}"))
  
  out<-c(out,"}")
  out.string<-paste(out,collapse = "\n")
  return(out.string)
}

semplate$parseAndPrintGenomicSEMResult <- function(resultDf = NULL){
  parsedSEMResult<-semplate$parseGenomicSEMResultAsDOTDataframes(resultDf)
  dot<-semplate$generateDOT(nodeDf=parsedSEMResult$nodeDf, edgeDf=parsedSEMResult$edgeDf)
  grViz(dot)
  return(dot)
}

# Needs and depends on the Genomic SEM package currently loaded - NOT FINISHED!
semplate$factorGWAS<-function(modelDefinitionDf, covstruc, SNPs, estimation="ML", GC){
  modelDefinitionAndResults<-modelDefinitionDf
  nModelDefinition<-nrow(modelDefinitionAndResults)
  for(iModelDefinition in 1:nModelDefinition){
    cModelDefinition<-modelDefinitionAndResults[iModelDefinition,]
    
    modelDefinitionAndResults[iModelDefinition,"indexFactorGWAS"]<-iModelDefinition
    
    #Amend models without SNP effects with SNP effect factors, and construct sub
    cSub<-c()
    for(iLatentFactor in 1:cModelDefinition$nLatentFactor){
      nSub<-paste0("F",iLatentFactor,"~SNP")
      cSub<-c(cSub,nSub)
      cModelDefinition$definition<-paste0(cModelDefinition$definition,"\n",nSub)
    }
    
    modelDefinitionAndResults[iModelDefinition,"definition.full"]<-cModelDefinition$definition
    
    #uncomment this
    #factorGWASResult<-userGWAS(covstruc = covstruc, SNPs = SNPs, estimation = estimation, model = cModelDefinition$definition, modelchi = FALSE, printwarn = TRUE, sub=cSub, GC=GC)
    
    #TEST
    cat("Model ",iModelDefinition,":",cModelDefinition$definition)
    
    #modelDefinitionAndResults[iModelDefinition,"factorGWASResult"]<-list(factorGWASResult)
    
  }
  
  return(modelDefinitionAndResults)
  
}



