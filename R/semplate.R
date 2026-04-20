semplate <-c()
semplate$powers.of.two32bit <- 2^(0:(32 - 1))
semplate$bit.one<-intToBits(1)[1]
semplate$bit.zero<-intToBits(0)[1]

#' @title Generate indicator loading patterns from factor loadings
#' @note Updated from the PhD project version to provide a more conveniently formatted result.
#' 
#' @usage semplate$generateIndicatorLoadingPatternsFromFactorLoadings(factorLoadings,increment,forceOneIndicatorLoading,valueCutoffStart)
#' @param factorLoadings data frame, factors as columns and indicators as rows
#' @return data frame, loading patterns, boolean values per indicator/factor combination
#' @export
semplate$generateIndicatorLoadingPatternsFromFactorLoadings<-function(factorLoadings, label='loadings', increment=0.0005, forceOneIndicatorLoading=T,valueCutoffStart=0.05, labelFactorsInOrder=TRUE){
  #factorLoadings = fa.combined.dimensions.g.5f$loadings
  
  nIndicators<-nrow(factorLoadings)
  nFactors<-ncol(factorLoadings)
  
  indicatorLoadingPatterns<-c()
  factorLoadings<-base::abs(factorLoadings)
  
  if(forceOneIndicatorLoading){
    forced<-as.vector(apply(factorLoadings,MARGIN = 1,FUN = which.max))
  }
  
  cValue<-valueCutoffStart
  iLoop<-1
  while(cValue < max(factorLoadings)){
    
    toadd<-(factorLoadings>cValue)
    #if(!all(toadd) | !any(toadd)){
    if(forceOneIndicatorLoading) {
        for(iIndicator in 1 : nIndicators){
          toadd[iIndicator,forced[iIndicator]]<-TRUE
      }
    }
    
    if(labelFactorsInOrder){
      colnames(toadd)<-1:nFactors
    }
    
    #indicatorLoadingPatterns<-rbind(indicatorLoadingPatterns,as.vector(x = unlist(as.data.frame(toadd)), mode = 'logical')) #unique works across dataframes also!!
    #indicatorLoadingPatterns[[iLoop]]<-as.data.frame(toadd)
    indicatorLoadingPatterns[[iLoop]]<-as.matrix(toadd)
    cValue<-(cValue+increment)
    iLoop<-iLoop+1
  }
  
  indicatorLoadingPatterns<-unique(indicatorLoadingPatterns)
  
  indicatorLoadingPatterns.df<-data.frame(
    label=rep(label,length(indicatorLoadingPatterns)),
    nIndicators=rep(nIndicators,length(indicatorLoadingPatterns)),
    nFactors=rep(nFactors,length(indicatorLoadingPatterns))
    #indicatorLoadingPatterns=indicatorLoadingPatterns
    )
  
  indicatorLoadingPatterns.df$indicatorLoadingPatterns<-NA
  
  for(iDf in 1:nrow(indicatorLoadingPatterns.df)){
    #iDf<-1
    indicatorLoadingPatterns.df[[iDf,c("indicatorLoadingPatterns")]]<-list(indicatorLoadingPatterns[[iDf]])
  }
  
  return(indicatorLoadingPatterns.df)
}

#deprecated
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


#defaults
# fix_loading.table.indicator_factor=NULL
# fix_correlation.table.factor_factor=NULL
# fixResidualVariance_v=NULL
# orthogonal=FALSE
# indicatorArgs=NULL
# universalResidualLimitMin=0.001
# universalIndicatorLoadingAbsoluteLimitMax=1
# universalCorrelationLimitMax=NA

#test
# allow_loading.table.indicator_factor = cIndicatorLoadings
# orthogonal = (cCorrelation=="ORT")
# universalCorrelationLimitMax = ifelse((cCorrelation=="OBL"),0.3,NA)

semplate$generateLavaanCFAModel<-function(
    allow_loading.table.indicator_factor,
    fix_loading.table.indicator_factor=NULL,
    fix_correlation.table.factor_factor=NULL, 
    fixResidualVariance_v=NULL,
    orthogonal=FALSE,
    indicatorArgs=NULL,
    universalResidualLimitMin=0.001,
    universalIndicatorLoadingAbsoluteLimitMax=1,
    universalCorrelationLimitMax=NA
    ){
  #allow_loading.table.indicator_factor<-cIndicatorLoadings
  
  
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
  if(is.null(indicatorArgs$indicatorLoadingAbsoluteLimitMax)){indicatorArgs$indicatorLoadingAbsoluteLimitMax=universalIndicatorLoadingAbsoluteLimitMax}
  
  
  indicatorArgs[which(is.na(indicatorArgs$residualLimitMin)),c("residualLimitMin")]<-universalResidualLimitMin
  
  correlationPatternCoefficientLabels=c()
  
  indicatorResidualPatternCoefficientLabels=c()
  #indicatorResidualPatternCoefficientLabelDefinitions=c()
  
  cond.indicatorLoadingSizeLimit=""
  
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
            vFixLoading=paste0("lf",iFactor,indicatorArgs$code[iIndicator],"*") else
              vFixLoading=paste0("NA*lf",iFactor,indicatorArgs$code[iIndicator],"*")
            
          #loading constraints - uses standardised
          if(!is.na(indicatorArgs$indicatorLoadingAbsoluteLimitMax[iIndicator])){
            cond.indicatorLoadingSizeLimit=paste0(cond.indicatorLoadingSizeLimit,'abs(lf',iFactor,indicatorArgs$code[iIndicator],')','<',indicatorArgs$indicatorLoadingAbsoluteLimitMax[iIndicator],'
                        ')
          }
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
  
  cond=paste0(cond,cond.indicatorLoadingSizeLimit) #add in indicator loading size limit rules (from the previous factor-indicator loops) to general condition string
  
  for(iIndicator in 1:nIndicators){
    if(is.null(fixResidualVariance_v)){
    indicatorResidualPatternCoefficientLabels[iIndicator]=paste0('r',iIndicator)
    #indicatorResidualPatternCoefficientLabelDefinitions[iIndicator]=paste0('label(',indicatorResidualPatternCoefficientLabels[iIndicator],')')
    
    lds.residuals=paste0(lds.residuals,indicatorArgs$code[iIndicator],"~~",indicatorResidualPatternCoefficientLabels[iIndicator],"*",indicatorArgs$code[iIndicator],"
                         ")
    } else {
      lds.residuals=paste0(lds.residuals,indicatorArgs$code[iIndicator],"~~(",fixResidualVariance_v[iIndicator],")*",indicatorArgs$code[iIndicator],"
                         ")
    }
  }
  
  #cross join to get all factor combinations
  indexFactor<-data.frame(index=c(1:nFactors))
  indexFactor<-indexFactor %>%
    dplyr::inner_join(indexFactor, by = character())
  
  
  lds.factorOther=""
  
  iCorrelationPatternCoefficient<-0
  for(iComb in 1:nrow(indexFactor)) {
    
    row<-indexFactor[iComb,]
    if(row$index.x==row$index.y || row$index.y < row$index.x)
      next()
    
    if(orthogonal){
      lds.factorOther=paste0(lds.factorOther,"
                           F",row$index.x,"~~0*F",row$index.y) 
    } else {
      iCorrelationPatternCoefficient<-iCorrelationPatternCoefficient+1
      correlationPatternCoefficientLabels[iCorrelationPatternCoefficient]<-paste0('c',iComb)
      if(is.null(fix_correlation.table.factor_factor)==FALSE){
        #use fixed correlation coefficients
        lds.factorOther=paste0(lds.factorOther,"
                           F",row$index.x,"~~",ifelse(is.na(fix_correlation.table.factor_factor[row$index.x,row$index.y]),"",paste0(fix_correlation.table.factor_factor[row$index.x,row$index.y],"*")),"F",row$index.y)
      } else {
        #do not use fixed correlation coefficients
        lds.factorOther=paste0(lds.factorOther,"
                           F",row$index.x,"~~",ifelse(is.na(universalCorrelationLimitMax),"",paste0(correlationPatternCoefficientLabels[iCorrelationPatternCoefficient],"*")),"F",row$index.y)
      }
    }
  }
  
  cond.correlationSizeLimit=""
  if(length(correlationPatternCoefficientLabels)>0 & !is.na(universalCorrelationLimitMax)){
    for(iCorrelationPatternCoefficient in 1:length(correlationPatternCoefficientLabels)){
      cond.correlationSizeLimit=paste0(cond.correlationSizeLimit,"abs(",correlationPatternCoefficientLabels[iCorrelationPatternCoefficient],")", "<",universalCorrelationLimitMax,"
                                        ")
    }
    cond=paste0(cond,cond.correlationSizeLimit)
  }
  
  
  if(is.null(fixResidualVariance_v)){
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

#Most of the code here was probably stolen from this example here: https://rpubs.com/tjmahr/sem_diagrammer
semplate$parseGenomicSEMResult <- function(resultDf = NULL) {
  #resultDf<-genomicSEMResult
  
  paths <- resultDf %>%
    dplyr::select(lhs, op, rhs, Unstand_Est, Unstand_SE, STD_Genotype, STD_Genotype_SE, p_value)
  
  #fix variable types
  paths$Unstand_SE<-as.numeric(paths$Unstand_SE)
  paths$STD_Genotype_SE<-as.numeric(paths$STD_Genotype_SE)
  paths$p_value<-as.numeric(paths$p_value)
  
  # Latent variables: left-hand side of "=~" paths
  var.latent <- paths %>%
    dplyr::filter(op == "=~") %>%
    dplyr::select(nodes = lhs) %>%
    dplyr::distinct() %>%
    dplyr::mutate(type="latent")
  
  # Manifest variables: not latent variables
  `%not_in%` <- Negate(`%in%`)
  var.manifest <- paths %>%
    dplyr::filter(op != "~1", lhs %not_in% var.latent$nodes) %>%
    dplyr::select(nodes = lhs) %>%
    dplyr::distinct() %>%
    dplyr::mutate(type="manifest")
  
  # Residuals: "~~" paths with 
  var.residual <- paths %>%
    dplyr::filter(op == "~~", lhs==rhs, lhs %in% var.manifest$nodes) %>%
    dplyr::select(nodes = lhs) %>%
    dplyr::distinct() %>%
    dplyr::mutate(type="residual")
  
  #all_nodes <- DiagrammeR::combine_edfs(var.latent, var.manifest, var.residual)
  variable <- rbind(var.latent, var.manifest, var.residual) %>%
    dplyr::mutate(id=dplyr::row_number())
  
  # Edges, labeled by the factor loading estimates
  edges <- paths %>%
    dplyr::filter(op != "~1") %>%
    dplyr::left_join(variable[variable$type!="residual",c("nodes","id")],by = c("lhs" = "nodes")) %>%
    dplyr::mutate(from=id) %>%
    dplyr::select(-id) %>%
    dplyr::left_join(variable[variable$type!="residual",c("nodes","id")],by = c("rhs" = "nodes")) %>%
    dplyr::mutate(to=id) %>%
    dplyr::select(-id) %>%
    dplyr::filter(from %not_in% variable$id[which(variable$type=="residual")]) %>%
    dplyr::filter(to %not_in% variable$id[which(variable$type=="residual")]) %>%
    dplyr::left_join(variable[variable$type=="residual",c("nodes","id")],by = c("lhs" = "nodes")) %>%
    dplyr::mutate(tofrom.residual=id) %>%
    dplyr::select(-id) #%>%
    #select(-Unstand_Est)
  
  #may need more filters here
  # factor loadings
  pattern <- edges %>%
    dplyr::filter(op == "=~") %>%
    dplyr::distinct() %>%
    dplyr::mutate(type="pattern")
  
  # Regressions: "~" lines
  regression <- edges %>%
    dplyr::filter(op == "~") %>%
    dplyr::distinct() %>%
    dplyr::mutate(type="regression")
  
  # Covariances: ~~ for non manifest variable residuals
  covariance <- edges %>%
    dplyr::filter(op == "~~") %>%
    dplyr::filter(is.na(tofrom.residual)) %>%
    dplyr::distinct() %>%
    dplyr::mutate(type="covariance")
  
  # Residual loadings: ~~ for non manifest variable residuals
  residual <- edges %>%
    dplyr::filter(op == "~~") %>%
    dplyr::filter(!is.na(tofrom.residual)) %>%
    dplyr::distinct() %>%
    dplyr::mutate(type="residual")
  
  # Residual self loadings: ~~ for non manifest variable residuals
  # residual_self_loading <- edges %>%
  #   filter(op == "~~") %>%
  #   filter(!is.na(tofrom.residual))
  
  loading<-rbind(pattern, regression, covariance, residual) %>%
    dplyr::mutate(id=dplyr::row_number())
  
  toReturn<-list(variable=variable, loading=loading)
  return(toReturn)
  
}

semplate$parseGenomicSEMResultAsMatrices <- function(resultDf){
  #resultDf<-genomicSEMResult
  
  gsemResult<-semplate$parseGenomicSEMResult(resultDf)
  factors<-gsemResult$variable[which(gsemResult$variable$type=="latent"),]
  manifestVariables<-gsemResult$variable[which(gsemResult$variable$type=="manifest"),]
  residualVariables<-gsemResult$variable[which(gsemResult$variable$type=="residual"),]
  
  patternCoefficients<-gsemResult$loading[which(gsemResult$loading$type=="pattern"),]
  coVariances<-gsemResult$loading[which(gsemResult$loading$type=="covariance"),]
  residualVariances<-gsemResult$loading[which(gsemResult$loading$type=="residual"),]
  
  nFactors<-nrow(factors)
  nManifestVariables<-nrow(manifestVariables)
  patternCoefficients.matrix<-matrix(data = NA_real_, nrow = nManifestVariables, ncol = nFactors)
  patternCoefficients.SE.matrix<-matrix(data = NA_real_, nrow = nManifestVariables, ncol = nFactors)
  patternCoefficientsSTDGenotype.matrix<-matrix(data = NA_real_, nrow = nManifestVariables, ncol = nFactors)
  patternCoefficientsSTDGenotype.SE.matrix<-matrix(data = NA_real_, nrow = nManifestVariables, ncol = nFactors)
  patternCoefficients.p.matrix<-matrix(data = NA_real_, nrow = nManifestVariables, ncol = nFactors)
  covariances.matrix<-matrix(data = NA_real_, nrow = nFactors, ncol = nFactors)
  covariances.SE.matrix<-matrix(data = NA_real_, nrow = nFactors, ncol = nFactors)
  covariancesSTDGenotype.matrix<-matrix(data = NA_real_, nrow = nFactors, ncol = nFactors)
  covariancesSTDGenotype.SE.matrix<-matrix(data = NA_real_, nrow = nFactors, ncol = nFactors)
  covariancesSTDGenotype.p.matrix<-matrix(data = NA_real_, nrow = nFactors, ncol = nFactors)
  
  #read in matrix elements
  for(nFactor in 1:nFactors){
    #nFactor<-1
    factor.from<-factors$id[nFactor]
    for(nFactor2 in 1:nFactors){
      #nFactor2<-1
      factor.to<-factors$id[nFactor2]
      toadd<-coVariances[which(coVariances$from==factor.from & coVariances$to==factor.to),c("Unstand_Est")]
      if(!is.null(toadd)&length(toadd)>0) covariances.matrix[nFactor2,nFactor]<-toadd
      if(!is.null(toadd)&length(toadd)>0) covariances.matrix[nFactor,nFactor2]<-toadd
      
      toadd<-coVariances[which(coVariances$from==factor.from & coVariances$to==factor.to),c("Unstand_SE")]
      if(!is.null(toadd)&length(toadd)>0) covariances.SE.matrix[nFactor2,nFactor]<-toadd
      if(!is.null(toadd)&length(toadd)>0) covariances.SE.matrix[nFactor,nFactor2]<-toadd
      
      toadd<-coVariances[which(coVariances$from==factor.from & coVariances$to==factor.to),c("STD_Genotype")]
      if(!is.null(toadd)&length(toadd)>0) covariancesSTDGenotype.matrix[nFactor2,nFactor]<-toadd
      if(!is.null(toadd)&length(toadd)>0) covariancesSTDGenotype.matrix[nFactor,nFactor2]<-toadd
      
      toadd<-coVariances[which(coVariances$from==factor.from & coVariances$to==factor.to),c("STD_Genotype_SE")]
      if(!is.null(toadd)&length(toadd)>0) covariancesSTDGenotype.SE.matrix[nFactor2,nFactor]<-toadd
      if(!is.null(toadd)&length(toadd)>0) covariancesSTDGenotype.SE.matrix[nFactor,nFactor2]<-toadd
      
      toadd<-coVariances[which(coVariances$from==factor.from & coVariances$to==factor.to),c("p_value")]
      if(!is.null(toadd)&length(toadd)>0) covariancesSTDGenotype.p.matrix[nFactor2,nFactor]<-toadd
      if(!is.null(toadd)&length(toadd)>0) covariancesSTDGenotype.p.matrix[nFactor,nFactor2]<-toadd
    }
    
    for(nManifestVariable in 1:nManifestVariables) {
      #nManifestVariable<-1
      manifestVariable.to<-manifestVariables$id[nManifestVariable]
      
      toadd<-patternCoefficients[which(patternCoefficients$from==factor.from & patternCoefficients$to==manifestVariable.to),c("Unstand_Est")]
      if(!is.null(toadd)&length(toadd)>0) patternCoefficients.matrix[nManifestVariable,nFactor]<-toadd
      
      toadd<-patternCoefficients[which(patternCoefficients$from==factor.from & patternCoefficients$to==manifestVariable.to),c("Unstand_SE")]
      if(!is.null(toadd)&length(toadd)>0) patternCoefficients.SE.matrix[nManifestVariable,nFactor]<-toadd
      
      toadd<-patternCoefficients[which(patternCoefficients$from==factor.from & patternCoefficients$to==manifestVariable.to),c("STD_Genotype")]
      if(!is.null(toadd)&length(toadd)>0) patternCoefficientsSTDGenotype.matrix[nManifestVariable,nFactor]<-toadd
      
      toadd<-patternCoefficients[which(patternCoefficients$from==factor.from & patternCoefficients$to==manifestVariable.to),c("STD_Genotype_SE")]
      if(!is.null(toadd)&length(toadd)>0) patternCoefficientsSTDGenotype.SE.matrix[nManifestVariable,nFactor]<-toadd
      
      toadd<-patternCoefficients[which(patternCoefficients$from==factor.from & patternCoefficients$to==manifestVariable.to),c("p_value")]
      if(!is.null(toadd)&length(toadd)>0) patternCoefficients.p.matrix[nManifestVariable,nFactor]<-toadd
    }
  }
  rownames(covariances.matrix)<-factors$nodes
  colnames(covariances.matrix)<-factors$nodes
  rownames(covariances.SE.matrix)<-factors$nodes
  colnames(covariances.SE.matrix)<-factors$nodes
  rownames(covariancesSTDGenotype.matrix)<-factors$nodes
  colnames(covariancesSTDGenotype.matrix)<-factors$nodes
  rownames(covariancesSTDGenotype.SE.matrix)<-factors$nodes
  colnames(covariancesSTDGenotype.SE.matrix)<-factors$nodes
  rownames(covariancesSTDGenotype.p.matrix)<-factors$nodes
  colnames(covariancesSTDGenotype.p.matrix)<-factors$nodes
  
  rownames(patternCoefficients.matrix)<-manifestVariables$nodes
  colnames(patternCoefficients.matrix)<-factors$nodes
  rownames(patternCoefficients.SE.matrix)<-manifestVariables$nodes
  colnames(patternCoefficients.SE.matrix)<-factors$nodes
  rownames(patternCoefficientsSTDGenotype.matrix)<-manifestVariables$nodes
  colnames(patternCoefficientsSTDGenotype.matrix)<-factors$nodes
  rownames(patternCoefficientsSTDGenotype.SE.matrix)<-manifestVariables$nodes
  colnames(patternCoefficientsSTDGenotype.SE.matrix)<-factors$nodes
  rownames(patternCoefficients.p.matrix)<-manifestVariables$nodes
  colnames(patternCoefficients.p.matrix)<-factors$nodes
  
  manifestOrder<-dplyr::row_number(manifestVariables$id[which(manifestVariables$id==residualVariances$from)])
  residualVariances<-residualVariances[manifestOrder,]
  residualVariances[,c("residualVariable")]<-residualVariables$nodes[which(residualVariables$id==residualVariances$tofrom.residual)]
  residualVariances.matrix<-as.matrix(residualVariances[,c("Unstand_Est")])
  residualVariances.SE.matrix<-as.matrix(residualVariances[,c("Unstand_SE")])
  residualVariancesSTDGenotype.matrix<-as.matrix(residualVariances[,c("STD_Genotype")])
  residualVariancesSTDGenotype.SE.matrix<-as.matrix(residualVariances[,c("STD_Genotype_SE")])
  residualVariances.p.matrix<-as.matrix(residualVariances[,c("p_value")])
  rownames(residualVariances.matrix)<-residualVariances$residualVariable
  rownames(residualVariances.SE.matrix)<-residualVariances$residualVariable
  rownames(residualVariancesSTDGenotype.matrix)<-residualVariances$residualVariable
  rownames(residualVariancesSTDGenotype.SE.matrix)<-residualVariances$residualVariable
  rownames(residualVariances.p.matrix)<-residualVariances$residualVariable
  
  #calculate variance explained by each latent factor and total
  factorVariance<-(patternCoefficients.matrix^2) #For latent factors with variance 1. Otherwise *Vfactor.
  totalIndicatorVariance<-residualVariances.matrix+rowSums(factorVariance, na.rm = T)
  relativeFactorVariance<-factorVariance/rep(totalIndicatorVariance,ncol(factorVariance)) #this only works for orthogonal factor models.
  meanRelativeFactorVariance <- colMeans(relativeFactorVariance, na.rm = T)
  
  
  modelFit<-data.frame()
  modelFit[1,c("totalVarianceExplained")]<-(1-mean(residualVariancesSTDGenotype.matrix))
  
  
  return(list(
    patternCoefficients.matrix=patternCoefficients.matrix,
    patternCoefficients.SE.matrix=patternCoefficients.SE.matrix,
    patternCoefficientsSTDGenotype.matrix=patternCoefficientsSTDGenotype.matrix,
    patternCoefficientsSTDGenotype.SE.matrix=patternCoefficientsSTDGenotype.SE.matrix,
    patternCoefficients.p.matrix=patternCoefficients.p.matrix,
    covariances.matrix=covariances.matrix,
    covariances.SE.matrix=covariances.SE.matrix,
    covariancesSTDGenotype.matrix=covariancesSTDGenotype.matrix,
    covariancesSTDGenotype.SE.matrix=covariancesSTDGenotype.SE.matrix,
    covariancesSTDGenotype.p.matrix=covariancesSTDGenotype.p.matrix,
    residualVariances.matrix=residualVariances.matrix,
    residualVariances.SE.matrix=residualVariances.SE.matrix,
    residualVariancesSTDGenotype.matrix=residualVariancesSTDGenotype.matrix,
    residualVariancesSTDGenotype.SE.matrix=residualVariancesSTDGenotype.SE.matrix,
    residualVariances.p.matrix=residualVariances.p.matrix,
    meanTotalRelativeVarianceExplained=(1-mean(residualVariancesSTDGenotype.matrix)),
    relativeVarianceExplainedPerFactor=relativeFactorVariance,
    meanRelativeVarianceExplainedPerFactor=meanRelativeFactorVariance
    ))
}

#new DOT generating function
#this tries to use lavaan variable names
semplate$parseSEMMatricesAsDOTDataframes <- function(mLambda,mPsi,mTheta, mLambda.se=NULL, mPsi.se=NULL,mTheta.se=NULL, mLambda.p=NULL, mTheta.p=NULL, mPsi.p=NULL) {
  # mLambda<-lavTech(cSem,what = "est",T)$lambda[sortCodes,]
  # mLambda.std<-lavTech(cSem,what = "std",T)$lambda[sortCodes,]
  # mLambda.se<-lavTech(cSem,what = "se",T)$lambda[sortCodes,]
  # mPsi<-lavTech(cSem,what = "est",T)$psi
  # mPsi.se<-lavTech(cSem,what = "se",T)$psi
  # mTheta<-lavTech(cSem,what = "est",T)$theta[sortCodes,sortCodes]
  # mTheta.se<-lavTech(cSem,what = "se",T)$theta[sortCodes,sortCodes]
  
  if(is.null(mLambda.se)) {
    mLambda.se<-mLambda
    mLambda.se<-NA_real_
  }
  
  if(is.null(mPsi.se)) {
    mPsi.se<-mPsi
    mPsi.se<-NA_real_
  }
  
  if(is.null(mTheta.se)) {
    mTheta.se<-mTheta
    mTheta.se<-NA_real_
  }
  
  lThetaDiagonal<-diag(mTheta,names = T)
  lTheta.seDiagonal<-diag(mTheta.se,names = T)
  
  #gsemResult<-semplate$parseGenomicSEMResult(resultDf) #legacy
  
  #nodes
  
  latentRaw <- data.frame(
    nodes=colnames(mLambda),
    type='latent',
    label=colnames(mLambda)
    )
  
  manifestRaw <- data.frame(
    nodes=rownames(mLambda),
    type='manifest',
    label=rownames(mLambda)
  )
  
  residualRaw <- data.frame(
    nodes=rownames(mLambda),
    type='residual',
    label = paste0("U(",rownames(mLambda),")") 
  )
  
  # # Latent variables: left-hand side of "=~" paths
  # latent <- gsemResult$variable[which(gsemResult$variable$type=="latent"),] %>%
  #   mutate(
  #     shape = "oval",
  #     label = nodes,
  #     fontname = 'Garamond',
  #     fontsize = '10',
  #     fixedsize = 'true',
  #     width = '1.7',
  #     fillcolor = 'moccasin',
  #     rank = 'min'
  #   )
  # 
  # # Manifest variables: not latent variables
  # manifest <- gsemResult$variable[which(gsemResult$variable$type=="manifest"),] %>%
  #   mutate(
  #     shape = "circle", #since already not original measurements
  #     label = nodes,
  #     fontname = 'Garamond',
  #     fontsize = '10',
  #     fixedsize = 'true',
  #     width = '0.8',
  #     fillcolor = 'GhostWhite',
  #     rank = '1'
  #   )
  # 
  # # Residuals: "~~" paths with 
  # residual <- gsemResult$variable[which(gsemResult$variable$type=="residual"),]  %>%
  #   mutate(
  #     shape = "circle",
  #     label = paste0("U(",nodes,")"),
  #     fontname = 'Garamond',
  #     fontsize = '8',
  #     fixedsize = 'true',
  #     width = '0.6',
  #     fillcolor = 'aliceblue',
  #     rank = 'max'
  #   )
  
  # Nodes are prepared
  node_set <- combine_edfs(latentRaw, manifestRaw, residualRaw) #assigns id
  
  #edges
  
  mLambdaFrom <- mLambda
  mLambdaFrom[,colnames(mLambda)]<-matrix(node_set[node_set$nodes %in% colnames(mLambda) & node_set$type == 'latent',]$id,ncol = ncol(mLambda), nrow(mLambda),byrow = TRUE)
  
  mLambdaTo <- mLambda
  mLambdaTo[,colnames(mLambda)]<-matrix(node_set[node_set$nodes %in% rownames(mLambda) & node_set$type == 'manifest',]$id,ncol = ncol(mLambda), nrow(mLambda),byrow = FALSE)
  
  #mLambdaCross <- t(mLambda) %*% mLambda
  
  # dtMLambda<-as.data.frame(mLambda)
  # setDT(dtMLambda)
  
  #this does not distinguish between factor loadings or regressions (yet)
  loadingsRaw <- data.frame(
    from=unlist(as.data.frame(mLambdaFrom)),
    to=unlist(as.data.frame(mLambdaTo)),
    type='loading',
    dir="forward",
    rel = "leading_to",
    values = unlist(as.data.frame(mLambda)),
    values.se = unlist(as.data.frame(mLambda.se))
  )
  
  
  # # factor loadings
  # loadings <- gsemResult$loading[gsemResult$loading$type=="pattern",] %>%
  #   mutate(type="loading", dir="forward", rel = "leading_to", minlen="5", headport="n", weight="2")
  # 
  # # Regressions: "~" lines
  # regressions <- gsemResult$loading[which(gsemResult$loading$type=="regression"),] %>%
  #   mutate(type="regression", dir="forward", rel = "leading_to", headport="n", weight="2")
  
  
  mPsiFrom <- mPsi
  mPsiFrom[,colnames(mLambda)]<-matrix(node_set[node_set$nodes %in% colnames(mLambda) & node_set$type == 'latent',]$id,ncol = ncol(mLambda), ncol(mLambda),byrow = TRUE)
  
  mPsiTo <- mPsi
  mPsiTo[,colnames(mLambda)]<-matrix(node_set[node_set$nodes %in% colnames(mLambda) & node_set$type == 'latent',]$id,ncol = ncol(mLambda), ncol(mLambda),byrow = FALSE)
  
  covariancesRaw <- data.frame(
    from=unlist(as.data.frame(mPsiFrom)),
    to=unlist(as.data.frame(mPsiTo)),
    type='covariance',
    dir="both",
    rel = "related",
    values = unlist(as.data.frame(mPsi)),
    values.se = unlist(as.data.frame(mPsi.se))
  )
  
  
  # # Covariances: ~~ for non manifest variable residuals
  # covariances <- gsemResult$loading[which(gsemResult$loading$type=="covariance"),] %>%
  #   mutate(type="covariance", rel = "related", dir="both")
  # 
  
  lThetaDiagonalFrom<-lThetaDiagonal
  lThetaDiagonalFrom<-node_set[node_set$type == 'residual',]$id
  
  lThetaDiagonalTo<-lThetaDiagonal
  lThetaDiagonalTo<-node_set[node_set$type == 'manifest',]$id
  
  residualLoadingsRaw <- data.frame(
    from=lThetaDiagonalFrom,
    to=lThetaDiagonalTo,
    type='residual_loading',
    dir="forward",
    rel = "related",
    values = 1,
    values.se = 0
  )
  
  # # Residual loadings: ~~ for non manifest variable residuals
  # residual_loadings <- gsemResult$loading[which(gsemResult$variable$type=="residual"),] %>%
  #   mutate(type="residual_loading", rel = "related", dir="forward", from=tofrom.residual, label="1", headport="s", tailport="n")
  # 
  
  residualSelfLoadingsRaw <- data.frame(
    from=lThetaDiagonalFrom,
    to=lThetaDiagonalFrom,
    type='residual_self_loading',
    dir="both",
    rel = "related",
    values = lThetaDiagonal,
    values.se = lTheta.seDiagonal
  )
  
  # # Residual self loadings: ~~ for non manifest variable residuals
  # residual_self_loadings <- gsemResult$loading[which(gsemResult$variable$type=="residual"),] %>%
  #   mutate(type="residual_loading", rel = "related", dir="both", from=tofrom.residual, to=tofrom.residual, headport="s", tailport="s")
  
  
  edge_set <- combine_edfs(loadingsRaw, covariancesRaw, residualLoadingsRaw, residualSelfLoadingsRaw)
  
  edge_set$z<-NA_real_
  edge_set$p<-NA_real_
  
  #edge statistics - p-values
  if(!is.null(mLambda.p)){
    mLambda.z<-sign(mLambda)*qnorm(p = mLambda.p/2, lower.tail = FALSE)
    edge_set[edge_set$type=="loading",c("z")]<-unlist(as.data.frame(mLambda.z))
    edge_set[edge_set$type=="loading",c("p")]<-unlist(as.data.frame(mLambda.p))
    
    mPsi.z<-sign(mPsi)*qnorm(p = mPsi.p/2, lower.tail = FALSE)
    edge_set[edge_set$type=="covariance",c("z")]<-unlist(as.data.frame(mPsi.z))
    edge_set[edge_set$type=="covariance",c("p")]<-unlist(as.data.frame(mPsi.p))
    
    mTheta.z<-sign(mTheta)*qnorm(p = mTheta.p/2, lower.tail = FALSE)
    edge_set[edge_set$type=="residual_self_loading",c("z")]<-diag(mTheta.z)
    edge_set[edge_set$type=="residual_self_loading",c("p")]<-diag(mTheta.p)
    
  } else if(!all(is.na(mLambda.se))){ #inferred from values and SE
    edge_set[,c("z")]<-ifelse(edge_set[,c("values.se")]>0,edge_set[,c("values")]/edge_set[,c("values.se")],NA_real_)
    edge_set[,c("p")]<-2*pnorm(abs(edge_set[,c("z")]),lower.tail = FALSE)
  }
  
  #edge labels
  if(all(is.na(edge_set$p))){
    edge_set[,c("label")]<-
      paste0("",round(edge_set[,c("values")],2))
  } else {
    edge_set[,c("label")]<-
      paste0("",round(edge_set[,c("values")],2)," (",round(edge_set[,c("values.se")],2),")")
  }
  

  
  # edge_set <- combine_edfs(loadings, regressions, covariances, residual_loadings, residual_self_loadings) %>% mutate(
  #   values = round(STD_Genotype, 2),
  #   values.se = round(STD_Genotype_SE, 2),
  #   label = paste0(values," (",values.se,")")) %>%
  #   mutate(
  #     tofrom.residual=id, 
  #     style = if_else(p_value<0.05, "solid","dashed")
  #   )
  
  
  
  
  
  # # Set longer length between latent factors
  # covariances <- covariances %>%
  #   mutate(minlen = if_else( lhs %in% node_set$nodes[which(node_set$type=="latent")],"5","0.5"))
  # 
  # # Set variance port
  # covariances[which(covariances$lhs==covariances$rhs),c("headport")]<-"n"
  # covariances[which(covariances$lhs==covariances$rhs),c("tailport")]<-"n"
  
  
  toReturn<-list(nodeDf=node_set, edgeDf=edge_set)
  return(toReturn)
  
}

#may have an error still. output seems off.
semplate$parseGenomicSEMResultAsDOTDataframes <- function(resultDf = NULL) {
  #resultDf<-genomicSEMResult
  
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
  loadings <- gsemResult$loading[gsemResult$loading$type=="pattern",] %>%
    mutate(type="loading", dir="forward", rel = "leading_to", minlen="5", headport="n", weight="2")
  
  # Regressions: "~" lines
  regressions <- gsemResult$loading[which(gsemResult$loading$type=="regression"),] %>%
    mutate(type="regression", dir="forward", rel = "leading_to", headport="n", weight="2")
  
  # Covariances: ~~ for non manifest variable residuals
  covariances <- gsemResult$loading[which(gsemResult$loading$type=="covariance"),] %>%
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

semplate$runPatternGenomicSEM <-function(indicatorLoadingPatternsDf, indexToRun, covstruct.mvLDSC, estimation = "ML", CFIcalc = FALSE, verbal=TRUE){
  # indicatorLoadingPatternsDf<-indicatorLoadingPatternsAll
  # indexToRun<-1
  # covstruct.mvLDSC<-ldsc.1kg.ldscpp.results
  # #covstruct.mvLDSC<-ldsc.debug.results
  
  originalTraitNames<-colnames(covstruct.mvLDSC$S)
  newTraitNames<-paste0('t',originalTraitNames)
  
  colnames(covstruct.mvLDSC$S)<-rownames(covstruct.mvLDSC$S)<-newTraitNames

  
  resultColumnNames<-c("chisq","df","p_chisq","AIC","CFI","SRMR")
  
  recordCode<-paste0(indicatorLoadingPatternsDf[indexToRun,]$label,"_",indicatorLoadingPatternsDf[indexToRun,]$nFactors,"_",indicatorLoadingPatternsDf[indexToRun,]$nIndicators,"_",indexToRun)
  
  
  allow_loading.table.indicator_factor<-indicatorLoadingPatternsDf[indexToRun,]$indicatorLoadingPatterns[[1]][[1]]
  
  originalTraitNamesPatterns<-rownames(allow_loading.table.indicator_factor)
  newTraitNamesPatterns<-paste0('t',originalTraitNamesPatterns)
  rownames(allow_loading.table.indicator_factor)<-newTraitNamesPatterns
  
  lmodelString<- semplate$generateLavaanCFAModel(
      allow_loading.table.indicator_factor = allow_loading.table.indicator_factor,
      #indicatorArgs = p$sumstats.sel[,c("code","residualSizeLimitMax")],
      #universalResidualLimitMin = 0.0001,
      orthogonal = FALSE
      #universalCorrelationLimitMax = 0.3
    )
  # 
  # lmodelString<-"F1 =~ NA*ltNA*tNA+ltAN*tAN+tHO+tEM+tAG+tNEU
  # F2 =~ NA*tAN+tDE+tDI+tPS+tHO
  # F3 =~ NA*tDE+tEM+tEX+tOP
  # F4 =~ NA*tAN+tDE+tDI+tEM+tCO+tEX+tNEU
  # 
  # F1~~F2
  #                          F1~~F3
  #                          F1~~F4
  #                          F2~~F3
  #                          F2~~F4
  #                          F3~~F4
  # 
  #                           F1~~1*F1
  #                           F2~~1*F2
  #                           F3~~1*F3
  #                           F4~~1*F4
  #                           
  #                         abs(ltNA)<1
  #                         abs(ltAN)<1
  #                         
  # 
  #                         tNA~~r20*tNA
  #                         tAN~~r21*tAN
  #                         tDE~~r22*tDE
  #                         tDI~~r23*tDI
  #                         tPS~~r24*tPS
  #                         tHO~~r25*tHO
  #                         tEM~~r26*tEM
  #                         tEX~~r27*tEX
  #                         tAG~~r28*tAG
  #                         tCO~~r29*tCO
  #                         tOP~~r30*tOP
  #                         tNEU~~r31*tNEU
  # 
  #                         r20>0.001
  #                                      r21>0.001
  #                                      r22>0.001
  #                                      r23>0.001
  #                                      r24>0.001
  #                                      r25>0.001
  #                                      r26>0.001
  #                                      r27>0.001
  #                                      r28>0.001
  #                                      r29>0.001
  #                                      r30>0.001
  #                                      r31>0.001
  # 
  #                        " #working
  
  #evaluate lavaan model in GenomicSEM
  if(verbal) cat("\nEvaluating new model:\t",recordCode,"\n")
  cModelResults = tryCatch(
    usermodel(
      covstruc = covstruct.mvLDSC,
      model = lmodelString,
      estimation = estimation,
      #std.lv=TRUE, #unit variance identification
      #fix_resid = TRUE,
      CFIcalc = CFIcalc
    ), error= function(e) e
  )
  
  resultDataFrame<-as.data.frame(matrix(data = NA,nrow = 0,ncol = 0))
  resultDataFrame[1,c('index','code','lmodel','gsem')]<-c(indexToRun,recordCode,lmodelString,NA)
  
  if(!inherits(cModelResults, "try-error") && !is.null(cModelResults$modelfit)){
    if(nrow(cModelResults$modelfit)>0 && any(resultColumnNames %in% colnames(cModelResults$modelfit))) {
      if(verbal) print(cModelResults$modelfit)
      #record results even though not fitting
      resultDataFrame[[1,c('gsem')]]<-list(cModelResults)
      cRescolnames<-intersect(resultColumnNames,colnames(cModelResults$modelfit))
      resultDataFrame[1,cRescolnames]<-cModelResults$modelfit[1,cRescolnames]
      if(!is.na(as.numeric(cModelResults$modelfit$chisq))){
        #This is considered a fitting model
        if(verbal) cat("\nFITTING!")
      }
    } else {
      if(verbal) cat("\nThe model did not yield correct results.")
    }
  } else {
    if(verbal) cat("\nThe model did not converge.")
  }
  
  return(resultDataFrame)
  
}

semplate$generateDOT <- function(nodeDf, edgeDf, filterLoadingSTDValueHideLT=0.4, filterLoadingSTDValueExcludeLT=0.0, filterLoadingPValueLT=0.05, allowExclude=FALSE, labelStyle="headlabel", penwidthScaleFactor=20) {
  
  # nodeDf=parsedSEMResult$nodeDf
  # edgeDf=parsedSEMResult$edgeDf
  
  # edgeDf[edgeDf$type=='loading' & !is.na(edgeDf$values) & is.finite(edgeDf$values) & edgeDf$values<filterLoadingSTDValueExcludeLT, c("type","from","to","label")]<-NA
  # 
  # edgeDf<-edgeDf[!is.na(edgeDf$type) & !is.na(edgeDf$from) & !is.na(edgeDf$to) & !is.na(edgeDf$label),]
  # 
  
  

  #landscape = 'true'
  #splines = 'line',
  #splines = 'spline', #this leads to oom issues
  #outputorder = 'edgesfirst',
  out<-c("digraph {")
  out<-c(out,"graph [layout = 'neato',
       rankdir = 'TB',
       outputorder = 'nodesfirst',
       bgcolor = 'white',
       splines = 'spline',
       ranksep = '3.0',
       nodesep = '0.3']

node [fontname = 'Garamond',
      fontsize = '30',
      shape = 'circle',
      fixedsize = 'true',
      width = '0.5',
      penwidth = '7.0',
      style = 'filled',
      fillcolor = 'aliceblue',
      color = 'gray70',
      fontcolor = 'gray20']

edge [fontname = 'Helvetica',
     fontsize = '30',
     color = 'black',
     minlen = '1.0',
     len = '2.0',
     penwidth = '4.0',
     arrowsize = '0.5']") 
  
  
  setDT(nodeDf)
  setDT(edgeDf)
  
  #cosmetics
  nodeDf[type=='latent',c("shape","width","height"):=list("oval","6.0","2.0")]
  nodeDf[type=='manifest',c("shape","width","height"):=list("square","2.0","2.0")]
  nodeDf[type=='residual',c("shape","width"):=list("circle","1.0")]
  
  #add dummy nodes for covariances
  maxExistingNodeId<-max(nodeDf$id)
  nodeDfCovariance<-edgeDf[type=='covariance' & from!=to,.(type,from,to,values,values.se,p,label)][,id:=.I+eval(maxExistingNodeId)][,nodes:=paste0("cov",id)]
  nodeDf<-rbindlist(list(nodeDf,nodeDfCovariance),use.names = TRUE, fill = TRUE)
  nodeDf[type=='covariance',c("shape","width","fillcolor","penwidth"):=list("circle","1.0","white","0.0")]
  
  
  #node layout - latent factor layout 1 - only works with neato and fdp
  ##source: https://observablehq.com/@magjac/placing-graphviz-nodes-in-fixed-positions
  layoutLF1_lfN<-nrow(nodeDf[type=='latent',])
  layoutLF1_iN<-nrow(nodeDf[type=='manifest',])
  
  layoutLF1_lfWidth<-6
  layoutLF1_lfHeight<-2
  
  layoutLF1_iWidthHeight<-2
  layoutLF1_rWidthHeight<-1
  
  layoutLF1_lfDeltaXRelW<-1.1
  layoutLF1_lfDeltaYRelH<-0.4
  
  layoutLF1_iDeltaXRelW<-1.8
  
  layoutLF1_rfDeltaYRelLFH<-0.2
  
  layoutLF1_iRanksepYRelH<-(-5.0)
  
  layoutLF1_rRanksepYRelH<-(-3.0)
  
  
  nodeDf[type=='latent',latentIndex:=(.I-1)]
  nodeDf[type=='manifest',indicatorIndex:=(.I-1)]
  nodeDf[type=='residual',residualIndex:=(.I-1)]
  
  nodeDf[type=='latent',layoutXPosUnscaled:=(latentIndex*eval(layoutLF1_lfWidth)*eval(layoutLF1_lfDeltaXRelW))]
  nodeDf[type=='latent',layoutYPosUnscaled:=0]
  #nodeDf[type=='latent',layoutYPosUnscaled:=((layoutLF1_lfN-abs(latentIndex-(layoutLF1_lfN-1)/2))*eval(layoutLF1_lfHeight)*eval(layoutLF1_lfDeltaYRelH))] #curved layout
  
  nodeDfInfo<-nodeDf
  colnames(nodeDfInfo)<-paste0(colnames(nodeDfInfo),"_info")
  nodeDf[nodeDfInfo,on=.(from=id_info),layoutXPosFrom:=i.layoutXPosUnscaled_info]
  nodeDf[nodeDfInfo,on=.(to=id_info),layoutXPosTo:=i.layoutXPosUnscaled_info]
  nodeDf[type=='covariance',layoutXPosUnscaled:=(layoutXPosFrom+layoutXPosTo)/2]
  nodeDf[type=='covariance',layoutYPosUnscaled:=abs(layoutXPosFrom-layoutXPosTo)*eval(layoutLF1_rfDeltaYRelLFH)]
  
  nodeDf[type=='manifest',layoutXPosUnscaled:=(indicatorIndex*eval(layoutLF1_iWidthHeight)*eval(layoutLF1_iDeltaXRelW))]
  nodeDf[type=='manifest',layoutYPosUnscaled:=eval(layoutLF1_iWidthHeight)*eval(layoutLF1_iRanksepYRelH)]
  
  nodeDf[type=='residual',layoutXPosUnscaled:=(residualIndex*eval(layoutLF1_iWidthHeight)*eval(layoutLF1_iDeltaXRelW))] #should be same as indicators, in the same order
  nodeDf[type=='residual',layoutYPosUnscaled:=eval(layoutLF1_iWidthHeight)*eval(layoutLF1_iRanksepYRelH) + eval(layoutLF1_rWidthHeight)*eval(layoutLF1_rRanksepYRelH)]
  
  layoutWidth<-max(nodeDf$layoutXPosUnscaled)
  lfWidth<-max(nodeDf[type=='latent',]$layoutXPosUnscaled)
  iWidth<-max(nodeDf[type=='manifest',]$layoutXPosUnscaled)
  
  nodeDf[type=='latent' | type=='covariance',layoutXPos:=layoutXPosUnscaled-eval(lfWidth/2)]
  nodeDf[type=='manifest' | type=='residual',layoutXPos:=layoutXPosUnscaled-eval(iWidth/2)]
  
  nodeDf[!is.na(layoutXPosUnscaled) & !is.na(layoutYPosUnscaled),pos:=paste0("",round(layoutXPos,2),",",round(layoutYPosUnscaled,2),"!")]
  
  
  
  edgeDf[type=='loading',c("headport","tailport","minlen"):=list("n","s","2")]
  edgeDf[type=='covariance',c("headport","tailport","minlen"):=list("n","n","2")]
  edgeDf[type=='covariance' & from==to,c("headport","tailport","minlen"):=list("n","n","1")] #variances
  edgeDf[type=='residual_loading',c("headport","tailport","minlen"):=list("s","n","1")]
  edgeDf[type=='residual_self_loading',c("headport","tailport","minlen"):=list("s","s","1")]
  
  #add dummy edges for covariances
  maxExistingEdgeId<-max(edgeDf$id)
  edgeDfCovariance_left<-nodeDf[type=='covariance',.(values, values.se, p, dummyNodeId=id, dummyNodeFrom=from, dummyNodeTo=to)][,type:="dummy_covariance_left"]
  edgeDfCovariance_right<-nodeDf[type=='covariance',.(values, values.se, p, dummyNodeId=id, dummyNodeFrom=from, dummyNodeTo=to)][,type:="dummy_covariance_right"]
  edgeDfCovariance<-rbindlist(list(edgeDfCovariance_left,edgeDfCovariance_right),use.names = TRUE, fill = TRUE)[,id:=.I+eval(maxExistingEdgeId)]
  
  edgeDf<-rbindlist(list(edgeDf,edgeDfCovariance),use.names = TRUE, fill = TRUE)
  edgeDf[type=='dummy_covariance_left', c("dir","rel","headport","tailport","minlen","from","to"):=list("forward","related","n","w","2",dummyNodeId,ifelse(dummyNodeFrom<dummyNodeTo,dummyNodeFrom,dummyNodeTo))]
  edgeDf[type=='dummy_covariance_right', c("dir","rel","headport","tailport","minlen","from","to"):=list("forward","related","n","e","2",dummyNodeId,ifelse(dummyNodeFrom>dummyNodeTo,dummyNodeFrom,dummyNodeTo))]
  
  edgeDf[(type=='loading' | type=='covariance' | type=='dummy_covariance_left' | type=='dummy_covariance_right') & values>0 & from!=to,c("color"):=list("red")]
  edgeDf[(type=='loading' | type=='covariance' | type=='dummy_covariance_left' | type=='dummy_covariance_right') & values<0 & from!=to,c("color"):=list("blue")]
  
  edgeDf[,significantValue:=p<eval(filterLoadingPValueLT)]
  nodeDf[,significantValue:=p<eval(filterLoadingPValueLT)]
  
  edgeDf[!is.na(label) & significantValue==TRUE,label:=paste0(label,"*")]
  nodeDf[!is.na(label) & significantValue==TRUE,label:=paste0(label,"*")]
  
  if(!any(colnames(edgeDf)=="xlabel")) edgeDf[,c("xlabel"):=list(NA_character_)]
  if(labelStyle=="xlabel") {
    edgeDf[!is.na(label),xlabel:=label][,label:=NA]
  } else if(labelStyle=="headlabel"){
    edgeDf[!is.na(label),c("headlabel","labeldistance"):=list(label,'5.0')][,label:=NA] #legacy labeldistance
    edgeDfLabeldistance<-edgeDf[type=='loading' & !is.na(headlabel) & !is.na(values),]
    setorder(edgeDfLabeldistance, to, -values, na.last = TRUE)
    edgeDfLabeldistance[,labeldistanceIndex:=.I]
    edgeDf[edgeDfLabeldistance, on=c("id"), c("labeldistanceIndex"):=list(i.labeldistanceIndex)]
    
    edgeDf[!is.na(labeldistanceIndex) & labeldistanceIndex %% 2 == 1,c("labeldistance"):=list("15.0")] #even loadings get the larger value
    
    edgeDf[type=='covariance' & from!=to & !is.na(labeldistance),c("labeldistance"):=list("20.0")]
  }
  
  edgeDf[,penwidthScale:=abs(values)*penwidthScaleFactor]
  edgeDf[is.finite(values), c("penwidth"):=list(paste0("",round(penwidthScale,1)))]
  #edgeDf[(type=='loading' | type=='covariance') & from!=to & !is.na(values), c("penwidth"):=list(paste0("",round(penwidthScale,1)))]
 
  #do not use
  #edgeDf[,constraint:='false'] #make nodes not change place due to edges
  
  #set initial included/excluded
  edgeDf[,include:=TRUE]
  if(allowExclude){
    edgeDf[is.na(values) | !is.finite(values) | values<filterLoadingSTDValueExcludeLT, c("include"):=list(FALSE)]
    edgeDf[type=='residual_self_loading',include:=TRUE]
  }
  
  edgeDf[type=='covariance' & from!=to, include:=FALSE] #exclude legacy covariance edges as we are using dummy covariances
  
  #set initial hidden
  edgeDf[,hide:=FALSE]
  nodeDf[,hide:=FALSE]
  edgeDf[is.na(values) | !is.finite(values) | values<filterLoadingSTDValueHideLT & from!=to, c("hide"):=list(TRUE)]
  nodeDf[type=='covariance' & from!=to & (is.na(values) | !is.finite(values) | values<filterLoadingSTDValueHideLT), c("hide"):=list(TRUE)]

  edgeDf[,suppress:=FALSE]
  edgeDf[!is.na(values) & is.finite(values) & significantValue==FALSE, c("suppress"):=list(TRUE)]
  
  
  edgeDfLoadings<-edgeDf[type=='loading',][,absValues:=abs(values)]
  setorder(edgeDfLoadings, -significantValue, -values, na.last = TRUE)
  
  #make sure every factor has one loading
  edgeDfLoadings_maxloadingFactor <- edgeDfLoadings[, .(id = head(.SD, 1)$id, values = head(.SD, 1)$values), by = from]
  edgeDf[id %in% edgeDfLoadings_maxloadingFactor$id, c("include","hide"):=list(TRUE,FALSE)]
  
  
  #make sure every indicator has one loading
  edgeDfLoadings_maxloadingIndicator <- edgeDfLoadings[, .(id = head(.SD, 1)$id, values = head(.SD, 1)$values), by = to]
  edgeDf[id %in% edgeDfLoadings_maxloadingIndicator$id, c("include","hide"):=list(TRUE,FALSE)]
  
  
  
  
  #suppressed style
  edgeDf[suppress==TRUE, c("style"):=list("dashed")]
  
  #hidden style - overrides suppressed style
  #edgeDf[hide==TRUE, c("penwidth","label","xlabel","headlabel","color"):=list("1.0","",NA,"","gray90")]
  edgeDf[hide==TRUE, c("style","label","xlabel","headlabel"):=list("invis",NA,NA,NA)]
  nodeDf[hide==TRUE, c("style","label","xlabel"):=list("invis",NA,NA)]
  
  
  # 
  # edgeDf[is.na(values) | !is.finite(values) | values<filterLoadingSTDValueExcludeLT, c("label","xlabel","color"):=list("",NA,"gray90")]
  
  #setorder(edgeDf, -hide, significantValue, values, na.last = FALSE) #order for rendering
  
  nodeDf<-as.data.frame(nodeDf)
  edgeDf<-as.data.frame(edgeDf)
  
  
  
  
  minRankNodes<-c()
  midRankNodes<-c()
  maxRankNodes<-c()
  
  # add/process nodes
  if(!is.null(nodeDf )) {
    nNode<-nrow(nodeDf)
    for(iNode in 1 : nNode) {
      cNode<-nodeDf[iNode,]
      nodeString<-paste0("\'",cNode[c("id")],"\'"," [")
      
      if("id" %in% colnames(cNode) && !is.na(cNode[c("id")]))
        nodeString<-paste0(nodeString,"id=","\'",cNode[c("id")],"\'",",")
      
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
      
      if("height" %in% colnames(cNode) && !is.na(cNode[c("height")]))
        nodeString<-paste0(nodeString,"height=","\'",cNode[c("height")],"\'",",")
      
      if("color" %in% colnames(cNode) && !is.na(cNode[c("color")]))
        nodeString<-paste0(nodeString,"color=","\'",cNode[c("color")],"\'",",")
      
      if("fontcolor" %in% colnames(cNode) && !is.na(cNode[c("fontcolor")]))
        nodeString<-paste0(nodeString,"fontcolor=","\'",cNode[c("fontcolor")],"\'",",")
      
      if("fillcolor" %in% colnames(cNode) && !is.na(cNode[c("fillcolor")]))
        nodeString<-paste0(nodeString,"fillcolor=","\'",cNode[c("fillcolor")],"\'",",")
      
      if("penwidth" %in% colnames(cNode) && !is.na(cNode[c("penwidth")]))
        nodeString<-paste0(nodeString,"penwidth=","\'",cNode[c("penwidth")],"\'",",")
      
      if("pos" %in% colnames(cNode) && !is.na(cNode[c("pos")]))
        nodeString<-paste0(nodeString,"pos=","\'",cNode[c("pos")],"\'",",")
      
      
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
      
      if(cEdge$include==FALSE) next #skip edge on instruction
      
      edgeString<-paste0("\'",cEdge[c("from")],"\'","->","\'",cEdge[c("to")],"\'"," [")
      
      if("id" %in% colnames(cEdge) && !is.na(cEdge[c("id")]))
        edgeString<-paste0(edgeString,"id=","\'",cEdge[c("id")],"\'",",")
      
      if("dir" %in% colnames(cEdge) && !is.na(cEdge[c("dir")]))
        edgeString<-paste0(edgeString,"dir=","\'",cEdge[c("dir")],"\'",",")
      
      if("weight" %in% colnames(cEdge) && !is.na(cEdge[c("weight")]))
        edgeString<-paste0(edgeString,"weight=","\'",cEdge[c("weight")],"\'",",")
      
      if("label" %in% colnames(cEdge) && !is.na(cEdge[c("label")]))
        edgeString<-paste0(edgeString,"label=","\'",cEdge[c("label")],"\'",",")
      
      if("xlabel" %in% colnames(cEdge) && !is.na(cEdge[c("xlabel")]))
        edgeString<-paste0(edgeString,"xlabel=","\'",cEdge[c("xlabel")],"\'",",")
      
      if("headlabel" %in% colnames(cEdge) && !is.na(cEdge[c("headlabel")]))
        edgeString<-paste0(edgeString,"headlabel=","\'",cEdge[c("headlabel")],"\'",",")
      
      if("labeldistance" %in% colnames(cEdge) && !is.na(cEdge[c("labeldistance")]))
        edgeString<-paste0(edgeString,"labeldistance=","\'",cEdge[c("labeldistance")],"\'",",")
      
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
      
      if("constraint" %in% colnames(cEdge) && !is.na(cEdge[c("constraint")]))
        edgeString<-paste0(edgeString,"constraint=","\'",cEdge[c("constraint")],"\'",",")
      
      
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


#requires
# library(DiagrammeR)
# library(DiagrammeRsvg)
# library(rsvg)
# library(xml2)
semplate$parseAndPrintSEMResult <- function(mLambda,mPsi,mTheta, mLambda.se=NULL,mPsi.se=NULL,mTheta.se=NULL,  mLambda.p=NULL,mPsi.p=NULL,mTheta.p=NULL, dotFilePath=NULL, vectorFilePath=NULL, rasterFilePath=NULL, doPrint=TRUE, newIndicatorSortOrderCodes=NULL){
  # mLambda<-lavTech(cSem,what = "est",T)$lambda[sortCodes,]
  # mLambda.se<-lavTech(cSem,what = "se",T)$lambda[sortCodes,]
  # mPsi<-lavTech(cSem,what = "est",T)$psi
  # mPsi.se<-lavTech(cSem,what = "se",T)$psi
  # mTheta<-lavTech(cSem,what = "est",T)$theta[sortCodes,sortCodes]
  # mTheta.se<-lavTech(cSem,what = "se",T)$theta[sortCodes,sortCodes]
  # dotFilePath<-"generated.dot"
  # vectorFilePath<-"generated.svg"
  # rasterFilePath = "generated.png"
  # newIndicatorSortOrderCodes<-sortCodes
  
  if(!is.null(newIndicatorSortOrderCodes)){
    mLambda<-mLambda[newIndicatorSortOrderCodes,]
    mLambda.se<-mLambda.se[newIndicatorSortOrderCodes,]
    mLambda.p<-mLambda.p[newIndicatorSortOrderCodes,]

    mTheta<-mTheta[newIndicatorSortOrderCodes,newIndicatorSortOrderCodes]
    mTheta.se<-mTheta.se[newIndicatorSortOrderCodes,newIndicatorSortOrderCodes]
    mTheta.p<-mTheta.se[newIndicatorSortOrderCodes,newIndicatorSortOrderCodes]
  }
  
  
  
  parsedSEMResult<-semplate$parseSEMMatricesAsDOTDataframes(mLambda,mPsi,mTheta,mLambda.se = mLambda.se,mPsi.se = mPsi.se, mTheta.se = mTheta.se, mLambda.p=mLambda.p, mPsi.p=mPsi.p, mTheta.p=mTheta.p)
  #parsedSEMResult<-semplate$parseSEMMatricesAsDOTDataframes(mLambda,mPsi,mTheta,mLambda.se = mLambda.se,mPsi.se = mPsi.se, mTheta.se = mTheta.se)
  dot<-semplate$generateDOT(nodeDf=parsedSEMResult$nodeDf, edgeDf=parsedSEMResult$edgeDf)
  
  if(!is.null(dotFilePath)){
    dotFile <- file(dotFilePath,open = "w",encoding = "UTF-8")
    write(dot,file = dotFile)
    close(dotFile)
  }
  #grViz(diagram = "generated.dot")
  #grViz(diagram = "test.dot")
  graphWidget<-grViz(dot)
  
  if(doPrint) graphWidget
  if(!is.null(rasterFilePath)){
    export_svg(graphWidget) %>% charToRaw %>% rsvg_png(rasterFilePath)
  }
  if(!is.null(vectorFilePath)){
    export_svg(graphWidget) %>% read_xml() %>% write_xml(vectorFilePath) #not sure this works correctly
  }
  
  return(dot)
}

#wrapper for a lavaan object
semplate$parseAndPrintSEMResultLavaan <- function(lavaanResultObject, dotFilePath=NULL, vectorFilePath=NULL, rasterFilePath=NULL, doPrint=TRUE, newIndicatorSortOrderCodes=NULL, useSTD=TRUE) {
  #lavaanResultObject<-cSem
  #lavaanResultObject<-tSEMRes
  
  mLambda<-lavTech(lavaanResultObject,what = "est",T)$lambda
  mLambda.std<-lavTech(lavaanResultObject,what = "std",T)$lambda
  mLambda.se<-lavTech(lavaanResultObject,what = "se",T)$lambda
  mLambda.z<-mLambda/mLambda.se
  mLambda.p<-2*pnorm(abs(mLambda.z),lower.tail = FALSE)
  
  mPsi<-lavTech(lavaanResultObject,what = "est",T)$psi
  mPsi.std<-lavTech(lavaanResultObject,what = "std",T)$psi
  mPsi.se<-lavTech(lavaanResultObject,what = "se",T)$psi
  mPsi.z<-mPsi/mPsi.se
  mPsi.p<-2*pnorm(abs(mPsi.z),lower.tail = FALSE)
  
  mTheta<-lavTech(lavaanResultObject,what = "est",T)$theta
  mTheta.std<-lavTech(lavaanResultObject,what = "std",T)$theta
  mTheta.se<-lavTech(lavaanResultObject,what = "se",T)$theta
  mTheta.z<-mTheta/mTheta.se
  mTheta.p<-2*pnorm(abs(mTheta.z),lower.tail = FALSE)
  
  #View(mLambda.std)
  #View(mLambda*ratio.standardisation)
  
  if(useSTD){
    mLambda.se.std<-(mLambda/mLambda.std)*mLambda.se
    mLambda.se.std[!is.finite(mLambda.se.std)]<-0
    mPsi.se.std<-(mPsi/mPsi.std)*mPsi.se
    mPsi.se.std[!is.finite(mPsi.se.std)]<-0
    mTheta.se.std<-(mTheta/mTheta.std)*mTheta.se
    mTheta.se.std[!is.finite(mTheta.se.std)]<-0
    
    #set original to STD versions
    mLambda<-mLambda.std
    mLambda.se<-mLambda.se.std
    mPsi<-mPsi.std
    mPsi.se<-mPsi.se.std
    mTheta<-mTheta.std
    mTheta.se<-mTheta.se.std
  }
  
  return(semplate$parseAndPrintSEMResult(mLambda,mPsi,mTheta, mLambda.se=mLambda.se,mPsi.se=mPsi.se,mTheta.se=mTheta.se,mLambda.p = mLambda.p, mPsi.p = mPsi.p, mTheta.p = mTheta.p, dotFilePath=dotFilePath, vectorFilePath=vectorFilePath, rasterFilePath=rasterFilePath, doPrint=doPrint, newIndicatorSortOrderCodes=newIndicatorSortOrderCodes))
  
}

#wrapper for a parsed Genomic SEM matrix object
semplate$parseAndPrintSEMResultGSEMMatrix <- function(parsedGenomicSEMMatrices, dotFilePath=NULL, vectorFilePath=NULL, rasterFilePath=NULL, doPrint=TRUE, newIndicatorSortOrderCodes=NULL, useSTD=TRUE) {
  #parsedGenomicSEMMatrices<-parsedGenomicSEMResults.best6_12
  
  if(useSTD){
    mLambda<-parsedGenomicSEMMatrices$patternCoefficientsSTDGenotype.matrix
    mLambda.se<-parsedGenomicSEMMatrices$patternCoefficientsSTDGenotype.SE.matrix
    mLambda.p<-parsedGenomicSEMMatrices$patternCoefficients.p.matrix
    mPsi<-parsedGenomicSEMMatrices$covariancesSTDGenotype.matrix
    mPsi.se<-parsedGenomicSEMMatrices$covariancesSTDGenotype.SE.matrix
    mPsi.p<-parsedGenomicSEMMatrices$covariancesSTDGenotype.p.matrix
    mTheta<-matrix(data=0,nrow = nrow(mLambda),ncol = nrow(mLambda))
    diag(mTheta)<-parsedGenomicSEMMatrices$residualVariancesSTDGenotype.matrix[,1]
    rownames(mTheta)<-colnames(mTheta)<-rownames(mLambda)
    mTheta.se<-matrix(data=0,nrow = nrow(mLambda),ncol = nrow(mLambda))
    diag(mTheta.se)<-parsedGenomicSEMMatrices$residualVariancesSTDGenotype.SE.matrix[,1]
    rownames(mTheta.se)<-colnames(mTheta.se)<-rownames(mLambda)
    mTheta.p<-matrix(data=0,nrow = nrow(mLambda),ncol = nrow(mLambda))
    diag(mTheta.p)<-parsedGenomicSEMMatrices$residualVariances.p.matrix
    rownames(mTheta.p)<-colnames(mTheta.p)<-rownames(mLambda)
  } else {
    mLambda<-parsedGenomicSEMMatrices$patternCoefficients.matrix
    mLambda.se<-parsedGenomicSEMMatrices$patternCoefficients.SE.matrix
    mLambda.p<-parsedGenomicSEMMatrices$patternCoefficients.p.matrix
    mPsi<-parsedGenomicSEMMatrices$covariances.matrix
    mPsi.se<-parsedGenomicSEMMatrices$covariances.SE.matrix
    mPsi.p<-parsedGenomicSEMMatrices$covariancesSTDGenotype.p.matrix
    mTheta<-matrix(data=0,nrow = nrow(mLambda),ncol = nrow(mLambda))
    diag(mTheta)<-parsedGenomicSEMMatrices$residualVariances.matrix[,1]
    rownames(mTheta)<-colnames(mTheta)<-rownames(mLambda)
    mTheta.se<-matrix(data=0,nrow = nrow(mLambda),ncol = nrow(mLambda))
    diag(mTheta.se)<-parsedGenomicSEMMatrices$residualVariances.SE.matrix[,1]
    rownames(mTheta.se)<-colnames(mTheta.se)<-rownames(mLambda)
    mTheta.p<-matrix(data=0,nrow = nrow(mLambda),ncol = nrow(mLambda))
    diag(mTheta.p)<-parsedGenomicSEMMatrices$residualVariances.p.matrix
    rownames(mTheta.p)<-colnames(mTheta.p)<-rownames(mLambda)
  }
  
  return(semplate$parseAndPrintSEMResult(mLambda,mPsi,mTheta, mLambda.se=mLambda.se,mPsi.se=mPsi.se,mTheta.se=mTheta.se, mLambda.p=mLambda.p, mPsi.p=mPsi.p, mTheta.p=mTheta.p, dotFilePath=dotFilePath, vectorFilePath=vectorFilePath, rasterFilePath=rasterFilePath, doPrint=doPrint, newIndicatorSortOrderCodes=newIndicatorSortOrderCodes))
  
}

#broken due to broken parseGenomicSEMResultAsDOTDataframes
semplate$parseAndPrintGenomicSEMResult <- function(resultDf = NULL){
  #resultDf<-genomicSEMResult
  parsedSEMResult<-semplate$parseGenomicSEMResultAsDOTDataframes(resultDf)
  dot<-semplate$generateDOT(nodeDf=parsedSEMResult$nodeDf, edgeDf=parsedSEMResult$edgeDf)
  grViz(dot)
  return(dot)
}

semplate$reliability.H <- function(factorLoadingsVector, singleVariableValues=F){
  if(singleVariableValues) return(sqrt(1/(1+var(factorLoadingsVector,na.rm = T))))
  factorLoadings2<-factorLoadingsVector^2
  factorLoading2Sum<-sum(factorLoadings2/(1-factorLoadings2))
  1/(1+1/factorLoading2Sum)
}

#amend Genomic SEM LDSC results with LDSC++ standard objects
semplate$amendGsemLDSC <- function(gsemLDSC, nBlocks=200){
  
  #test
  #gsemLDSC<-p$mvLD$covstruct.GSEMmvLDSC.1kg
  
  #fix naming
  rownames(gsemLDSC$S)<-colnames(gsemLDSC$S)
  
  #S.E. of S
  r<-nrow(gsemLDSC$S)
  gsemLDSC$S.SE<-matrix(1, r, r)
  gsemLDSC$S.SE[lower.tri(gsemLDSC$S.SE,diag=TRUE)] <-sqrt(diag(gsemLDSC$V))
  gsemLDSC$S.SE[upper.tri(gsemLDSC$S.SE)]<-t(gsemLDSC$S.SE)[upper.tri(gsemLDSC$S.SE)]
  colnames(gsemLDSC$S.SE)<-colnames(gsemLDSC$S)
  rownames(gsemLDSC$S.SE)<-colnames(gsemLDSC$S)
  gsemLDSC$S.SE.unsigned<-gsemLDSC$S.SE
  
  #these end up on the liability scale!
  gsemLDSC$cov.p.liab<-matrix(1, r, r)
  gsemLDSC$cov.p.liab[lower.tri(gsemLDSC$cov.p.liab,diag=TRUE)]<-2 * pnorm(abs(gsemLDSC$S[lower.tri(gsemLDSC$S,diag=TRUE)] / gsemLDSC$S.SE[lower.tri(gsemLDSC$S.SE,diag=TRUE)]), lower.tail = FALSE)
  gsemLDSC$cov.p.liab[upper.tri(gsemLDSC$cov.p.liab)]<-t(gsemLDSC$cov.p.liab)[upper.tri(gsemLDSC$cov.p.liab)]
  gsemLDSC$cov.p.liab[is.na(gsemLDSC$cov.p.liab)]<-1 #fallback for NA
  colnames(gsemLDSC$cov.p.liab)<-colnames(gsemLDSC$S)
  rownames(gsemLDSC$cov.p.liab)<-colnames(gsemLDSC$S)
  
  gsemLDSC$cov.blocks<-matrix(nBlocks, r, r)
  colnames(gsemLDSC$cov.blocks)<-colnames(gsemLDSC$S)
  rownames(gsemLDSC$cov.blocks)<-colnames(gsemLDSC$S)
  
  #enter SEs from diagonal of standardized V
  gsemLDSC$S_Stand.SE<-matrix(1, r, r)
  gsemLDSC$S_Stand.SE[lower.tri(gsemLDSC$S_Stand.SE,diag=TRUE)] <- gsemLDSC$S_Stand.SE[upper.tri(gsemLDSC$S_Stand.SE,diag=TRUE)] <- sqrt(diag(gsemLDSC$V_Stand))
  colnames(gsemLDSC$S_Stand.SE) <- colnames(gsemLDSC$S)
  rownames(gsemLDSC$S_Stand.SE) <- colnames(gsemLDSC$S)
  
  gsemLDSC$S_Stand.SE.unsigned <- gsemLDSC$S_Stand.SE
  
  return(gsemLDSC)
  
}

