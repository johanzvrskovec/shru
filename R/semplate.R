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
  
  indicatorArgs[which(is.na(indicatorArgs$residualLimitMin)),c("residualLimitMin")]<-universalResidualLimitMin
  
  correlationPatternCoefficientLabels=c()
  
  indicatorResidualPatternCoefficientLabels=c()
  #indicatorResidualPatternCoefficientLabelDefinitions=c()
  
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
    dplyr::left_join(variable[,c("nodes","id")],by = c("lhs" = "nodes")) %>%
    dplyr::mutate(from=id) %>%
    dplyr::select(-id) %>%
    dplyr::left_join(variable[,c("nodes","id")],by = c("rhs" = "nodes")) %>%
    dplyr::mutate(to=id) %>%
    dplyr::select(-id) %>%
    dplyr::filter(from %not_in% variable$id[which(variable$type=="residual")]) %>%
    dplyr::filter(to %not_in% variable$id[which(variable$type=="residual")]) %>%
    dplyr::left_join(variable[which(variable$type=="residual"),c("nodes","id")],by = c("lhs" = "nodes")) %>%
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

#test
#resultDf<-p$CFA$models.selected[iSelected,]$gsemResults[[1]][[1]]$results

semplate$parseGenomicSEMResultAsMatrices <- function(resultDf){
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

semplate$runPatternGenomicSEM <-function(indicatorLoadingPatternsDf, indexToRun, covstruct.mvLDSC, estimation = "ML", CFIcalc = FALSE, verbal=TRUE){
  #indicatorLoadingPatternsDf<-indicatorLoadingPatternsAll
  #indexToRun<-31
  #covstruct.mvLDSC<-ldsc.1kg.ldscpp.results
  #covstruct.mvLDSC<-ldsc.debug.results
  
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
  
  # lmodelString<-"                      F1 =~ NA*ADHD06+ANOR02+ANXI03+AUTI09+BIPO03+BORD02+DEPR05+NEUR05+OPEN02+PTSD11+SCHI06+HO+AG+OP+NEU
  #                     F2 =~ NA*ADHD06+AUTI09+BODY14+COAD03+EDUC03+LONG04+PTSD04+PTSD11+DE+DI+PS+HO+EX
  #                     F3 =~ NA*AGAB02+NEUR05+NA+AN+DE+DI+HO+EM+AG+CO+OP+NEU
  #                     F4 =~ NA*AUTI09+BIPO03+COAD03+CONS02+DEPR05+EDUC03+EXTR02+NEUR05+SCHI06+AN+DE+DI+EM+EX+NEU
  #                     F5 =~ NA*ANOR02+DE+DI+PS+HO+EM+CO+OP
  # ADHD06~~r1*ADHD06
  #                        AGAB02~~r2*AGAB02
  #                        ANOR02~~r3*ANOR02
  #                        ANXI03~~r4*ANXI03
  #                        AUTI09~~r5*AUTI09
  #                        BIPO03~~r6*BIPO03
  #                        BODY14~~r7*BODY14
  #                        BORD02~~r8*BORD02
  #                        COAD03~~r9*COAD03
  #                        CONS02~~r10*CONS02
  #                        DEPR05~~r11*DEPR05
  #                        EDUC03~~r12*EDUC03
  #                        EXTR02~~r13*EXTR02
  #                        LONG04~~r14*LONG04
  #                        NEUR05~~r15*NEUR05
  #                        OPEN02~~r16*OPEN02
  #                        PTSD04~~r17*PTSD04
  #                        PTSD11~~r18*PTSD11
  #                        SCHI06~~r19*SCHI06
  #                        NA~~r20*NA
  #                        AN~~r21*AN
  #                        DE~~r22*DE
  #                        DI~~r23*DI
  #                        PS~~r24*PS
  #                        HO~~r25*HO
  #                        EM~~r26*EM
  #                        EX~~r27*EX
  #                        AG~~r28*AG
  #                        CO~~r29*CO
  #                        OP~~r30*OP
  #                        NEU~~r31*NEU
  #                        
  #                           F1~~1*F1
  #                           F2~~1*F2
  #                           F3~~1*F3
  #                           F4~~1*F4
  #                           F5~~1*F5
  #                          F1~~F2
  #                          F1~~F3
  #                          F1~~F4
  #                          F1~~F5
  #                          F2~~F3
  #                          F2~~F4
  #                          F2~~F5
  #                          F3~~F4
  #                          F3~~F5
  #                          F4~~F5
  # r1>0.001
  #                                     r2>0.001
  #                                     r3>0.001
  #                                     r4>0.001
  #                                     r5>0.001
  #                                     r6>0.001
  #                                     r7>0.001
  #                                     r8>0.001
  #                                     r9>0.001
  #                                     r10>0.001
  #                                     r11>0.001
  #                                     r12>0.001
  #                                     r13>0.001
  #                                     r14>0.001
  #                                     r15>0.001
  #                                     r16>0.001
  #                                     r17>0.001
  #                                     r18>0.001
  #                                     r19>0.001
  #                                     r20>0.001
  #                                     r21>0.001
  #                                     r22>0.001
  #                                     r23>0.001
  #                                     r24>0.001
  #                                     r25>0.001
  #                                     r26>0.001
  #                                     r27>0.001
  #                                     r28>0.001
  #                                     r29>0.001
  #                                     r30>0.001
  #                                     r31>0.001"
  # 
  # lmodelString<-"F1 =~ NA*NA+AN+HO+EM+AG+NEU
  # F2 =~ NA*AN+DE+DI+PS+HO
  # F3 =~ NA*DE+EM+EX+OP
  # F4 =~ NA*AN+DE+DI+EM+EX+NEU
  # F5 =~ NA*AN+DE+DI+EM+CO
  # 
  # F1~~c1*F2
  #                          F1~~c2*F3
  #                          F1~~c3*F4
  #                          F1~~c4*F5
  #                          F2~~c5*F3
  #                          F2~~c6*F4
  #                          F2~~c7*F5
  #                          F3~~c8*F4
  #                          F3~~c9*F5
  #                          F4~~c10*F5
  #                     
  #                           F1~~1*F1
  #                           F2~~1*F2
  #                           F3~~1*F3
  #                           F4~~1*F4
  #                           F5~~1*F5
  #                           
  #                           
  #                           c1<=1
  #                           c2<=1
  #                           c3<=1
  #                           c4<=1
  #                           c5<=1
  #                           c6<=1
  #                           c7<=1
  #                           c8<=1
  #                           c9<=1
  #                           c10<=1
  #                        "
  # 
  # lmodelString<-"F1 =~ NA*tNA+tAN+tHO+tEM+tAG+tNEU
  # F2 =~ NA*tAN+tDE+tDI+tPS+tHO
  # F3 =~ NA*tDE+tEM+tEX+tOP
  # F4 =~ NA*tAN+tDE+tDI+tEM+tEX+tNEU
  # F5 =~ NA*tAN+tDE+tDI+tEM+tCO
  # 
  # F1~~c1*F2
  #                          F1~~c2*F3
  #                          F1~~c3*F4
  #                          F1~~c4*F5
  #                          F2~~c5*F3
  #                          F2~~c6*F4
  #                          F2~~c7*F5
  #                          F3~~c8*F4
  #                          F3~~c9*F5
  #                          F4~~c10*F5
  #                     
  #                           F1~~1*F1
  #                           F2~~1*F2
  #                           F3~~1*F3
  #                           F4~~1*F4
  #                           F5~~1*F5
  #                           
  #                           
  #                           c1<1
  #                           c2<1
  #                           c3<1
  #                           c4<1
  #                           c5<1
  #                           c6<1
  #                           c7<1
  #                           c8<1
  #                           c9<1
  #                           c10<1
  #                        "
  # 
  # 
  # lmodelString<-"F1 =~ NA*tNA+tAN+tHO+tEM+tAG+tNEU
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
  # 
  #                        " #working
  # 
  # 
  # lmodelString<-"F1 =~ NA*tNA+tAN+tHO+tEM+tAG+tNEU
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
  resultDataFrame[1,c('index','code','lmodel')]<-c(indexToRun,recordCode,lmodelString)
  
  if(!inherits(cModelResults, "try-error") && !is.null(cModelResults$modelfit)){
    if(nrow(cModelResults$modelfit)>0 && any(resultColumnNames %in% colnames(cModelResults$modelfit))) {
      if(verbal) print(cModelResults$modelfit)
      #record results even though not fitting
      resultDataFrame[1,c('gsem')]<-list(cModelResults)
      cRescolnames<-intersect(resultColumnNames,colnames(cModelResults$modelfit))
      resultDataFrame[1,cRescolnames]<-cModelResults$modelfit[1,cRescolnames]
      if(is.numeric(cModelResults$modelfit$chisq)){
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

