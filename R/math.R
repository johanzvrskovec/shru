#This is a tweaked version of the p.adjust function from the stats package, to which all the credit belongs. The change involves being able to use n which is not fixed to the length of the p vector to allow for setting this based on a PCA analysis for example.
# https://stackoverflow.com/questions/30108510/p-adjust-with-n-than-number-of-tests
# theory: https://stats.libretexts.org/Bookshelves/Applied_Statistics/Biological_Statistics_(McDonald)/06%3A_Multiple_Tests/6.01%3A_Multiple_Comparisons#:~:text=Sometimes%20you%20will%20see%20a,rank%20of%20the%20smallest%20P

# p.mat<-mvldsc$cov.p.liab[lower.tri(mvldsc$cov.p.liab,diag = T)]
# method = "fdr"
# n = effectiveNumberOfTests
# isSymmetric=F

p.adjust2 <- function (p.mat, method = p.adjust.methods, n = NULL, isSymmetric=F) 
{
  
  #p.mat <- corrsForTable.p
  #method = "fdr"
  #n = NULL
  #isSymmetric=F
  
  #+++JZ: this is an added part for dealing with matrices and symmetry
  if(is.matrix(p.mat) & isSymmetric==T){
    p<-p.mat[lower.tri(x = p.mat,diag = T)]
    if(is.null(n)) n<- length(p)
  } else if(is.null(n)){
    p<-p.mat
    n<-length(p.mat)
  } else {
    p<-p.mat
  }
  
  method <- match.arg(method)
  if (method == "fdr") 
    method <- "BH"
  nm <- names(p)
  p <- as.numeric(p)
  p0 <- setNames(p, nm)
  if (all(nna <- !is.na(p))) 
    nna <- TRUE
  else p <- p[nna] #+++JZ: this is making the list effectively shorter in case of NA's
  lp <- length(p)
  #stopifnot(n >= lp) #+++JZ: Now we can run with n < lp !!!! Thank you Stack Overflow!
  if (n <= 1) 
    return(p0)
  if (n == 2 && method == "hommel") 
    method <- "hochberg"
  p0[nna] <- switch(method, bonferroni = pmin(1, n * p), holm = {
    i <- seq_len(lp)
    o <- order(p)
    ro <- order(o)
    pmin(1, cummax((n + 1L - i) * p[o]))[ro]
  }, hommel = {
    if (n > lp) p <- c(p, rep.int(1, n - lp))
    i <- seq_len(n)
    o <- order(p)
    p <- p[o]
    ro <- order(o)
    q <- pa <- rep.int(min(n * p/i), n)
    for (j in (n - 1L):2L) {
      ij <- seq_len(n - j + 1L)
      i2 <- (n - j + 2L):n
      q1 <- min(j * p[i2]/(2L:j))
      q[ij] <- pmin(j * p[ij], q1)
      q[i2] <- q[n - j + 1L]
      pa <- pmax(pa, q)
    }
    pmax(pa, p)[if (lp < n) ro[1L:lp] else ro]
  }, hochberg = {
    i <- lp:1L
    o <- order(p, decreasing = TRUE)
    ro <- order(o)
    pmin(1, cummin((n + 1L - i) * p[o]))[ro]
  }, BH = {
    i <- lp:1L
    o <- order(p, decreasing = TRUE)
    ro <- order(o)
    pmin(1, cummin(n/i * p[o]))[ro]
  }, BY = {
    i <- lp:1L
    o <- order(p, decreasing = TRUE)
    ro <- order(o)
    q <- sum(1/(1L:n))
    pmin(1, cummin(q * n/i * p[o]))[ro]
  }, none = p)
  
  p0 <- ifelse(p0 < p, p, p0) #+++JZ: Added this to avoid any p0 < p when using small n
  
  #+++JZ: this is an added part for dealing with matrices and symmetry
  if(is.matrix(p.mat)){
    if(isSymmetric==T){
      p0.mat <- matrix(NA,nrow = nrow(p.mat),ncol = ncol(p.mat))
      p0.mat[lower.tri(p0.mat,diag = T)]<-p0
      p0.mat[upper.tri(p0.mat,diag = F)]<-t(p0.mat)[upper.tri(p0.mat,diag = F)]
    } else {
      p0.mat <- matrix(p0,nrow = nrow(p.mat),ncol = ncol(p.mat))
    }
    colnames(p0.mat)<-colnames(p.mat)
    rownames(p0.mat)<-rownames(p.mat)
    
    p0.mat
    
  } else {
    p0
  }
}

getEffectiveNumberOfTests <- function(
    covarianceMatrix,
    symmetric = T,
    eigenSumLimit = 0.995
){
  pcaRes <- eigen(covarianceMatrix,symmetric = symmetric)
  eigenSum<-sum(abs(pcaRes$values))
  for(iEV in 1:length(pcaRes$values)){
    if(sum(abs(pcaRes$values[1:iEV]))/eigenSum > eigenSumLimit) break
  }
  
  return(iEV)
}

#m<-c(2,6,10)
symMatrixConform <- function(m){
  if(isSymmetric.matrix(m)) {
    return(m)
  } else if(is.data.frame(m)) {
    if(nrow(m)==ncol(m)){
      return(as.matrix(m))
    } else stop("Non-symmetic dataframe detected")
  } else {
    toReturn <- matrix(data = NA,ncol = length(m), nrow = length(m))
    diag(toReturn)<-m
    return(toReturn)
  }
}

#this test assumes per default that the variables are closely related/ NOT INDEPENDENT and have correlation ~ 1

#These tests have to include both covariance values and their standard errors as the covariances are used for testing the difference in standard errors.

# #Run examnple
# comparisonTest <- difftest.matrix(
#   mValues1 = mvldsc$S, #estimate matrix 1
#   mStandard_errors1 = mvldsc$S.SE, #standard error matrix 1
#   df1 = mvldsc$cov.blocks, #degrees of freedom matrix 1 (the number of samples)
#   mValues2 = mvldscComparison$S, #estimate matrix 2
#   mStandard_errors2 = mvldscComparison$S.SE, #standard error matrix 2
#   df2 = mvldscComparison$cov.blocks, #degrees of freedom matrix 2 (the number of samples)
#   symmetric = T, #are the matrices symmetric? (usually yes if they are cov/corr matrices)
#   mValueCovariances=mvldsc$V_Stand #calculate effective number of tests from a sampling covariance matrix
#   #OR
#   #effectiveNumberOfTests = effectiveNumberOfTests #force an effective number of tests for multiple testing correction
#   )

# t1<-p$mvLD$covstruct.mvLDSC.GSEMemulation.1kg.maf0_01
# t2<-semplate$amendGsemLDSC(  p$mvLD$covstruct.GSEMmvLDSC.1kg)
# 
# mValues1 = t1$S
# mStandard_errors1 = t1$S.SE
# mValues2 = t2$S
# mStandard_errors2 = t2$S.SE
# df1 = t1$cov.blocks
# df2 = t2$cov.blocks
# mValueCovariances <- t1$V_Stand
# mCorrelationEstimate.values = 0.95
# mCorrelationEstimate.standard_errors = 0.95
# effectiveNumberOfTests=NULL
# #mValueCovariances=NULL
# symmetric = T
# eigenSumLimit = 0.995
# #mCorrelationEstimate.values=1

# testRes <- difftest.matrix(
#   mValues1 = p$mvLD$covstruct.mvLDSC.1kg.vbcs.varblock.winfo.altcw$S,
#   mStandard_errors1 = p$mvLD$covstruct.mvLDSC.1kg.vbcs.varblock.winfo.altcw$S.SE,
#   #mStandard_errors1.std = p$mvLD$covstruct.mvLDSC.1kg.vbcs.varblock.winfo.altcw$S.SE.std,
#   df1 = p$mvLD$covstruct.mvLDSC.1kg.vbcs.varblock.winfo.altcw$cov.blocks,
#   mValues2 = p$mvLD$covstruct.mvLDSC.1kg.vbcs.varblock$S,
#   mStandard_errors2 = p$mvLD$covstruct.mvLDSC.1kg.vbcs.varblock$S.SE,
#   #mStandard_errors2.std = p$mvLD$covstruct.mvLDSC.1kg.vbcs.varblock$S.SE.std,
#   df2 = p$mvLD$covstruct.mvLDSC.1kg.vbcs.varblock$cov.blocks,
#   mValueCovariances = p$mvLD$covstruct.mvLDSC.1kg.vbcs.varblock.winfo.altcw$V_Stand
#                 )


difftest.matrix <- function(
    mValues1,
    mStandard_errors1,
    #mStandard_errors1.std=NULL,
    mValues2,
    mStandard_errors2,
    #mStandard_errors2.std=NULL,
    df1 = NULL,
    df2 = NULL,
    #mN1,
    #mN2,
    mCorrelationEstimate.values = 0.95, #assume correlation of 0.95 for more conservative test
    mCorrelationEstimate.standard_errors = 0.95, #assume correlation of 0.95 for more conservative test
    effectiveNumberOfTests=NULL,
    mValueCovariances=NULL,
    symmetric = T,
    eigenSumLimit = 0.995
){
  
  #check if values are matrix - convert to matrix with values on the diagonal otherwise
  mValues1 <- symMatrixConform(mValues1)
  mStandard_errors1 <- symMatrixConform(mStandard_errors1)
  # if(!is.null(mStandard_errors1.std)){
  #   mStandard_errors1.std <- symMatrixConform(mStandard_errors1.std)
  # }
  
  mValues2 <- symMatrixConform(mValues2)
  mStandard_errors2 <- symMatrixConform(mStandard_errors2)
  # if(!is.null(mStandard_errors2.std)){
  #   mStandard_errors2.std <- symMatrixConform(mStandard_errors2.std)
  # }
  
  
  if(is.null(effectiveNumberOfTests) & !is.null(mValueCovariances)){
    effectiveNumberOfTests<-getEffectiveNumberOfTests(covarianceMatrix = mValueCovariances,symmetric = symmetric, eigenSumLimit = eigenSumLimit)
  }
  
  #old function call
  
  # if(!symmetric){
  #   t<-test.meta.deltacovariance(
  #     values1 = mValues1,standard_errors1 = mStandard_errors1,values2 = mValues2,standard_errors2 = mStandard_errors2,n1 = mN1,n2 = mN2,effectiveNumberOfTests = effectiveNumberOfTests,fullyDependent = fullyDependent)
  # } else {
  #   t<-test.meta.deltacovariance(
  #     values1 = mValues1.lm,standard_errors1 = standard_errors1.lm,values2 = mValues2.lm,standard_errors2 = standard_errors2.lm,n1 = mN1.lm,n2 = mN2.lm,effectiveNumberOfTests = effectiveNumberOfTests,fullyDependent = fullyDependent)
  
  #let's not do this yet to keep proper track of the standardisation required for the genetic covariance cutoff in the test of differences in genetic covariance standard error
  # mValues1.lm <- mValues1[lower.tri(mValues1,diag = T)]
  # standard_errors1.lm <- mStandard_errors1[lower.tri(mStandard_errors1,diag = T)]
  # mN1.lm <- mN1[lower.tri(mN1,diag = T)]
  # mValues2.lm <- mValues2[lower.tri(mValues2,diag = T)]
  # standard_errors2.lm <- mStandard_errors2[lower.tri(mStandard_errors2,diag = T)]
  # mN2.lm <- mN2[lower.tri(mN2,diag = T)]
  
  #begin perform test
    
  #test covariance difference
  mTest.values<-abs(mValues1-mValues2)
  
  sde<-sqrt((mStandard_errors1^2)+(mStandard_errors2^2) - 2*mCorrelationEstimate.values*mStandard_errors1*mStandard_errors2) #confirmed to correspond to the formula from the supplement
  
  #the test is double sided overall though, so 2*p!!!!!!!
  
  pTest.values <- 2*pnorm(
                       q = mTest.values,
                       sd = sde,
                       lower.tail = F
                     )
    
  pTest.values.adj<-p.adjust2(pTest.values, method = "fdr", n = effectiveNumberOfTests, isSymmetric = symmetric)
  rownames(pTest.values.adj)<-rownames(pTest.values)
  colnames(pTest.values.adj)<-colnames(pTest.values)
  
  pTest.standard_errors<-NULL
  pTest.standard_errors.adj<-NULL
  
  if(!is.null(mStandard_errors1) & !is.null(mStandard_errors2) & !is.null(df1)){
    
    if(is.null(df2)) df2<-df1
    #this has now been corrected to correspond to the derived formula from the supplement
    
    tVar <- (mStandard_errors1^2)-(mStandard_errors2^2)
    vare <- df1*(mStandard_errors1^2) + df2*(mStandard_errors2^2) - 2*sqrt(2*df1*(mStandard_errors1^2))*sqrt(2*df2*(mStandard_errors2^2))*mCorrelationEstimate.standard_errors #this is simplified and approximates df-1~df to avoid weird behaviour at large df

    #old
    # var_bar <-(mStandard_errors1^2+mStandard_errors2^2)/2
    # vare <- var_bar^2 * (2/(df1-1) + 2/(df2-1) - 2*mCorrelationEstimate.standard_errors*sqrt(2/(df1-1))*sqrt(2/(df2-1)))

    pTest.standard_errors<-2*pnorm(q = abs(tVar), mean =  0, sd = sqrt(vare),lower.tail = F)

    pTest.standard_errors.adj<- p.adjust2(pTest.standard_errors, method = "fdr", n = effectiveNumberOfTests, isSymmetric = symmetric)
    rownames(pTest.standard_errors.adj)<-rownames(pTest.standard_errors)
    colnames(pTest.standard_errors.adj)<-colnames(pTest.standard_errors)

  }
  
  # if(!is.null(mStandard_errors1.std) & !is.null(mStandard_errors2.std)){
  #   
  #   tVar <- abs((mStandard_errors1.std^2)-(mStandard_errors2.std^2))
  #   
  #   sde<-sqrt((1^2)+(1^2) - 2*mCorrelationEstimate.standard_errors*1*1)
  # 
  #   pTest.standard_errors<-2*pnorm(q = tVar, mean =  0, sd = sde,lower.tail = F)
  #   
  #   pTest.standard_errors.adj<-matrix(data = p.adjust2(pTest.standard_errors, method = "fdr", n = effectiveNumberOfTests), nrow = nrow(pTest.standard_errors), ncol = ncol(pTest.standard_errors))
  #   rownames(pTest.standard_errors.adj)<-rownames(pTest.standard_errors)
  #   colnames(pTest.standard_errors.adj)<-colnames(pTest.standard_errors)
  #   
  # }
  
  return(
    list(
      pTest.values=pTest.values,
      pTest.values.adj=pTest.values.adj,
      pTest.standard_errors=pTest.standard_errors,
      pTest.standard_errors.adj=pTest.standard_errors.adj
    )
  )
  
}
